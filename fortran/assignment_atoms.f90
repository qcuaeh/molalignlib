! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module assignment_atoms
use parameters
use types_basic
use random
use euclidean
use permutation
use lap_hungarian
use lap_jv_sparse
use options
implicit none
private
public assign_atoms
public assign_atoms_pruned
public assign_atoms_nearest

contains

subroutine assign_atoms(atomtypes, costs, atomperm1)
!------------------------------------------------------------------------
! Finds the optimal mapping between points with fixed orientation.
! Uses assndx (Hungarian algorithm) per partition block.
! Assumes num_items1 == num_items2 for all parts.
!------------------------------------------------------------------------
   type(partition_t), target, intent(in) :: atomtypes
   type(real_matrix), dimension(:), intent(in) :: costs
   integer(ik), dimension(:), allocatable, intent(inout) :: atomperm1

   type(partition_part_t), pointer :: part
   integer(ik), dimension(:), allocatable :: partperm
   real(rk), allocatable :: a(:,:)
   real(rk) :: lapcost
   integer(ik) :: h, m

   m = maxval(atomtypes%parts%num_items1)
   allocate (partperm(m), a(m, m))

   do h = 1, atomtypes%num_parts
      part => atomtypes%parts(h)
      m = part%num_items1
      ! assndx modifies its cost matrix in-place, so copy costs into a
      ! working array. assndx uses a(row, col) = a(item1, item2) with
      ! MODE=1, which minimises Sum_i a(i, k(i))
      a(1:m, 1:m) = costs(h)%a(1:m, 1:m)
      call assndx(1, a(1:m, 1:m), m, m, partperm(1:m), lapcost)
      ! partperm(i) is the item2 index assigned to item1 i.
      atomperm1(part%items1) = part%items2(partperm(1:m))
   end do
end subroutine

subroutine assign_atoms_pruned(atomtypes, coords1, coords2, prunes, atomperm1)
!------------------------------------------------------------------------
! Finds the optimal mapping between points with fixed orientation.
! Uses jovosap (Jonker-Volgenant sparse algorithm) with a pruned cost matrix.
!------------------------------------------------------------------------
   type(partition_t), target, intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(bool_matrix), dimension(:), intent(in) :: prunes
   integer(ik), dimension(:), allocatable, intent(out) :: atomperm1

   type(partition_part_t), pointer :: part
   integer(ik), dimension(:), allocatable :: partperm
   real(rk) :: dist
   integer(ik) :: h

   allocate (partperm(maxval(atomtypes%parts%num_items1)))

   ! Initialise atomperm1 as the identity permutation.
   call init_identity_permutation(size(coords1, 2), atomperm1)

   do h = 1, atomtypes%num_parts
      part => atomtypes%parts(h)
      call solve_lap_pruned(part%num_items1, part%items1, part%items2, &
                            coords1, coords2, prunes(h)%a, partperm, dist)
      atomperm1(part%items1) = part%items2(partperm(1:part%num_items1))
   end do
end subroutine

subroutine assign_atoms_nearest(atomtypes, coords1, coords2, atomperm1)
!------------------------------------------------------------------------
! Finds the optimal mapping between points with fixed orientation.
! Uses jovosap (Jonker-Volgenant sparse algorithm) restricted to nearest neighbours.
!------------------------------------------------------------------------
   type(partition_t), target, intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer(ik), dimension(:), allocatable, intent(out) :: atomperm1

   type(partition_part_t), pointer :: part
   integer(ik), dimension(:), allocatable :: partperm
   real(rk) :: dist
   integer(ik) :: h

   allocate (partperm(maxval(atomtypes%parts%num_items1)))

   ! Initialise atomperm1 as the identity permutation.
   call init_identity_permutation(size(coords1, 2), atomperm1)

   do h = 1, atomtypes%num_parts
      part => atomtypes%parts(h)
      call solve_lap_nearest(part%num_items1, part%items1, part%items2, &
                             coords1, coords2, partperm, dist)
      atomperm1(part%items1) = part%items2(partperm(1:part%num_items1))
   end do
end subroutine

subroutine solve_lap_pruned(n, s1, s2, x1, x2, pruned, partperm, dist)
!--------------------------------------------------------------------
! Interface to jovosap for minimum-distance atom mapping with a
! pruned (sparse) cost matrix.
!
! partperm(i) = j  means  x1(:, s1(i)) <--> x2(:, s2(j))
! dist = sum_i  ||x1(:,s1(i)) - x2(:,s2(partperm(i)))||^2
!--------------------------------------------------------------------
   integer(ik), intent(in) :: n
   integer(ik), intent(in) :: s1(n), s2(n)
   real(rk),    intent(in) :: x1(3, *), x2(3, *)
   logical(lk), intent(in) :: pruned(n, n)

   integer(ik), intent(out) :: partperm(n)
   real(rk),    intent(out) :: dist

   ! Sparse matrix storage (CSR-style)
   integer(ik), allocatable :: kk(:), first(:), y(:)
   integer(int64), allocatable :: cc(:), u(:), v(:)
   integer(ik) :: i, j, k, sz
   integer(int64) :: h
   real(rk), parameter :: scale = 1.0e6_rk

   sz = n*n - count(pruned)

   allocate (kk(sz), cc(sz), first(n+1), y(n), u(n), v(n))

   ! Build row-pointer array.
   first(1) = 1
   do i = 1, n
      first(i+1) = first(i) + n - count(pruned(:, i))
   end do

   ! Fill sparse cost matrix (scaled squared distances).
   do i = 1, n
      k = first(i)
      do j = 1, n
         if (.not. pruned(j, i)) then
            cc(k) = nint(sum((x1(:, s1(i)) - x2(:, s2(j)))**2) * scale, int64)
            kk(k) = j
            k = k + 1
         end if
      end do
   end do

   call jovosap(n, sz, cc, kk, first, partperm, y, u, v, h)

   if (h < 0) then
      ! jovosap returns h=-1 when the initial guess was already optimal;
      ! recompute the objective from the stored cost entries.
      h = 0_int64
      do i = 1, n
         j = first(i)
         do while (kk(j) /= partperm(i))
            j = j + 1
         end do
         h = h + cc(j)
      end do
   end if

   if (DEBUG_TESTS) then
      if (.not. is_permutation(partperm)) then
         error stop 'Assignment is not a permutation'
      end if
   end if

   dist = real(h, rk) / scale
end subroutine

subroutine solve_lap_nearest(n, s1, s2, x1, x2, partperm, dist)
!--------------------------------------------------------------------
! Interface to jovosap restricted to the maxnei nearest neighbours
! of each atom (sparse approximation to the full assignment problem).
!
! Adapted from GMIN: A program for finding global minima
! Copyright (C) 1999-2006 David J. Wales
! Original interface by Tomas Oppelstrup, Jul 10, 2003.
!
! partperm(i) = j  means  x1(:, s1(i)) <--> x2(:, s2(j))
! dist = sum_i  ||x1(:,s1(i)) - x2(:,s2(partperm(i)))||^2
!--------------------------------------------------------------------
   integer(ik), intent(in) :: n
   integer(ik), intent(in) :: s1(n), s2(n)
   real(rk),    intent(in) :: x1(3, *), x2(3, *)
   integer(ik), parameter  :: maxnei = 20

   integer(ik), intent(out) :: partperm(n)
   real(rk),    intent(out) :: dist

   integer(ik), allocatable :: kk(:), first(:), y(:)
   integer(int64), allocatable :: cc(:), u(:), v(:)
   integer(ik) :: m, i, j, k, l, l2, a, sz, t
   integer(int64) :: d_int, hswap, h
   real(rk), parameter :: scale = 1.0e6_rk

   m = min(n, maxnei)
   sz = m * n

   allocate (kk(sz), cc(sz), first(n+1), y(n), u(n), v(n))

   first(1) = 1
   do i = 1, n
      first(i+1) = first(i) + m
   end do

   if (m == n) then
      ! Full matrix case.
      do i = 1, n
         k = first(i)
         do j = 1, n
            cc(k) = nint(sum((x1(:, s1(i)) - x2(:, s2(j)))**2) * scale, int64)
            kk(k) = j
            k = k + 1
         end do
      end do

   else
      ! Sparse case: keep only the m nearest neighbours per row using
      ! a max-heap of size m so each row costs O(n log m).
      do i = 1, n
         k = first(i) - 1

         ! Seed the heap with the first m neighbours.
         do j = 1, m
            cc(k+j) = nint(sum((x1(:, s1(i)) - x2(:, s2(j)))**2) * scale, int64)
            kk(k+j) = j
            ! Sift up.
            l = j
            do while (l > 1)
               l2 = l / 2
               if (cc(k+l2) < cc(k+l)) then
                  hswap    = cc(k+l2); cc(k+l2) = cc(k+l); cc(k+l) = hswap
                  t        = kk(k+l2); kk(k+l2) = kk(k+l); kk(k+l) = t
                  l = l2
               else
                  exit
               end if
            end do
         end do

         ! Process remaining neighbours; replace heap root if closer.
         do j = m+1, n
            d_int = nint(sum((x1(:, s1(i)) - x2(:, s2(j)))**2) * scale, int64)
            if (d_int < cc(k+1)) then
               cc(k+1) = d_int
               kk(k+1) = j
               ! Sift down.
               l = 1
               do
                  l2 = 2 * l
                  if (l2 + 1 <= m) then
                     a = merge(k+l2+1, k+l2, cc(k+l2+1) > cc(k+l2))
                  else if (l2 <= m) then
                     a = k + l2
                  else
                     exit
                  end if
                  if (cc(a) > cc(k+l)) then
                     hswap  = cc(a); cc(a) = cc(k+l); cc(k+l) = hswap
                     t      = kk(a); kk(a) = kk(k+l); kk(k+l) = t
                     l      = a - k
                  else
                     exit
                  end if
               end do
            end if
         end do
      end do
   end if

   call jovosap(n, sz, cc, kk, first, partperm, y, u, v, h)

   if (h < 0) then
      h = 0_int64
      do i = 1, n
         j = first(i)
         do while (kk(j) /= partperm(i))
            j = j + 1
         end do
         h = h + cc(j)
      end do
   end if

   if (DEBUG_TESTS) then
      if (.not. is_permutation(partperm)) then
         error stop 'Assignment is not a permutation'
      end if
   end if

   dist = real(h, rk) / scale
end subroutine

end module
