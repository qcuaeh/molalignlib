! MolAlignLib
! Copyright (C) 2022 José M. Vásquez

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
! Topology-unaware atom assignment at fixed orientation (linear sum
! assignment within each atom type), used for atom clusters and as the
! starting point of the isomer search
use parameters
use flags
use common_types
use permutation
use lap_jv_sparse
use lap_hungarian
use euclidean
use random
use error_codes
implicit none
private
public assign_atoms
public assign_atoms_pruned

contains

subroutine assign_atoms( atomtypes, costs, mapping1)
! Assignment minimizing the sum of costs, solved independently for each
! atom type with the Hungarian algorithm (assndx). costs(h) is the cost
! matrix of atom type h.
   type(partition_t), target, intent(in) :: atomtypes
   type(real_matrix), dimension(:), intent(in) :: costs
   integer(ik), dimension(:), intent(out) :: mapping1
   ! Local variables
   integer(ik), dimension(:), allocatable :: k
   real(rk), dimension(:,:), allocatable :: a
   integer(ik) :: h, i, j, n, m
   real(rk) :: s

   n =  maxval(atomtypes%parts%n_items1)
   m =  maxval(atomtypes%parts%n_items2)
   allocate (k(n))
   allocate (a(n, m))

   ! Every atom belongs to one block, so mapping1 is fully assigned
   do h = 1, atomtypes%n_parts
      n = atomtypes%parts(h)%n_items1
      m = atomtypes%parts(h)%n_items2
      a(1:n, 1:m) = costs(h)%a(1:n, 1:m)
      call assndx(1, a, n, m, k, s)
      mapping1(atomtypes%parts(h)%items1) = atomtypes%parts(h)%items2(k(1:n))
   end do
end subroutine

subroutine assign_atoms_pruned( atomtypes, coords1, coords2, prunes, mapping1, error_code)
! Assignment minimizing the sum of squared distances, solved independently
! for each atom type with the sparse Jonker-Volgenant algorithm; pairs
! marked in prunes(h) are excluded
   type(partition_t), target, intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(bool_matrix), dimension(:), intent(in) :: prunes
   ! Allocated by the caller with size(atomtypes%itemdir1) elements
   integer(ik), dimension(:), intent(out) :: mapping1
   integer(ik), intent(out) :: error_code
   ! Local variables
   integer(ik), dimension(:), allocatable :: submap1
   integer(ik) :: h, n
   real(rk) :: dist

   error_code = MOLALIGN_SUCCESS
   allocate (submap1(maxval(atomtypes%parts%n_items1)))

   ! Every atom belongs to one atom type, so mapping1 is fully assigned
   do h = 1, atomtypes%n_parts
      n = atomtypes%parts(h)%n_items1
      call solve_lap_pruned(n, atomtypes%parts(h)%items1, atomtypes%parts(h)%items2, &
            coords1, coords2, prunes(h)%a, submap1, dist, error_code)
      if (error_code /= MOLALIGN_SUCCESS) return
      mapping1(atomtypes%parts(h)%items1) = atomtypes%parts(h)%items2(submap1(1:n))
   end do
end subroutine

subroutine solve_lap_pruned(n, s1, s2, x1, x2, pruned, submap1, dist, error_code)
! Minimum sum of squared distances between atoms x1(:,s1(i)) and
! x2(:,s2(submap1(i))), i = 1..n, over the permutations submap1 that avoid
! the pruned pairs (pruned(j,i) excludes pairing s1(i) with s2(j)). The
! sparse cost matrix is solved with jovosap in fixed point (scale).
! Adapted from GMIN (Copyright (C) 1999-2006 David J. Wales), interface by
! Tomas Oppelstrup.
   integer(ik), intent(in) :: n
   integer(ik), intent(in) :: s1(n), s2(n)
   real(rk), intent(in) :: x1(3, *), x2(3, *)
   logical(lk), intent(in) :: pruned(n, n)
   real(rk), parameter :: scale = 1.0e6_rk ! Fixed-point precision
   integer(ik), intent(out) :: submap1(n)
   real(rk), intent(out) :: dist
   integer(ik), intent(out) :: error_code

   ! Sparse cost matrix: row i has columns kk(first(i):first(i+1)-1) with
   ! costs cc(first(i):first(i+1)-1)
   integer(ik) :: first(n+1), y(n)
   integer(ik) :: i, j, k, sz
   integer(int64) :: u(n), v(n), h
   integer(ik), allocatable :: kk(:)
   integer(int64), allocatable :: cc(:)
   logical(lk) :: found, col_used(n)

   error_code = MOLALIGN_SUCCESS

   ! A fully pruned row or column leaves no feasible assignment. Catch it
   ! before calling jovosap, which does not report infeasibility and whose
   ! result on such a matrix is undefined.
   do i = 1, n
      if (all(pruned(:, i)) .or. all(pruned(i, :))) then
         error_code = MOLALIGN_ERROR_PRUNED_ASSIGNMENT_FAILED
         return
      end if
   end do

   sz = n*n - count(pruned)

   allocate (kk(sz))
   allocate (cc(sz))

   first(1) = 1
   do i = 1, n
      first(i+1) = first(i) + n - count(pruned(:, i))
   end do

   do i = 1, n
      k = first(i)
      do j = 1, n
         if (.not. pruned(j, i)) then
            cc(k) = nint(sum((x1(:, s1(i)) - x2(:, s2(j)))**2)*scale, int64)
            kk(k) = j
            k = k + 1
         end if
      end do
   end do

   call jovosap(n, sz, cc, kk, first, submap1, y, u, v, h)

   ! Validate the assignment: every pair must be unpruned and no column used
   ! twice, otherwise no valid assignment exists under the pruning. The cost
   ! is recomputed here because jovosap leaves h negative when its initial
   ! guess is already optimal.
   col_used = .false.
   h = 0
   do i = 1, n
      if (submap1(i) < 1 .or. submap1(i) > n) then
         error_code = MOLALIGN_ERROR_PRUNED_ASSIGNMENT_FAILED
         return
      end if
      if (col_used(submap1(i))) then
         error_code = MOLALIGN_ERROR_PRUNED_ASSIGNMENT_FAILED
         return
      end if
      col_used(submap1(i)) = .true.
      found = .false.
      do k = first(i), first(i+1) - 1
         if (kk(k) == submap1(i)) then
            h = h + cc(k)
            found = .true.
            exit
         end if
      end do
      if (.not. found) then
         ! Assignment uses a pruned pair: pruning tolerance might be too tight
         error_code = MOLALIGN_ERROR_PRUNED_ASSIGNMENT_FAILED
         return
      end if
   end do

   dist = real(h, rk) / scale
end subroutine

end module
