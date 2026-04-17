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
use parameters
use flags
use types_basic
use permutation
use lap_jv_sparse
use lap_hungarian
use euclidean
use random
implicit none
private
public assign_atoms
public assign_atoms_pruned

contains

subroutine assign_atoms( atomtypes, costs, atomperm1)
! ----------------------------------------------------------------
! Finds the optimal mapping between points with fixed orientation
! ----------------------------------------------------------------
   type(partition_t), target, intent(in) :: atomtypes
   type(real_matrix), dimension(:), intent(in) :: costs
   integer(ik), dimension(:), intent(out) :: atomperm1
   ! Local variables
   integer(ik), dimension(:), allocatable :: k
   real(rk), dimension(:,:), allocatable :: a
   integer(ik) :: h, i, j, n, m
   real(rk) :: s

   n =  maxval(atomtypes%parts%num_items1)
   m =  maxval(atomtypes%parts%num_items2)
   allocate (k(n))
   allocate (a(n, m))

   ! Optimize atomperm1 for each block
   do h = 1, atomtypes%num_parts
      n = atomtypes%parts(h)%num_items1
      m = atomtypes%parts(h)%num_items2
      a(1:n, 1:m) = costs(h)%a(1:n, 1:m)
      call assndx(1, a, n, m, k, s)
      atomperm1(atomtypes%parts(h)%items1) = atomtypes%parts(h)%items2(k(1:n))
   end do
end subroutine

subroutine assign_atoms_pruned( atomtypes, coords1, coords2, prunes, atomperm1, error_code)
! ----------------------------------------------------------------
! Finds the optimal mapping between points with fixed orientation
! ----------------------------------------------------------------
   type(partition_t), target, intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(bool_matrix), dimension(:), intent(in) :: prunes
   integer(ik), dimension(:), allocatable, intent(out) :: atomperm1
   integer(ik), intent(out) :: error_code
   ! Local variables
   integer(ik), dimension(:), allocatable :: perm1
   integer(ik) :: h, n
   real(rk) :: dist

   error_code = 0
   allocate (perm1(maxval(atomtypes%parts%num_items1)))
   allocate (atomperm1(sum(atomtypes%parts%num_items1)))

   ! Optimize atomperm1 for each block
   do h = 1, atomtypes%num_parts
      n = atomtypes%parts(h)%num_items1
      call solve_lap_pruned(n, atomtypes%parts(h)%items1, atomtypes%parts(h)%items2, &
            coords1, coords2, prunes(h)%a, perm1, dist, error_code)
      if (error_code /= 0) return
      atomperm1(atomtypes%parts(h)%items1) = atomtypes%parts(h)%items2(perm1(1:n))
   end do
end subroutine

subroutine solve_lap_pruned(n, s1, s2, x1, x2, pruned, perm1, dist, error_code)
! Adapted from GMIN: A program for finding global minima
! Copyright (C) 1999-2006 David J. Wales

!   Interface to spjv.f for calculating minimum distance
!   of two atomic configurations with respect to
!   particle permutations.
!   The function permdist determines the distance or weight function,
!
!       Tomas Oppelstrup, Jul 10, 2003
!       tomaso@nada.kth.se
!

!   This is the main routine for minimum distance calculation.
!   Given two coordinate vectors x1,x2 of particles each, return
!   the minimum distance in dist, and the permutation in perm1.
!   perm1 is an integer vector such that
!     x1(i) <--> x2(perm1(i))
!   i.e.
!     sum(i=1,n) permdist(x1(i), x2(perm1(i))) == dist

!   Input
!     n  : System size
!     x1,x2: Coordinate vectors (n particles)

   integer(ik), intent(in) :: n
   integer(ik), intent(in) :: s1(n), s2(n)
   real(rk), intent(in) :: x1(3, *), x2(3, *)
   logical(lk), intent(in) :: pruned(n, n)
   real(rk), parameter :: scale = 1.0e6_rk ! Precision

!   Output
!     perm1: Permutation so that x1(i) <--> x2(perm1(i))
!     dist: Minimum attainable distance
!   We have
   integer(ik), intent(out) :: perm1(n)
   real(rk), intent(out) :: dist
   integer(ik), intent(out) :: error_code

!   Internal variables
!   cc, kk, first:
!     Sparse matrix of distances
!   first(i):
!     Beginning of row i in data,index vectors
!   kk(first(i)..first(i+1)-1):
!     Column indexes of existing elements in row i
!   cc(first(i)..first(i+1)-1):
!     Matrix elements of row i
   integer(ik) :: first(n+1), y(n)
!   integer :: m, i, j, k, l, l2, a, sz, t
   integer(ik) :: i, j, k, sz
   integer(int64) :: u(n), v(n), h
   integer(ik), allocatable :: kk(:)
   integer(int64), allocatable :: cc(:)

   error_code = 0
   sz = n*n - count(pruned)

   allocate (kk(sz))
   allocate (cc(sz))

   first(1) = 1
   do i = 1, n
      first(i+1) = first(i) + n - count(pruned(:, i))
   end do

!  Compute the sparse cost matrix...

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

!   Call bipartite matching routine
   call jovosap(n, sz, cc, kk, first, perm1, y, u, v, h)

   if (h < 0) then
!   If initial guess correct, deduce solution distance
!   which is not done in jovosap
      h = 0
      do i = 1, n
         j = first(i)
30       if (j > sz) then
            write (stderr, '(a)') 'Error: Assignment failed'
            error_code = 3
            return
         end if
         if (kk(j) /= perm1(i)) then
            j = j + 1
            goto 30
         end if
         h = h + cc(j)
      end do
   end if

   if (DEBUG_TESTS) then
      if (.not. is_permutation(perm1)) then
         error stop 'Assignment is not a permutation'
      end if
   end if

   dist = real(h, rk) / scale
end subroutine

end module
