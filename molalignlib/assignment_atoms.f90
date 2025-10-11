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
use options
use derived_types
use random
use euclidean
use permutation
use lap_jvc_dense
use lap_jvc_sparse
!use lap_hungarian
implicit none
private
public assign_atoms_biased
public assign_atoms_pruned
public assign_atoms_nearest

contains

subroutine assign_atoms_biased( atomtypes, biases, coords1, coords2, atomperm)
!------------------------------------------------------------------------
! Finds the optimal mapping between points with fixed orientation
! Uses JVC algorithm (assumes num_items1 == num_items2 for all parts)
!------------------------------------------------------------------------
   type(partition_t), target, intent(in) :: atomtypes
   type(int_matrix), dimension(:), intent(in) :: biases
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, dimension(:), allocatable, intent(out) :: atomperm
   ! Local variables
   type(partition_part_t), pointer :: part
   integer, dimension(:), allocatable :: perm
   real(rk), dimension(:,:), allocatable :: costs
   integer :: maxnum_items
   real(rk) :: lapcost
   integer :: h

   ! Since num_items1 == num_items2, we only need one size
   maxnum_items = maxval(atomtypes%parts%num_items1)

   allocate (perm(maxnum_items))
   allocate (atomperm(sum(atomtypes%parts%num_items1)))
   allocate (costs(maxnum_items, maxnum_items))

   ! Optimize atomperm for each block
   do h = 1, atomtypes%num_parts
      part => atomtypes%parts(h)

      ! Build cost matrix: distances + biases
      costs(1:part%num_items1, 1:part%num_items1) = &
         distance_matrix(part, coords1, coords2) + biases(h)%ee

      ! Solve assignment problem using JVC algorithm
      call jvc_dense(costs, part%num_items1, perm, lapcost)

      ! Map the solution back to original atom indices
      atomperm(part%items1) = part%items2(perm(1:part%num_items1))
   end do
end subroutine

subroutine assign_atoms_pruned( atomtypes, coords1, coords2, prunes, atomperm)
!------------------------------------------------------------------------
! Finds the optimal mapping between points with fixed orientation
!------------------------------------------------------------------------
   type(partition_t), target, intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(bool_matrix), dimension(:), intent(in) :: prunes
   integer, dimension(:), allocatable, intent(out) :: atomperm
   ! Local variables
   type(partition_part_t), pointer :: part
   integer, dimension(:), allocatable :: perm
   real(rk) :: dist
   integer :: h

   allocate (perm(maxval(atomtypes%parts%num_items1)))
   allocate (atomperm(sum(atomtypes%parts%num_items1)))

   ! Optimize atomperm for each block
   do h = 1, atomtypes%num_parts
      part => atomtypes%parts(h)
      call solve_lap_pruned(part%num_items1, part%items1, part%items2, coords1, coords2, prunes(h)%ee, perm, dist)
      atomperm(part%items1) = part%items2(perm(1:part%num_items1))
   end do
end subroutine

subroutine solve_lap_pruned(n, s1, s2, x1, x2, prun, perm, dist)
!--------------------------------------------------------------------
! Interface to JVC sparse algorithm for calculating minimum distance
! of two atomic configurations with respect to
! particle permutations.
!
! This is the main routine for minimum distance calculation.
! Given two coordinate vectors x1,x2 of particles each, return
! the minimum distance in dist, and the permutation in perm.
! perm is an integer vector such that
!   x1(i) <--> x2(perm(i))
! i.e.
!   sum(i=1,n) distance(x1(i), x2(perm(i))) == dist
!
! Input
!   n  : System size
!   x1,x2: Coordinate vectors (n particles)
!--------------------------------------------------------------------

!  Input
!   n: System size
!   x1,x2: Coordinate vectors (n particles)
   integer, intent(in) :: n
   integer, intent(in) :: s1(n), s2(n)
   real(rk), intent(in) :: x1(3, *), x2(3, *)
   logical, intent(in) :: prun(n, n)

!  Output
!   perm: Permutation so that x1(i) <--> x2(perm(i))
!   dist: Minimum attainable distance
   integer, intent(out) :: perm(n)
   real(rk), intent(out) :: dist

!  Local variables
!   cc, kk, first:
!     Sparse matrix of distances
!   first(i):
!     Beginning of row i in data,index vectors
!   kk(first(i)..first(i+1)-1):
!     Column indexes of existing elements in row i
!   cc(first(i)..first(i+1)-1):
!     Matrix elements of row i
   integer :: first(n+1)
   integer :: i, j, k, sz, ierr
   real(rk), allocatable :: cc(:)
   integer, allocatable :: kk(:)

   ! Calculate size of sparse matrix (excluding pruned elements)
   sz = n*n - count(prun)

   allocate (kk(sz))
   allocate (cc(sz))

   ! Build first array (row pointers)
   first(1) = 1
   do i = 1, n
      first(i+1) = first(i) + n - count(prun(:, i))
   end do

   ! Build sparse cost matrix (squared distances, no scaling needed)
   do i = 1, n
      k = first(i)
      do j = 1, n
         if (.not. prun(j, i)) then
            cc(k) = sum((x1(:, s1(i)) - x2(:, s2(j)))**2)
            kk(k) = j
            k = k + 1
         end if
      end do
   end do

   ! Call JVC sparse bipartite matching routine
   call jvc_sparse(n, sz, cc, kk, first, perm, dist, ierr)

   if (ierr /= 0) then
      write (stderr, '(a)') 'Error: Assignment failed'
      stop
   end if

   if (DEBUGGING) then
      if (.not. is_permutation(perm)) then
         error stop 'Assignment is not a permutation'
      end if
   end if
end subroutine

subroutine assign_atoms_nearest( atomtypes, coords1, coords2, atomperm)
!------------------------------------------------------------------------
! Finds the optimal mapping between points with fixed orientation
!------------------------------------------------------------------------
   type(partition_t), target, intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, dimension(:), allocatable, intent(out) :: atomperm
   ! Local variables
   type(partition_part_t), pointer :: part
   integer, dimension(:), allocatable :: perm
   real(rk) :: dist
   integer :: h

   allocate (perm(maxval(atomtypes%parts%num_items1)))
   allocate (atomperm(sum(atomtypes%parts%num_items1)))

   ! Fill distance matrix for each block

   do h = 1, atomtypes%num_parts
      part => atomtypes%parts(h)
      call solve_lap_nearest(part%num_items1, part%items1, part%items2, coords1, coords2, perm, dist)
      atomperm(part%items1) = part%items2(perm(1:part%num_items1))
   end do
end subroutine

subroutine solve_lap_nearest(n, s1, s2, x1, x2, perm, dist)
!--------------------------------------------------------------------
! Adapted from GMIN: A program for finding global minima
! Copyright (C) 1999-2006 David J. Wales
!
! Interface to JVC sparse algorithm for calculating minimum distance
! of two atomic configurations with respect to
! particle permutations.
! The function permdist determines the distance or weight function,
!
!     Tomas Oppelstrup, Jul 10, 2003
!     tomaso@nada.kth.se
!
! This is the main routine for minimum distance calculation.
! Given two coordinate vectors x1,x2 of particles each, return
! the minimum distance in dist, and the permutation in perm.
! perm is an integer vector such that
!   x1(i) <--> x2(perm(i))
! i.e.
!   sum(i=1,n) distance(x1(i), x2(perm(i))) == dist
!--------------------------------------------------------------------

!  Input
!   n: System size
!   x1,x2: Coordinate vectors (n particles)
   integer, intent(in) :: n
   integer, intent(in) :: s1(n), s2(n)
   real(rk), intent(in) :: x1(3, *), x2(3, *)
   integer, parameter :: maxnei = 20 ! Maximum number of closest neighbours

!  Output
!   perm: Permutation so that x1(i) <--> x2(perm(i))
!   dist: Minimum attainable distance
   integer, intent(out) :: perm(n)
   real(rk), intent(out) :: dist

!  Local variables
!   cc, kk, first:
!     Sparse matrix of distances
!   first(i):
!     Beginning of row i in data,index vectors
!   kk(first(i)..first(i+1)-1):
!     Column indexes of existing elements in row i
!   cc(first(i)..first(i+1)-1):
!     Matrix elements of row i
   integer :: first(n+1)
   integer :: m, i, j, k, l, l2, a, sz, t, ierr
   real(rk) :: d, h
   real(rk), allocatable :: cc(:)
   integer, allocatable :: kk(:)

   if (n <= maxnei) then
      m = n
   else
      m = maxnei
   end if

   sz = m*n

   allocate (kk(sz))
   allocate (cc(sz))

   first(1) = 1
   do i = 1, n
      first(i+1) = first(i) + m
   end do

   if (m == n) then

!  Compute the full matrix (no scaling needed)...

      do i = 1, n
         k = first(i)
         do j = 1, n
            cc(k) = sum((x1(:, s1(i)) - x2(:, s2(j)))**2)
            kk(k) = j
            k = k + 1
         end do
      end do

   else

!  We need to store the distances of the maxnei closest neighbors
!  of each particle. The following builds a heap to keep track of
!  the maxnei closest neighbours seen so far. It might be more
!  efficient to use quick-select instead... (This is definitely
!  true in the limit of infinite systems.)

      do i = 1, n
         k = first(i) - 1
         do j = 1, m
            d = sum((x1(:, s1(i)) - x2(:, s2(j)))**2)
            cc(k+j) = d
            kk(k+j) = j
            l = j
10             if (l <= 1) goto 11
            l2 = l/2
            if (cc(k+l2) < cc(k+l)) then
               h = cc(k+l2)
               cc(k+l2) = cc(k+l)
               cc(k+l) = h
               t = kk(k+l2)
               kk(k+l2) = kk(k+l)
               kk(k+l) = t
               l = l2
               goto 10
            end if
11       end do
         do j = m+1, n
            d = sum((x1(:, s1(i)) - x2(:, s2(j)))**2)
            if (d < cc(k+1)) then
               cc(k+1) = d
               kk(k+1) = j
               l = 1
20                l2 = 2*l
               if (l2+1 > m) goto 21
               if (cc(k+l2+1) > cc(k+l2)) then
                  a = k+l2+1
               else
                  a = k+l2
               end if
               if (cc(a) > cc(k+l)) then
                  h = cc(a)
                  cc(a) = cc(k+l)
                  cc(k+l) = h
                  t = kk(a)
                  kk(a) = kk(k+l)
                  kk(k+l) = t
                  l = a-k
                  goto 20
               end if
21             if (l2 <= m) then ! split if statements to avoid a segmentation fault
                  if (cc(k+l2) > cc(k+l)) then
                     h = cc(k+l2)
                     cc(k+l2) = cc(k+l)
                     cc(k+l) = h
                     t = kk(k+l2)
                     kk(k+l2) = kk(k+l)
                     kk(k+l) = t
                  end if
               end if
            end if
         end do
!      PRINT '(A,I6,A)','atom ',i,' nearest neighbours and distances:'
!      PRINT '(20I6)',kk(m*(i-1)+1:m*i)
!      PRINT '(12I15)',cc(m*(i-1)+1:m*i)
      end do

   end if

!   Call JVC sparse bipartite matching routine
   call jvc_sparse(n, sz, cc, kk, first, perm, dist, ierr)

   if (ierr /= 0) then
      write (stderr, '(a)') 'Error: Assignment failed'
      stop
   end if

   if (DEBUGGING) then
      if (.not. is_permutation(perm)) then
         error stop 'Assignment is not a permutation'
      end if
   end if
end subroutine

end module
