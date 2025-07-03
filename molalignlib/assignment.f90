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

module assignment
use parameters
use globals
use derived_types
use lcrs_tree
use permutation
use lap_jv
use hungarian
use random

implicit none

private
public assign_atoms
public assign_atoms_biased
public assign_atoms_pruned
public assign_atoms_nearest

contains

subroutine assign_atoms_nearest( atomtypes, coords1, coords2, atomperm)
! Find best correspondence between points sets with fixed orientation

   type(partition_t), target, intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, dimension(:), intent(out) :: atomperm
   ! Local variables
   integer :: h, num_items1
   integer, dimension(:), allocatable :: auxperm
   integer, dimension(:), pointer :: items1, items2
   real(rk) :: dist

   allocate (auxperm(maxval(atomtypes%parts%num_items1)))

   ! Fill distance matrix for each block

   do h = 1, atomtypes%num_parts
      num_items1 = atomtypes%parts(h)%num_items1
      items1 => atomtypes%parts(h)%items1
      items2 => atomtypes%parts(h)%items2
      call solve_lap_nearest(num_items1, items1, items2, coords1, coords2, auxperm, dist)
      atomperm(items1) = items2(auxperm(:num_items1))
   end do
end subroutine

subroutine assign_atoms_pruned( atomtypes, coords1, coords2, prunes, atomperm)
! Find best correspondence between points sets with fixed orientation

   type(partition_t), target, intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(bool_matrix), dimension(:), intent(in) :: prunes
   integer, dimension(:), intent(out) :: atomperm
   ! Local variables
   integer :: h, num_items1
   integer, dimension(:), allocatable :: auxperm
   integer, dimension(:), pointer :: items1, items2
   real(rk) :: dist

   allocate (auxperm(maxval(atomtypes%parts%num_items1)))

   ! Optimize atomperm for each block
   do h = 1, atomtypes%num_parts
      num_items1 = atomtypes%parts(h)%num_items1
      items1 => atomtypes%parts(h)%items1
      items2 => atomtypes%parts(h)%items2
      call solve_lap_pruned(num_items1, items1, items2, coords1, coords2, prunes(h)%ee, auxperm, dist)
      atomperm(items1) = items2(auxperm(:num_items1))
   end do
end subroutine

subroutine assign_atoms( atomtypes, coords1, coords2, atomperm, dist)
! Find best correspondence between points sets with fixed orientation

   type(partition_t), intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, dimension(:), intent(out) :: atomperm
   real(rk), intent(out) :: dist
   ! Local variables
   integer :: h
   integer, dimension(:), allocatable :: auxperm

   allocate (auxperm(maxval(atomtypes%parts%num_items1)))

   ! Optimize atomperm for each block
   do h = 1, atomtypes%num_parts
      call solve_lap(atomtypes%parts(h), coords1, coords2, auxperm, dist)
      atomperm(atomtypes%parts(h)%items1) = atomtypes%parts(h)%items2(auxperm(:atomtypes%parts(h)%num_items1))
   end do
end subroutine

subroutine assign_atoms_biased( atomtypes, coords1, coords2, biases, atomperm)
! Find best correspondence between points sets with fixed orientation

   type(partition_t), intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(int_matrix), dimension(:), intent(in) :: biases
   integer, dimension(:), intent(out) :: atomperm
   ! Local variables
   integer :: h
   integer, dimension(:), allocatable :: auxperm
   real(rk) :: dist

   allocate (auxperm(maxval(atomtypes%parts%num_items1)))

   ! Optimize atomperm for each block
   do h = 1, atomtypes%num_parts
      call solve_lap_biased(atomtypes%parts(h), coords1, coords2, biases(h)%ee, auxperm, dist)
      atomperm(atomtypes%parts(h)%items1) = atomtypes%parts(h)%items2(auxperm(:atomtypes%parts(h)%num_items1))
   end do
end subroutine

subroutine solve_lap(part, p, q, perm, dist)
   type(partition_part_t), intent(in) :: part
   real(rk), dimension(:,:), intent(in) :: p, q
   integer, dimension(:), intent(out) :: perm
   real(rk), intent(out) :: dist
   ! Local variables
   integer :: i, j
   real(rk), dimension(:,:), allocatable :: costs

   allocate (costs(part%num_items1, part%num_items2))

   do j = 1, part%num_items2
      do i = 1, part%num_items1
         costs(i, j) = sum((p(:, part%items1(i)) - q(:, part%items2(j)))**2)
      end do
   end do

   call assndx(1, costs, part%num_items1, part%num_items2, perm, dist)
end subroutine

subroutine solve_lap_biased(part, p, q, biases, perm, dist)
   type(partition_part_t), intent(in) :: part
   real(rk), dimension(:,:), intent(in) :: p, q
   integer, dimension(:,:), intent(in) :: biases
   integer, dimension(:), intent(out) :: perm
   real(rk), intent(out) :: dist
   ! Local variables
   real(rk), dimension(:,:), allocatable :: costs
   integer :: i, j, maxbias

   allocate (costs(part%num_items1, part%num_items2))
   maxbias = maxval(biases)

   do j = 1, part%num_items2
      do i = 1, part%num_items1
!         costs(i, j) = maxbias - biases(i, j) + bias_scale*sum((p(:, part%items1(i)) - q(:, part%items2(j)))**2)
         costs(i, j) = maxbias - biases(i, j) + random_standard_real()
      end do
   end do

   call assndx(1, costs, part%num_items1, part%num_items2, perm, dist)
end subroutine

subroutine solve_lap_pruned(n, s1, s2, p, q, prun, perm, dist)
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
!   Given two coordinate vectors p,q of particles each, return
!   the minimum distance in dist, and the permutation in perm.
!   perm is an integer vector such that
!     p(i) <--> q(perm(i))
!   i.e.
!     sum(i=1,n) permdist(p(i), q(perm(i))) == dist

!   Input
!     n  : System size
!     p,q: Coordinate vectors (n particles)

   integer, intent(in) :: n
   integer, intent(in) :: s1(n), s2(n)
   real(rk), intent(in) :: p(3, *), q(3, *)
   logical, intent(in) :: prun(n, n)
   real(rk), parameter :: scale = 1.0e6_rk ! Precision

!   Output
!     perm: Permutation so that p(i) <--> q(perm(i))
!     dist: Minimum attainable distance
!   We have
   integer, intent(out) :: perm(n)
   real(rk), intent(out) :: dist

!   Internal variables
!   cc, kk, first:
!     Sparse matrix of distances
!   first(i):
!     Beginning of row i in data,index vectors
!   kk(first(i)..first(i+1)-1):
!     Column indexes of existing elements in row i
!   cc(first(i)..first(i+1)-1):
!     Matrix elements of row i
   integer :: first(n+1), y(n)
!   integer :: m, i, j, k, l, l2, a, sz, t
   integer :: i, j, k, sz
   integer(int64) :: u(n), v(n), h
   integer, allocatable :: kk(:)
   integer(int64), allocatable :: cc(:)

   sz = n*n - count(prun)

   allocate (kk(sz))
   allocate (cc(sz))

   first(1) = 1
   do i = 1, n
      first(i+1) = first(i) + n - count(prun(:, i))
   end do

!  Compute the sparse cost matrix...

   do i = 1, n
      k = first(i)
      do j = 1, n
         if (.not. prun(j, i)) then
            cc(k) = scale * sum((p(:, s1(i)) - q(:, s2(j)))**2)
            kk(k) = j
            k = k + 1
         end if
      end do
   end do

!   Call bipartite matching routine
   call jovosap(n, sz, cc, kk, first, perm, y, u, v, h)

   if (h < 0) then
!   If initial guess correct, deduce solution distance
!   which is not done in jovosap
      h = 0
      do i = 1, n
         j = first(i)
30       if (j > sz) then
            write (stderr, '(a)') 'Error: Assignment failed'
            stop
         end if
         if (kk(j) /= perm(i)) then
            j = j + 1
            goto 30
         end if
         h = h + cc(j)
      end do
   end if

   if (.not. is_perm(perm)) then
      write (stderr, '(a)') 'Assignment is not a permutation'
      stop
   end if

   dist = real(h, rk) / scale
end subroutine

subroutine solve_lap_nearest(n, s1, s2, p, q, perm, dist)
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
!   Given two coordinate vectors p,q of particles each, return
!   the minimum distance in dist, and the permutation in perm.
!   perm is an integer vector such that
!     p(i) <--> q(perm(i))
!   i.e.
!     sum(i=1,n) permdist(p(i), q(perm(i))) == dist

!   Input
!     n  : System size
!     p,q: Coordinate vectors (n particles)

   integer, intent(in) :: n
   integer, intent(in) :: s1(n), s2(n)
   real(rk), intent(in) :: p(3, *), q(3, *)
   real(rk), parameter :: scale = 1.0e6_rk ! Precision
   integer, parameter :: maxnei = 20 ! Maximum number of closest neighbours

!   Output
!     perm: Permutation so that p(i) <--> q(perm(i))
!     dist: Minimum attainable distance
!   We have
   integer, intent(out) :: perm(n)
   real(rk), intent(out) :: dist

!   Internal variables
!   cc, kk, first:
!     Sparse matrix of distances
!   first(i):
!     Beginning of row i in data,index vectors
!   kk(first(i)..first(i+1)-1):
!     Column indexes of existing elements in row i
!   cc(first(i)..first(i+1)-1):
!     Matrix elements of row i
   integer :: first(n+1), y(n)
   integer :: m, i, j, k, l, l2, a, sz, t
   integer(int64) :: u(n), v(n), d, h
   integer, allocatable :: kk(:)
   integer(int64), allocatable :: cc(:)

!   Distance function
!    real(rk) permdist

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

!  Compute the full matrix...

      do i = 1, n
         k = first(i)
         do j = 1, n
            cc(k) = scale * sum((p(:, s1(i)) - q(:, s2(j)))**2)
            kk(k) = j
            k = k + 1
         end do
      end do

   else

!  We need to store the distances of the maxnei closeest neighbors
!  of each particle. The following builds a heap to keep track of
!  the maxnei closest neighbours seen so far. It might be more
!  efficient to use quick-select instead... (This is definately
!  true in the limit of infinite systems.)

      do i = 1, n
         k = first(i) - 1
         do j = 1, m
            d = scale * sum((p(:, s1(i)) - q(:, s2(j)))**2)
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
            d = scale * sum((p(:, s1(i)) - q(:, s2(j)))**2)
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

!  Create and maintain an ordered list, smallest to largest from kk(m*(i-1)+1:m*i) for atom i.
!  NOTE that there is no symmetry with respect to exchange of I and J!
!  This runs slower than the above heap algorithm.
!
!      cc(1:m*n) = huge(1_int64)
!      do i = 1, n
!         k = first(i) - 1
!         do j = 1, n
!            d = scale * sum((p(:, s1(i)) - q(:, s2(j)))**2)
!            if (d > cc(k+m)) cycle
!            do l = m, 2, -1
!               if (d > cc(k+l-1)) exit
!               cc(k+l) = cc(k+l-1)
!               kk(k+l) = kk(k+l-1)
!            end do
!            cc(k+l) = d
!            kk(k+l) = j
!         end do
!      end do

   end if

!   Call bipartite matching routine
   call jovosap(n, sz, cc, kk, first, perm, y, u, v, h)

   if (h < 0) then
!   If initial guess correct, deduce solution distance
!   which is not done in jovosap
      h = 0
      do i = 1, n
         j = first(i)
30       if (j > sz) then
            write (stderr, '(a)') 'Error: Assignment failed'
            stop
         end if
         if (kk(j) /= perm(i)) then
            j = j + 1
            goto 30
         end if
         h = h + cc(j)
      end do
   end if

   if (.not. is_perm(perm)) then
      error stop 'Assignment is not a permutation'
   end if

   dist = real(h, rk) / scale
end subroutine

end module
