module lap_sparse
use iso_fortran_env, only: int64
use kinds
use stdio
use permutation
implicit none
private
public minperm_nearest
public minperm_pruned

!   Parameters
!     scale : Precision
!     maxnei: Maximum number of closest neighbours
real(rk), parameter :: scale = 1.0e6_rk
integer, parameter :: maxnei = 20

contains

subroutine minperm_nearest(n, p, q, mask, perm, dist)

! Adapted from GMIN: A program for finding global minima
! Copyright (C) 1999-2006 David J. Wales

!   Interface to spjv.f for calculating minimum distance
!   of two atomic configurations with respect to
!   particle permutations.
!   The function permdist determines the distance or weight function,
!   minperm is the main routine.
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

   integer :: n
   real(rk) p(3, n), q(3, n)
   logical mask(n, n)

!   Output
!     perm: Permutation so that p(i) <--> q(perm(i))
!     dist: Minimum attainable distance
!   We have
   integer :: perm(n)
   real(rk) dist
   
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
            cc(k) = scale * sum((p(:, i) - q(:, j))**2)
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
            d = scale * sum((p(:, i) - q(:, j))**2)
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
            d = scale * sum((p(:, i) - q(:, j))**2)
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
!            d = scale * sum((p(:, i) - q(:, j))**2)
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

   if (.not. is_permutation(perm)) then
      write (stderr, '(a)') 'Error: Assignment returned by jovosap is not a permutation'
      stop
   end if

   dist = real(h, rk) / scale

end subroutine

subroutine minperm_pruned(n, p, q, mask, perm, dist)

! Adapted from GMIN: A program for finding global minima
! Copyright (C) 1999-2006 David J. Wales

!   Interface to spjv.f for calculating minimum distance
!   of two atomic configurations with respect to
!   particle permutations.
!   The function permdist determines the distance or weight function,
!   minperm is the main routine.
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

   integer :: n
   real(rk) p(3, n), q(3, n)
   logical mask(n, n)

!   Output
!     perm: Permutation so that p(i) <--> q(perm(i))
!     dist: Minimum attainable distance
!   We have
   integer :: perm(n)
   real(rk) dist
   
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
   integer(int64) :: u(n), v(n), d, h
   integer, allocatable :: kk(:)
   integer(int64), allocatable :: cc(:)

   sz = count(mask)

   allocate (kk(sz))
   allocate (cc(sz))

   first(1) = 1
   do i = 1, n
      first(i+1) = first(i) + count(mask(:, i))
   end do

!  Compute the maskd cost matrix...

   do i = 1, n
      k = first(i)
      do j = 1, n
         if (mask(j, i)) then
            cc(k) = scale * sum((p(:, i) - q(:, j))**2)
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

   if (.not. is_permutation(perm)) then
      write (stderr, '(a)') 'Error: Assignment returned by jovosap is not a permutation'
      stop
   end if

   dist = real(h, rk) / scale

end subroutine
   
subroutine jovosap(n,sz,cc,kk,first,x,y,u,v,h)

! Adapted from GMIN: A program for finding global minima
! Copyright (C) 1999-2006 David J. Wales

! This subroutine performs weighted bipartite matching for
! for a sparse non-negative integer weight matrix.
! The original source is
!     http://www.magiclogic.com/assignment.html
! A publication reference can be found on the above homepage and
! in a comment below

! THIS SUBROUTINE SOLVES THE SPARSE LINEAR ASSIGNMENT PROBLEM
! ACCORDING 
!
! "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear   
!  Assignment Problems," Computing 38, 325-340, 1987
! 
! by
! 
! R. Jonker and A. Volgenant, University of Amsterdam.

! INPUT PARAMETERS :
! N = NUMBER OF ROWS AND COLUMNS
! C = WEIGHT MATRIX

! OUTPUT PARAMETERS
! X = COL ASSIGNED TO ROW
! Y = ROW ASSIGNED TO COL
! U = DUAL ROW VARIABLE
! V = DUAL COLUMN VARIABLE
! H = VALUE OF OPTIMAL SOLUTION

! INITIALIZATION

!   Next line added by tomaso@nada.kth.se, to enable detection
!   of solutions being equivalent to the initial guess

! If Y(:) is initialised to zero then we see segmentation faults if 
! a Y element is unset, etc.

   integer :: n,sz
   integer :: kk(sz),first(n+1),x(n),y(n)
   integer :: i,i0,j,j0,j1,k,l,l0,t,t0,td,cnt
   integer :: lab(n),free(n),todo(n)
   integer(int64) :: cc(sz),u(n),v(n),d(n)
   integer(int64) :: h,v0,vj,dj,min
   logical :: ok(n)

   integer(int64), parameter :: bigint = 1000000000000_int64

   y(1:n) = 0
   x(1:n) = 0
   todo(1:n)=0
   h = -1
   do j=1,n
      v(j)=bigint
   end do
   do i=1,n
      x(i)=0
      do t=first(i),first(i+1)-1
         j=kk(t)
         if (cc(t) < v(j)) then
           v(j)=cc(t)
           y(j)=i
         end if
      end do
   end do
   do j=1,n
      j0=n-j+1
      i=y(j0)
      if (i == 0) then
!         print '(a,i6,a)','minperm> warning b - matching failed'
         return
      end if
      if (x(i) /= 0) then
        x(i)=-abs(x(i))
        y(j0)=0
      else
        x(i)=j0
      end if
   end do
   l=0
   do i=1,n
      if (x(i) == 0) then
        l=l+1
        free(l)=i
        cycle
      end if
      if (x(i) < 0) then
        x(i)=-x(i)
      else
        j1=x(i)
        min=bigint
        do t=first(i),first(i+1)-1
           j=kk(t)
           if (j == j1) cycle
           if (cc(t)-v(j) < min) min=cc(t)-v(j)
        end do
        v(j1)=v(j1)-min
      end if
   end do
! improve the initial solution
   cnt=0
   if (l == 0) goto 1000
41 l0=l
   k=1
   l=0
50 i=free(k)
   k=k+1
   v0=bigint
   vj=bigint
   do t=first(i),first(i+1)-1
      j=kk(t)
      h=cc(t)-v(j)
      if (h < vj) then
        if (h >= v0) then
          vj=h
          j1=j
        else
          vj=v0
          v0=h
          j1=j0
          j0=j
        end if
      end if
   end do
   i0=y(j0)
   if (v0 < vj) then
     v(j0)=v(j0)-vj+v0
   else
     if (i0 == 0) goto 43
     j0=j1
     i0=y(j1)
   end if
   if (i0 == 0) goto 43
   if (v0 < vj) then
     k=k-1
     free(k)=i0
   else
     l=l+1
     free(l)=i0
   end if
43 x(i)=j0
   y(j0)=i
   if (k <= l0) goto 50
   cnt=cnt+1
   if ((l > 0).and.(cnt < 2)) goto 41
! augmentation part
   l0=l
   do l=1,l0
      do j=1,n
         ok(j)=.false.
         d(j)=bigint
      end do
      min=bigint
      i0=free(l)
      td=n
      do t=first(i0),first(i0+1)-1
         j=kk(t)
         dj=cc(t)-v(j)
         d(j)=dj
         lab(j)=i0
         if (dj <= min) then
           if (dj < min) then
             min=dj
             k=1
             todo(1)=j
           else
             k=k+1
             todo(k)=j
           end if
         end if
      end do
      do j0=1,k
         j=todo(j0)
         if (j == 0) then
!            print '(a,i6,a)','minperm> warning c - matching failed'
            return
         end if
         if (y(j) == 0) goto 80
         ok(j)=.true.
      end do
! repeat until a free row has been found
60    if (k == 0) then
!         print '(a,i6,a)','minperm> warning d - matching failed'
         return
      end if
      j0=todo(k)
      k=k-1
      i=y(j0)
      todo(td)=j0
      td=td-1
      t0=first(i)
      t=t0-1
61    t=t+1
      if (kk(t) /= j0) goto 61
      h=cc(t)-v(j0)-min
      do t=t0,first(i+1)-1
         j=kk(t)
         if (.not. ok(j)) then
           vj=cc(t)-h-v(j)
           if (vj < d(j)) then
             d(j)=vj
             lab(j)=i
             if (vj == min) then
               if (y(j) == 0) goto 70
               k=k+1
               todo(k)=j
               ok(j)=.true.
             end if
           end if
         end if
      end do
      if (k /= 0) goto 60
      min=bigint-1
      do j=1,n
         if (d(j) <= min) then
           if (.not. ok(j)) then
             if (d(j) < min) then
               min=d(j)
               k=1
               todo(1)=j
             else
               k=k+1
               todo(k)=j
             end if
           end if
         end if
      end do
      do j0=1,k
         j=todo(j0)
         if (y(j) == 0) goto 70
         ok(j)=.true.
      end do
      goto 60
70    if (min == 0) goto 80
      do k=td+1,n
         j0=todo(k)
         v(j0)=v(j0)+d(j0)-min
      end do
80    i=lab(j)
      y(j)=i
      k=j
      j=x(i)
      x(i)=k
      if (i0 /= i) goto 80
   end do
   h=0
   do i=1,n
      j=x(i)
      t=first(i)-1
101   t=t+1
      if (t > sz) then
         print '(a,i6,a)','minperm> warning d - atom ',i,' not matched - maximum number of neighbours too small?'
         return
      end if
      if (kk(t) /= j) goto 101
      dj=cc(t)
      u(i)=dj-v(j)
      h=h+dj
   end do

1000 return

end subroutine

end module
