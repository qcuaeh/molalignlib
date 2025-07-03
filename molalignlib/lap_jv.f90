module lap_jv
use parameters
implicit none
contains

subroutine jovosap(n,sz,cc,kk,first,x,y,u,v,h)
! This subroutine performs weighted bipartite matching for
! for a sparse non-negative integer weight matrix.

! Adapted from GMIN: A program for finding global minima
! Copyright (C) 1999-2006 David J. Wales
! The original source is
!     http://www.magiclogic.com/assignment.html
! A publication reference can be found on the above homepage and
! in a comment below

   integer :: n,sz
   integer :: kk(sz),first(n+1),x(n),y(n)
   integer :: i,i0,j,j0,j1,k,l,l0,t,t0,td,cnt
   integer :: lab(n),free(n),todo(n)
   integer(int64) :: cc(sz),u(n),v(n),d(n)
   integer(int64) :: h,v0,vj,dj,min
   integer(int64), parameter :: bigint = 1000000000000_int64
   logical :: ok(n)

! THIS SUBROUTINE SOLVES THE SPARSE LINEAR ASSIGNMENT PROBLEM
! ACCORDING 
!
! "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear   
!  Assignment Problems," Computing 38, 325-340, 1987
! 
! by
! 
! R. Jonker and A. Volgenant, University of Amsterdam.
!
! INPUT PARAMETERS :
! N = NUMBER OF ROWS AND COLUMNS
! C = WEIGHT MATRIX
!
! OUTPUT PARAMETERS
! X = COL ASSIGNED TO ROW
! Y = ROW ASSIGNED TO COL
! U = DUAL ROW VARIABLE
! V = DUAL COLUMN VARIABLE
! H = VALUE OF OPTIMAL SOLUTION
!
! INITIALIZATION

!   Next line added by tomaso@nada.kth.se, to enable detection
!   of solutions being equivalent to the initial guess
! If Y(:) is initialised to zero then we see segmentation faults if 
! a Y element is unset, etc.

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
!         print '(a,i6,a)','jovosap> warning b - matching failed'
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
!            print '(a,i6,a)','jovosap> warning c - matching failed'
            return
         end if
         if (y(j) == 0) goto 80
         ok(j)=.true.
      end do
! repeat until a free row has been found
60    if (k == 0) then
!         print '(a,i6,a)','jovosap> warning d - matching failed'
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
         print '(a,i6,a)','jovosap> warning d - atom ',i,' not matched - maximum number of neighbours too small?'
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
