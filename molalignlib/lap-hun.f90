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

module lap
use parameters

implicit none

contains

subroutine minperm(n, p, q, pq, perm, dist)

   integer, intent(in) :: n
   real(wp), intent(in) :: p(3, n), q(3, n)
   real(wp), intent(in) :: pq(n, n)
   integer, intent(out) :: perm(n)
   real(wp), intent(out) :: dist

   integer :: i, j
!  real(wp) c(n, n)
   real(wp), allocatable :: c(:, :)
   allocate(c(n, n))

   do i = 1, n
      do j = 1, n
         c(i, j) = sum((p(:, i) - q(:, j))**2) + pq(i, j)
      end do
   end do

   call assndx(1, c, n, n, perm, dist)

end subroutine


subroutine assndx(mode, a, n, m, k, sum)

!https://wp.csiro.au/alanmiller/assndx.f90
! Code converted using TO_F90 by Alan Miller
! Date: 2002-03-06  Time: 08:36:31

! Converted with permission, from the F77 code in the CERN MATHLIB library.
! $Id: assndx.F,v 1.1.1.1 1996/04/01 15:02:49 mclareni Exp $

! $Log: assndx.F,v $
! Revision 1.1.1.1  1996/04/01 15:02:49  mclareni
! Mathlib gen/H (H301)
! Author: F. Bourgeois, 15 February 1994

! N.B. Arguments IDA, IW & IDW have been removed.

! If MODE = 1, then it finds k(1), k(2), ..., k(n) to minimize
!        S = Sum(i=1, .., n) a(i, k(i))
! If MODE = 2,  then it finds k(1), k(2), ..., k(m) to minimize
!        S = Sum(j=1, .., m) a(k(j), j)
! given the array a(n,m).

! References:
! Munkres, J. (1957) `Algorithms for the assignment and transportation problems',
!                    J. SIAM, vol.5, 32-38.
! Silver, R. (1960) `An algorithm for the assignment problem', Comm. ACM, vol.3,
!                   605-606.   The algorithm (CACM 27) is in Algol.

integer, intent(in)   :: mode
real(wp), intent(in out)  :: a(:,:)
integer, intent(in)   :: n
integer, intent(in)   :: m
integer, intent(out)  :: k(:)
real(wp), intent(out)     :: sum

logical  :: lsw
integer  :: i, icbl, icl, icl0, iflag, imax, imin, ipp, irl, irs, &
                  j, j1, jsv, new
real(wp) :: rmin
integer, allocatable  :: iw(:,:)

if (n < 1 .or. m < 1) then
   write(*, '(a, 2i8)') ' ** error in call to assndx; m, n = ', m, n
   return
end if

imax = max(n,m)
imin = min(n,m)
allocate( iw(imax,6) )
sum = 0.0
if (n <= m) then
   do  i = 1, n
      rmin = a(i,1)
      do  j = 1, m
         rmin = min(rmin, a(i,j))
      end do
      sum = sum + rmin
      a(i,1:m) = a(i,1:m) - rmin
   end do
end if
if (n >= m) then
   do  j = 1, m
      rmin = a(1,j)
      do  i = 1, n
         rmin = min(rmin,a(i,j))
      end do
      sum = sum + rmin
      a(1:n,j) = a(1:n,j) - rmin
   end do
end if

do  i = 1, imax
   k(i) = 0
   iw(i,1) = 0
end do

loop90:  do  i = 1, n
   do  j = 1, m
      if (a(i,j)+iw(j,1) == 0) then
         k(i) = j
         iw(j,1) = i
         cycle loop90
      end if
   end do
end do loop90

100 iflag = n
irl = 0
icl = 0
irs = 1

do  i = 1, n
   iw(i,5) = 0
   if (k(i) == 0) then
      irl = irl + 1
      iw(irl,6) = i
      iw(i,5) = -1
      iflag = iflag - 1
   end if
end do
if (iflag == imin) then
   if (mode == 2) k(1:imax) = iw(1:imax,1)
   return
end if

iw(1:m,4) = 0

140 i = iw(irs,6)
irs = irs + 1
do  j = 1, m
   if (a(i,j)+iw(j,4) == 0) then
      iw(j,4) = i
      icl = icl + 1
      iw(icl,2) = j
      new = iw(j,1)
      if (new == 0) then
         j1 = j
         do
            iw(j1,1) = iw(j1,4)
            i = iw(j1,4)
            if (k(i) == 0) then
               k(i) = j1
               go to 100
            end if
            jsv = j1
            j1 = k(i)
            k(i) = jsv
         end do
      end if
      irl = irl + 1
      iw(irl,6) = new
      iw(new,5) = i
   end if
end do
if (irs <= irl) go to 140

lsw = .true.
icl0 = icl
icbl = 0
do  j = 1, m
   if (iw(j,4) == 0) then
      icbl = icbl + 1
      iw(icbl,3) = j
   end if
end do
rmin = a(iw(1,6),iw(1,3))
do  i = 1, irl
   do  j = 1, icbl
      rmin = min(rmin, a(iw(i,6), iw(j,3)))
   end do
end do
sum = sum + rmin * (irl+icbl-imax)

do  i = 1, n
   if (iw(i,5) == 0) then
      do  ipp = 1, icl0
         a(i,iw(ipp,2)) = a(i,iw(ipp,2)) + rmin
      end do
      cycle
   end if
   do  ipp = 1, icbl
      new = iw(ipp,3)
      a(i,new) = a(i,new) - rmin
      if (lsw.and.a(i,new)+iw(new,4) == 0) then
         iw(new,4) = i
         if (iw(new,1) == 0) then
            j1 = new
            lsw = .false.
         else
            icl = icl + 1
            iw(icl,2) = new
            irl = irl + 1
            iw(irl,6) = iw(new,1)
         end if
      end if
   end do
end do

if (lsw) then
   do  i = icl0 + 1, icl
      iw(iw(iw(i,2),1),5) = iw(i,2)
   end do
   go to 140
else
   do
      iw(j1,1) = iw(j1,4)
      i = iw(j1,4)
      if (k(i) == 0) then
         k(i) = j1
         go to 100
      end if
      jsv = j1
      j1 = k(i)
      k(i) = jsv
   end do
end if

return
end subroutine assndx

end module
