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

module sorting
use kinds

implicit none

private
public sort
public sorted
public sorted_order

interface sort
   module procedure intquicksort
   module procedure realquicksort
end interface

interface sorted
   module procedure int_sorted
   module procedure int_sorted_partial
   module procedure real_sorted_partial
end interface

interface sorted_order
   module procedure int_sorted_order
   module procedure int_sorted_order_partial
   module procedure real_sorted_order_partial
   module procedure char_sorted_order_partial
end interface

contains

function int_sorted(x) result(y)
   integer, intent(in) :: x(:)
   integer, allocatable :: y(:)
   y = x(:)
   call intquicksort(y, 1, size(y))
end function

function int_sorted_partial(x, n) result(y)
   integer, intent(in) :: n
   integer, intent(in) :: x(:)
   integer :: y(n)
   y(1:n) = x(1:n)
   call intquicksort(y, 1, n)
end function

function real_sorted_partial(x, n) result(y)
   integer, intent(in) :: n
   real(wp), intent(in) :: x(:)
   real(wp), allocatable :: y(:)
   y = x(1:n)
   call realquicksort(y, 1, n)
end function

function int_sorted_order(x) result(o)
   integer, intent(in) :: x(:)
   integer :: i, o(size(x)), t((size(x)+1)/2)
   o = [(i, i=1, size(x))]
   call intmergesort(x, o, size(x), t)
end function

function int_sorted_order_partial(x, n) result(o)
   integer, intent(in) :: n
   integer, intent(in) :: x(:)
   integer :: i, o(n), t((n+1)/2)
   o = [(i, i=1, n)]
   call intmergesort(x, o, n, t)
end function

function real_sorted_order_partial(x, n) result(o)
   integer, intent(in) :: n
   real(wp), intent(in) :: x(:)
   integer :: i, o(n), t((n+1)/2)
   o = [(i, i=1, n)]
   call realmergesort(x, o, n, t)
end function

function char_sorted_order_partial(x, n) result(o)
   integer, intent(in) :: n
   character(*), intent(in) :: x(:)
   integer :: i, o(n), t((n+1)/2)
   o = [(i, i=1, n)]
   call charmergesort(x, o, n, t)
end function

recursive subroutine intquicksort(x, m, n)
   integer, intent(in) :: m, n
   integer, intent(inout) :: x(:)

   integer :: i, j
   integer :: xi, xp

   xp = x((m+n)/2)
   i = m
   j = n

   do
      do while (x(i) < xp)
         i = i + 1
      end do
      do while (xp < x(j))
         j = j - 1
      end do
      if (i >= j) exit
      xi = x(i); x(i) = x(j); x(j) = xi
      i = i + 1
      j = j - 1
   end do

   if (m < i-1) call intquicksort(x, m, i-1)
   if (j+1 < n) call intquicksort(x, j+1, n)
end subroutine

recursive subroutine realquicksort(x, m, n)
   integer, intent(in) :: m, n
   real(wp), intent(inout) :: x(:)

   integer :: i, j
   real(wp) :: xi, xp

   xp = x((m+n)/2)
   i = m
   j = n

   do
      do while (x(i) < xp)
         i = i + 1
      end do
      do while (xp < x(j))
         j = j - 1
      end do
      if (i >= j) exit
      xi = x(i); x(i) = x(j); x(j) = xi
      i = i + 1
      j = j - 1
   end do

   if (m < i-1) call realquicksort(x, m, i-1)
   if (j+1 < n) call realquicksort(x, j+1, n)
end subroutine

subroutine intmerge(x, a, na, b, nb, c, nc)
   integer, intent(in) :: na, nb, nc
   integer, intent(in) :: x(:)
   integer, intent(in) :: b(nb)
   integer, intent(inout) :: a(na), c(nc)
    
   integer :: i, j, k
    
   i = 1; j = 1; k = 1;
   do while(i <= na .and. j <= nb)
      if (x(a(i)) <= x(b(j))) then
         c(k) = a(i)
         i = i+1
      else
         c(k) = b(j)
         j = j+1
      end if
      k = k + 1
   end do
   do while (i <= na)
      c(k) = a(i)
      i = i + 1
      k = k + 1
   end do
end subroutine
 
subroutine realmerge(x, a, na, b, nb, c, nc)
   integer, intent(in) :: na, nb, nc
   real(wp), intent(in) :: x(:)
   integer, intent(in) :: b(nb)
   integer, intent(inout) :: a(na), c(nc)
    
   integer :: i, j, k
    
   i = 1; j = 1; k = 1;
   do while(i <= na .and. j <= nb)
      if (x(a(i)) <= x(b(j))) then
         c(k) = a(i)
         i = i+1
      else
         c(k) = b(j)
         j = j+1
      end if
      k = k + 1
   end do
   do while (i <= na)
      c(k) = a(i)
      i = i + 1
      k = k + 1
   end do
end subroutine
 
recursive subroutine realmergesort(x, o, n, t)
   integer, intent(in) :: n
   real(wp), intent(in) :: x(:)
   integer, intent(inout) :: o(n)
   integer, intent(out) :: t((n+1)/2)
    
   integer :: i, o1
    
   if (n < 2) return
   if (n == 2) then
      if (x(o(1)) > x(o(2))) then
         o1 = o(1); o(1) = o(2); o(2) = o1
      end if
      return
   end if      
   i = (n+1)/2
    
   call realmergesort(x, o, i, t)
   call realmergesort(x, o(i+1), n-i, t)
    
   if (x(o(i)) > x(o(i+1))) then
      t(1:i) = o(1:i)
      call realmerge(x, t, i, o(i+1), n-i, o, n)
   end if
end subroutine

recursive subroutine intmergesort(x, o, n, t)
   integer, intent(in) :: n
   integer, intent(in) :: x(:)
   integer, intent(inout) :: o(n)
   integer, intent(out) :: t((n+1)/2)
    
   integer :: i, o1
    
   if (n < 2) return
   if (n == 2) then
      if (x(o(1)) > x(o(2))) then
         o1 = o(1); o(1) = o(2); o(2) = o1
      end if
      return
   end if      
   i = (n+1)/2
    
   call intmergesort(x, o, i, t)
   call intmergesort(x, o(i+1), n-i, t)
    
   if (x(o(i)) > x(o(i+1))) then
      t(1:i) = o(1:i)
      call intmerge(x, t, i, o(i+1), n-i, o, n)
   end if
end subroutine

subroutine charmerge(x, a, na, b, nb, c, nc)
   integer, intent(in) :: na, nb, nc
   character(*), intent(in) :: x(:)
   integer, intent(in) :: b(nb)
   integer, intent(inout) :: a(na), c(nc)
    
   integer :: i, j, k
    
   i = 1; j = 1; k = 1;
   do while(i <= na .and. j <= nb)
      if (x(a(i)) <= x(b(j))) then
         c(k) = a(i)
         i = i+1
      else
         c(k) = b(j)
         j = j+1
      end if
      k = k + 1
   end do
   do while (i <= na)
      c(k) = a(i)
      i = i + 1
      k = k + 1
   end do
end subroutine
 
recursive subroutine charmergesort(x, o, n, t)
   integer, intent(in) :: n
   character(*), intent(in) :: x(:)
   integer, intent(inout) :: o(n)
   integer, intent(out) :: t((n+1)/2)
    
   integer :: i, o1
    
   if (n < 2) return
   if (n == 2) then
      if (x(o(1)) > x(o(2))) then
         o1 = o(1); o(1) = o(2); o(2) = o1
      end if
      return
   end if      
   i = (n+1)/2
    
   call charmergesort(x, o, i, t)
   call charmergesort(x, o(i+1), n-i, t)
    
   if (x(o(i)) > x(o(i+1))) then
      t(1:i) = o(1:i)
      call charmerge(x, t, i, o(i+1), n-i, o, n)
   end if
end subroutine

end module
