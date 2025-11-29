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
use parameters

implicit none

private
public quicksort
public sort_pairs

interface quicksort
   module procedure int_quicksort
   module procedure real_quicksort
   module procedure int_quicksort_all
   module procedure real_quicksort_all
end interface

contains

subroutine int_quicksort_all(x)
   integer, intent(inout) :: x(:)
   call int_quicksort(x, 1, size(x))
end subroutine

subroutine real_quicksort_all(x)
   real(rk), intent(inout) :: x(:)
   call real_quicksort(x, 1, size(x))
end subroutine

recursive subroutine int_quicksort(x, m, n)
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

   if (m < i-1) call int_quicksort(x, m, i-1)
   if (j+1 < n) call int_quicksort(x, j+1, n)
end subroutine

recursive subroutine real_quicksort(x, m, n)
   integer, intent(in) :: m, n
   real(rk), intent(inout) :: x(:)

   integer :: i, j
   real(rk) :: xi, xp

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

   if (m < i-1) call real_quicksort(x, m, i-1)
   if (j+1 < n) call real_quicksort(x, j+1, n)
end subroutine

subroutine sort_pairs(pairs)
! Sort pairs lexicographically: first by atom1, then by atom2
   integer, dimension(:,:), intent(inout) :: pairs
   integer :: n, i, j, temp1, temp2
   logical :: swapped

   n = size(pairs, 2)

   ! Simple bubble sort (can be replaced with quicksort if needed)
   do i = n, 2, -1
      swapped = .FALSE.
      do j = 1, i - 1
         if (pairs(1, j) > pairs(1, j + 1) .or. &
             (pairs(1, j) == pairs(1, j + 1) .and. pairs(2, j) > pairs(2, j + 1))) then
            temp1 = pairs(1, j)
            temp2 = pairs(2, j)
            pairs(1, j) = pairs(1, j + 1)
            pairs(2, j) = pairs(2, j + 1)
            pairs(1, j + 1) = temp1
            pairs(2, j + 1) = temp2
            swapped = .TRUE.
         end if
      end do
      if (.not. swapped) exit
   end do
end subroutine

end module
