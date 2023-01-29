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

module discrete
use parameters
use settings

implicit none

contains

function identityperm(n)
! Purpose: Get an identity permutation
   integer, intent(in) :: n
   integer :: identityperm(n)
   integer :: i

   do i = 1, n
      identityperm(i) = i
   end do

end function

function inverseperm(perm)
! Purpose: Get the inverse permutation
   integer, intent(in) :: perm(:)
   integer :: inverseperm(size(perm))
   integer :: i

   do i = 1, size(perm)
      inverseperm(perm(i)) = i
   end do

end function

logical function is_perm(arr, n) result(flag)
! Purpose: Check if array is a valid permutation
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: arr
   integer :: i, j

   flag = .true.
   do i = 1, n
      do j = 1, n
         if (arr(j) == i) exit
      end do
      if (j == n + 1) then
         flag = .false.
         exit
      end if
   end do

end function

subroutine union(n1, list1, n2, list2, n3, list3)
! Purpose: Find the union of two sorted lists without repeated elements
   integer, intent(in) :: n1, n2, list1(:), list2(:)
   integer, intent(out) :: n3, list3(:)
   integer :: i1, i2

   i1 = 1
   i2 = 1
   n3 = 0

   do while (i1 <= n1 .and. i2 <= n2)
      n3 = n3 + 1
      if (list1(i1) < list2(i2)) then
         list3(n3) = list1(i1)
         i1 = i1 + 1
      else if (list1(i1) > list2(i2)) then
         list3(n3) = list2(i2)
         i2 = i2 + 1
      else
         list3(n3) = list2(i2)
         i1 = i1 + 1
         i2 = i2 + 1
      end if
   end do

   do while (i1 <= n1)
      n3 = n3 + 1
      list3(n3) = list1(i1)
      i1 = i1 + 1
   end do

   do while (i2 <= n2)
      n3 = n3 + 1
      list3(n3) = list2(i2)
      i2 = i2 + 1
   end do

end subroutine

subroutine intersection(n1, list1, n2, list2, n3, list3)
! Purpose: Find the intersection of two sorted lists without repeated elements
   integer, intent(in) :: n1, n2, list1(:), list2(:)
   integer, intent(out) :: n3, list3(:)
   integer :: i1, i2

   i1 = 1
   i2 = 1
   n3 = 0

   do while (i1 <= n1 .and. i2 <= n2)
      if (list1(i1) < list2(i2)) then
         i1 = i1 + 1
      else if (list1(i1) > list2(i2)) then
         i2 = i2 + 1
      else
         n3 = n3 + 1
         list3(n3) = list2(i2)
         i1 = i1 + 1
         i2 = i2 + 1
      end if
   end do

end subroutine

end module
