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

module setutils

implicit none

contains

function intersection(list1, list2, hash_size)
! Find the intersection of two lists
   integer, intent(in) :: hash_size
   integer, intent(in) :: list1(:), list2(:)
   ! Local variables
   integer :: i, n
   logical, allocatable :: hash_table(:)
   integer, allocatable :: intersection(:)

   allocate (hash_table(hash_size))

   hash_table = .false.

   ! Add elements in list1 to hash table
   do i = 1, size(list1)
      hash_table(list1(i)) = .true.
   end do

   ! Count intersection elements
   n = 0
   do i = 1, size(list2)
      if (hash_table(list2(i))) then
         n = n + 1
      end if
   end do

   allocate (intersection(n))

   ! Store intersection elements
   n = 0
   do i = 1, size(list2)
      if (hash_table(list2(i))) then
         n = n + 1
         intersection(n) = list2(i)
      end if
   end do

end function

subroutine union_sorted(n1, list1, n2, list2, n3, list3)
! Find the union of two sorted lists without repeated elements
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

subroutine intersection_sorted(n1, list1, n2, list2, n3, list3)
! Find the intersection of two sorted lists without repeated elements
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
