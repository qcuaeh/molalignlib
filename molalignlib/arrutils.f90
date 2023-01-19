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

module arrutils
use parameters
use settings

implicit none

contains

function idenperm(n)
! Purpose: Get an identity permutation
   integer, intent(in) :: n
   integer :: idenperm(n)
   integer :: i

   do i = 1, n
      idenperm(i) = i
   end do

end function

function inverperm(perm)
! Purpose: Get the inverse permutation
   integer, intent(in) :: perm(:)
   integer :: inverperm(size(perm))
   integer :: i

   do i = 1, size(perm)
      inverperm(perm(i)) = i
   end do

end function

logical function isperm(arr, n) result(flag)
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

end module
