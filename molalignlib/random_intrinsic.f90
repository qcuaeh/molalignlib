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

module random
use kinds
use flags

implicit none

contains

subroutine random_initialize()

   if (test_flag) then
      call random_init(.true., .true.)
   else
      call random_init(.false., .true.)
   end if

end subroutine

function randvec() result(r)
   real(rk) :: r(3)
   call random_number(r)
end function

subroutine shuffle(a)
! Fisher-Yates shuffle
   integer, intent(inout) :: a(:)
   integer :: i, j, temp
   real(rk) :: r
   do i = size(a), 2, -1
      call random_number(r)
      j = int(r * i) + 1
      temp = a(j)
      a(j) = a(i)
      a(i) = temp
   end do
end subroutine

end module
