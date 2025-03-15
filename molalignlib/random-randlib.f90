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
use parameters
use globals
use randlib

implicit none

contains

subroutine random_initialize()

   if (test_flag) then
      call default_set_seeds()
   else
      call time_set_seeds()
   end if

end subroutine

function randvec() result(x)
   real(rk) :: x(3)
   x(1) = random_standard_uniform()
   x(2) = random_standard_uniform()
   x(3) = random_standard_uniform()
end function

subroutine shuffle(a)
! Fisher-Yates shuffle
   integer, dimension(:), intent(inout) :: a
   integer :: i, j, temp
   do i = size(a), 2, -1
!      j = int(random_standard_uniform() * i) + 1
      j = random_uniform_integer(1, i)
      temp = a(j)
      a(j) = a(i)
      a(i) = temp
   end do
end subroutine

function random_element(a) result(e)
   integer, dimension(:), intent(in) :: a
   integer :: i, e
   i = random_uniform_integer(1, size(a))
   e = a(i)
end function

end module
