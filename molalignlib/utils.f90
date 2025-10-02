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

module utils
use parameters
implicit none
private
public str
public uniform_weights
public lowercase
public uppercase

interface str
   module procedure str_int
   module procedure str_real
end interface

character(26), parameter :: LOWERCHARSET = 'abcdefghijklmnopqrstuvwxyz'
character(26), parameter :: UPPERCHARSET = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

contains

function str_int(x) result(str)
   integer, intent(in) :: x
   character(:), allocatable :: str
   character(10) :: temp
   write(temp, '(I0)') x
   str = trim(temp)
end function

function str_real(x) result(str)
   real(rk), intent(in) :: x
   character(:), allocatable :: str
   integer, parameter :: WIDTH=DECIMAL_PLACES+10
   character(WIDTH) :: temp
   write(temp, '(F'//str_int(WIDTH)//'.'//str_int(DECIMAL_PLACES)//')') x
   str = trim(adjustl(temp))
end function

function uniform_weights(x, n) result(xn)
   real(rk), intent(in) :: x
   integer, intent(in) :: n
   real(rk), dimension(:), allocatable :: xn
   allocate(xn(n), source=x)
end function

function lowercase(x) result(l)
   character(*), intent(in) :: x
   character(len(x)) :: l
   integer :: i, j
   l = x
   do j = 1, len(x)
      i = index(UPPERCHARSET, x(j:j))
      if (i > 0) l(j:j) = LOWERCHARSET(i:i)
   end do
end function

function uppercase(x) result(u)
   character(*), intent(in) :: x
   character(len(x)) :: u
   integer :: i, j
   u = x
   do j = 1, len(x)
      i = index(LOWERCHARSET, x(j:j))
      if (i > 0) u(j:j) = UPPERCHARSET(i:i)
   end do
end function

end module
