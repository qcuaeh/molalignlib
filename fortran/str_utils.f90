! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

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

module str_utils
use parameters
implicit none
private
public NUMCHAR
public ALPHACHAR
public LOWERCHAR
public UPPERCHAR
public str
public int
public lowercase
public uppercase

interface int
   module procedure int_str
end interface

interface str
   module procedure str_int
   module procedure str_real
end interface

character(*), parameter :: NUMCHAR = '1234567890'
character(*), parameter :: LOWERCHAR = 'abcdefghijklmnopqrstuvwxyz'
character(*), parameter :: UPPERCHAR = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
character(*), parameter :: ALPHACHAR = LOWERCHAR // UPPERCHAR

contains

elemental function lowercase(x) result(y)
   character(*), intent(in) :: x
   character(len(x)) :: y
   integer(ik) :: i, j
   do j = 1, len(x)
      i = index(UPPERCHAR, x(j:j))
      if (i > 0) then
         y(j:j) = LOWERCHAR(i:i)
      else
         y(j:j) = x(j:j)
      end if
   end do
end function

elemental function uppercase(x) result(y)
   character(*), intent(in) :: x
   character(len(x)) :: y
   integer(ik) :: i, j
   do j = 1, len(x)
      i = index(LOWERCHAR, x(j:j))
      if (i > 0) then
         y(j:j) = UPPERCHAR(i:i)
      else
         y(j:j) = x(j:j)
      end if
   end do
end function

function str_int(x) result(str)
   integer(ik), intent(in) :: x
   character(:), allocatable :: str
   character(10) :: temp
   write(temp, '(I0)') x
   str = trim(temp)
end function

function str_real(x) result(str)
   real(rk), intent(in) :: x
   character(:), allocatable :: str
   integer(ik), parameter :: WIDTH=DECIMAL_PLACES+10
   character(WIDTH) :: temp
   write(temp, '(F'//str_int(WIDTH)//'.'//str_int(DECIMAL_PLACES)//')') x
   str = trim(adjustl(temp))
end function

function int_str(x) result(intval)
   character(*), intent(in) :: x
   integer(ik) :: intval, stat
   read (x,*,iostat=stat) intval
   if (stat /= 0) then
      print *, x
      write (stderr, '(A,A)') 'String is not an integer'
      stop
   end if
end function

end module
