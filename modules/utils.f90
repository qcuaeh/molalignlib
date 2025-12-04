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

module utils
use parameters
implicit none
private
public str
public lowercase
public uppercase
public parse_path
public uniform_weights

interface str
   module procedure str_int
   module procedure str_real
end interface

character(26), parameter :: UPCHARSET = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
character(26), parameter :: LOWCHARSET = 'abcdefghijklmnopqrstuvwxyz'

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
      i = index(UPCHARSET, x(j:j))
      if (i > 0) l(j:j) = LOWCHARSET(i:i)
   end do
end function

function uppercase(x) result(u)
   character(*), intent(in) :: x
   character(len(x)) :: u
   integer :: i, j
   u = x
   do j = 1, len(x)
      i = index(LOWCHARSET, x(j:j))
      if (i > 0) u(j:j) = UPCHARSET(i:i)
   end do
end function

subroutine parse_path(filepath, filetype)
   character(*), intent(in) :: filepath
   character(:), allocatable, intent(out) :: filetype
   ! Local variables
   character(:), allocatable :: basename
   character(:), allocatable :: dirname, filename
   integer :: pos
   pos = index(filepath, '/', back=.TRUE.)
   if (pos /= 0) then
      dirname = filepath(:pos-1)
      basename = filepath(pos+1:)
      if (len(basename) == 0) then
         write (stderr, '(A,1X,A)') 'STOP File name is missing'
         stop
      end if
   else
      dirname = '.'
      basename = filepath
   end if
   pos = index(basename, '.', back=.TRUE.)
   if (pos /= 0) then
      filename = basename(:pos-1)
      filetype = basename(pos+1:)
      if (len(filename) == 0 .or. len(filetype) == 0) then
         write (stderr, '(A,1X,A)') 'STOP Invalid file name', basename
         stop
      end if
   else
      write (stderr, '(A,1X,A)') 'STOP File extension is missing', basename
      stop
   end if
end subroutine

end module
