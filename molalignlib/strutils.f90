! MolAlign
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

module strutils
use parameters
use settings

implicit none

private
public str
public lower
public upper
public basename
public getext

interface str
   module procedure int2str
   module procedure real2str
end interface

contains

function lower(str)
   character(*), intent(in) :: str
   character(26), parameter :: lowercase = 'abcdefghijklmnopqrstuvwxyz'
   character(26), parameter :: uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character(len(str)) :: lower
   integer :: i, j
   lower = str
   do j = 1, len(str)
      i = index(uppercase, str(j:j))
      if (i > 0) lower(j:j) = lowercase(i:i)
   end do
end function

function upper(str)
   character(*), intent(in) :: str
   character(26), parameter :: lowercase = 'abcdefghijklmnopqrstuvwxyz'
   character(26), parameter :: uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character(len(str)) :: upper
   integer :: i, j
   upper = str
   do j = 1, len(str)
      i = index(lowercase, str(j:j))
      if (i > 0) upper(j:j) = uppercase(i:i)
   end do
end function

function int2str(x) result(strx)
   integer, intent(in) :: x
   character(floor(log10(real(x, wp))) + 1) strx
   write (strx, '(i0)') x
end function

function real2str(x) result(strx)
   real(wp), intent(in) :: x
   character(floor(log10(x)) + 6) strx
   write (strx, '(f0.4)') x
end function

function basename(filepath)
   character(*), intent(in) :: filepath
   character(len(filepath)) basename
   integer pos
   pos = index(trim(filepath), '/', back=.true.)
   if (pos /= 0) then
      basename = filepath(pos+1:len_trim(filepath))
   else
      basename = trim(filepath)
   end if
end function

function getext(filename)
   character(*), intent(in) :: filename
   character(len(filename)) getext
   integer pos
   pos = index(trim(filename), '.', back=.true.)
   if (pos > 1 .and. pos < len_trim(filename)) then
      getext = filename(pos+1:len_trim(filename))
   else
      getext = ''
   end if
end function

end module
