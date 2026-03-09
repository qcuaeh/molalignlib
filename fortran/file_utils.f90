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

module file_utils
use parameters
use str_utils
use chemdata
implicit none
private
public open2read
public open2write

contains

subroutine split_path(path, dirname, filename)
   character(*), intent(in) :: path
   character(:), allocatable, intent(out) :: dirname, filename
   ! Local variables
   integer(ik) :: pos
   pos = index(path, '/', back=.TRUE.)
   if (pos /= 0) then
      dirname = path(:pos-1)
      filename = path(pos+1:)
      if (len(filename) == 0) then
         write (stderr, '(A,1X,A)') 'Invalid path'
         stop
      end if
   else
      dirname = ''
      filename = path
   end if
end subroutine

subroutine split_filename(filename, basename, extension)
   character(*), intent(in) :: filename
   character(:), allocatable, intent(out) :: basename, extension
   ! Local variables
   integer(ik) :: pos
   pos = index(filename, '.', back=.TRUE.)
   if (pos /= 0) then
      basename = filename(:pos-1)
      extension = filename(pos+1:)
      if (len(basename) == 0 .or. len(extension) == 0) then
         write (stderr, '(A,1X,A)') 'Invalid file name', filename
         stop
      end if
   else
      basename = ''
      extension = filename
   end if
end subroutine

subroutine open2read(path, extension, unit)
   character(*), intent(in) :: path
   character(:), allocatable, intent(out) :: extension
   integer(ik), intent(out) :: unit
   ! Local variables
   character(:), allocatable :: dirname, filename, basename
   integer(ik) :: stat

   call split_path(path, dirname, filename)
   call split_filename(filename, basename, extension)

   if (basename == '') then
      unit = stdout
      return
   end if

   open(newunit=unit, file=path, action='read', status='old', iostat=stat)
   if (stat /= 0) then
      write (stderr, '(A,1X,A,1X,A)') 'opening', path, 'for reading'
      stop
   end if
end subroutine

subroutine open2write(path, extension, unit)
   character(*), intent(in) :: path
   character(:), allocatable, intent(out) :: extension
   integer(ik), intent(out) :: unit
   ! Local variables
   character(:), allocatable :: dirname, filename, basename
   integer(ik) :: stat

   call split_path(path, dirname, filename)
   call split_filename(filename, basename, extension)

   if (basename == '') then
      unit = stdout
      return
   end if

   open(newunit=unit, file=path, action='write', status='replace', iostat=stat)
   if (stat /= 0) then
      write (stderr, '(A,1X,A,1X,A)') 'Can''t open', path, 'for writing'
      stop
   end if
end subroutine

end module
