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

module fileio
use parameters
use molecule
use adjacency
use readmol
use writemol

implicit none

contains

subroutine open2read(filepath, unit, extension)
   character(*), intent(in) :: filepath
   integer, intent(out) :: unit
   character(:), allocatable, intent(out) :: extension
   integer :: stat

   extension = get_extension(filepath)

   open(newunit=unit, file=filepath, action='read', status='old', iostat=stat)
   if (stat /= 0) then
      write (stderr, '(a,1x,a,1x,a)') 'Error opening', filepath, 'for reading'
      stop
   end if

end subroutine

subroutine open2write(filepath, unit, extension)
   character(*), intent(in) :: filepath
   integer, intent(out) :: unit
   character(:), allocatable, intent(out) :: extension
   integer :: stat

   extension = get_extension(filepath)

   open(newunit=unit, file=filepath, action='write', status='replace', iostat=stat)
   if (stat /= 0) then
      write (stderr, '(a,1x,a,1x,a)') 'Error opening', filepath, 'for writing'
      stop
   end if

end subroutine

subroutine read_file(unit, fmtin, mol)
   integer, intent(in) :: unit
   character(*), intent(in) :: fmtin
   type(mol_type), intent(out) :: mol

   select case (fmtin)
   case ('xyz')
      call readxyz(unit, mol)
   case ('mol2')
      call readmol2(unit, mol)
   case default
      write (stderr, '(a,1x,a)') 'Invalid format:', fmtin
      stop
   end select

end subroutine

subroutine write_file(unit, fmtout, mol)
   integer, intent(in) :: unit
   character(*), intent(in) :: fmtout
   type(mol_type), intent(in) :: mol

   select case (fmtout)
   case ('xyz')
      call writexyz(unit, mol)
   case ('mol2')
      call writemol2(unit, mol)
   case default
      write (stderr, '(a,1x,a)') 'Invalid format:', fmtout
      stop
   end select

   flush(stderr)

end subroutine

function get_extension(filepath) result(extension)
   character(*), intent(in) :: filepath
   character(:), allocatable :: filename, basename, extension
   integer :: pos
   pos = index(filepath, '/', back=.true.)
   filename = filepath(pos+1:)
   if (len(filename) == 0) then
      write (stderr, '(a,1x,a)') 'Missing file name'
      stop
   end if
   pos = index(filename, '.', back=.true.)
   if (pos /= 0) then
      basename = filename(:pos-1)
      extension = filename(pos+1:)
   else
      basename = filename
      extension = ''
   end if
   if (len(basename) == 0) then
      write (stderr, '(a,1x,a)') 'Missing base name of', filename
      stop
   end if
   if (len(extension) == 0) then
      write (stderr, '(a,1x,a)') 'Missing extension of', filename
      stop
   end if
end function

end module
