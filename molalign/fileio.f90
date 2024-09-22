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
use stdio
use types
use readmol
use writemol
use adjacency
use assorting
use tracking

implicit none

contains

subroutine open2read(filepath, unit, fileext)
   character(*), intent(in) :: filepath
   integer, intent(out) :: unit
   character(:), allocatable, intent(out) :: fileext
   integer :: stat

   if (len(filepath) == 0) then
      write (stderr, '(a)') 'Error: File path is empty'
      stop
   end if

   fileext = baseext(filepath)

   if (len(fileext) == 0) then
      write (stderr, '(a,1x,a)') 'Missing file extension:', filepath
      stop
   end if

   open(newunit=unit, file=filepath, action='read', status='old', iostat=stat)
   if (stat /= 0) then
      write (stderr, '(a,1x,a,1x,a)') 'Error opening', filepath, 'for reading'
      stop
   end if

end subroutine

subroutine open2write(filepath, unit, fileext)
   character(*), intent(in) :: filepath
   integer, intent(out) :: unit
   character(:), allocatable, intent(out) :: fileext
   integer :: stat

   fileext = baseext(filepath)

   if (len(fileext) == 0) then
      write (stderr, '(a,1x,a)') 'Missing file extension:', filepath
      stop
   end if

   open(newunit=unit, file=filepath, action='write', status='replace', iostat=stat)
   if (stat /= 0) then
      write (stderr, '(a,1x,a,1x,a)') 'Error opening', filepath, 'for writing'
      stop
   end if

end subroutine

subroutine readfile(unit, fmtin, mol)
   integer, intent(in) :: unit
   character(*), intent(in) :: fmtin
   type(t_mol), intent(out) :: mol

   select case (fmtin)
   case ('xyz')
      call readxyz(unit, mol)
   case ('mol2')
      call readmol2(unit, mol)
   case default
      write (stderr, '(a,1x,a)') 'Invalid format:', fmtin
      stop
   end select

   call set_bonds(mol)
   call find_molfrags(mol)
   call assort_atoms(mol)
   call set_equiv_atoms(mol)
   call assort_neighbors(mol)

end subroutine

subroutine writefile(unit, fmtout, mol)
!   integer, intent(in) :: unit, natom
!   integer, dimension(:), intent(in) :: elnums
!   real(rk), dimension(:, :), intent(in) :: coords
!   logical, dimension(:, :), intent(in) :: adjmat
!   character(*), intent(in) :: title, fmtout
   integer, intent(in) :: unit
   character(*), intent(in) :: fmtout
   type(t_mol), intent(in) :: mol

   integer :: i, j
   integer :: nbond
   integer :: bonds(2, mol%natom*maxcoord)

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

end module
