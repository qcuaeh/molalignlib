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

module writemol
use parameters
use molecule
use strutils
use chemdata

implicit none

contains

subroutine writexyz(unit, mol)
   integer, intent(in) :: unit
   type(mol_type) :: mol
   ! Local varibles
   type(atom_type), allocatable :: atoms(:)
   integer :: i

   atoms = mol%atoms

   write (unit, '(i0)') size(atoms)
   write (unit, '(a)') mol%title

   do i = 1, size(atoms)
      write (unit, '(a,3(2x,f12.6))') element_symbols(atoms(i)%elnum), atoms(i)%coords
   end do

end subroutine

subroutine writemol2(unit, mol)
   integer, intent(in) :: unit
   type(mol_type) :: mol
   ! Local variables
   real(rk), allocatable :: coords(:, :)
   type(bond_type), allocatable :: bonds(:)
   type(atom_type), allocatable :: atoms(:)
   integer :: i

   atoms = mol%atoms
   bonds = mol%get_bonds()
   coords = mol%get_coords()

   write (unit, '(a)') '@<TRIPOS>MOLECULE'
   if (mol%title /= '') then
      write (unit, '(a)') mol%title
   else
      write (unit, '(a)') 'Untitled'
   end if
   write (unit, '(5(i4,1x))') size(atoms), size(bonds), 0, 0, 0
   write (unit, '(a)') 'SMALL'
   write (unit, '(a)') 'NO_CHARGES'
   write (unit, '(a)') '@<TRIPOS>ATOM'

   do i = 1, size(atoms)
      write (unit, '(i4,2x,a2,3(1x,f12.6),2x,a4,1x,i2,1x,a4,1x,f7.3)') &
         i, element_symbols(atoms(i)%elnum), coords(:, i), atomtype(atoms(i)), 1, 'MOL1', 0.
   end do

   write (unit, '(a)') '@<TRIPOS>BOND'

   do i = 1, size(bonds)
      write (unit, '(i4,1x,2(1x,i4),1x,a2)') i, bonds(i)%atomidx1, bonds(i)%atomidx2, '1'
   end do

end subroutine

function atomtype(atom)
   type(atom_type), intent(in) :: atom
   ! Local variables
   integer :: nadj, hybnum
   character(4) :: elsym, atomtype

   nadj = size(atom%adjlist)
   elsym = element_symbols(atom%elnum)

   select case (atom%elnum)
   case (6)
      hybnum = nadj - 1
      atomtype = trim(elsym)//'.'//str(hybnum)
   case (7, 15)
      hybnum = min(nadj, 3)
      atomtype = trim(elsym)//'.'//str(hybnum)
   case (8, 16)
      hybnum = min(nadj + 1, 3)
      atomtype = trim(elsym)//'.'//str(hybnum)
   case default
      atomtype = elsym
   end select

end function

end module
