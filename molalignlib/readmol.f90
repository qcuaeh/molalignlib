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

module readmol
use stdio
use bounds
use strutils
use chemdata
use moltypes

implicit none

contains

subroutine readxyz(unit, mol)
   integer, intent(in) :: unit
   type(Molecule), intent(out) :: mol
   character(ll) :: buffer

   integer :: i

   read (unit, *, end=99) mol%natom
   mol%nbond = 0

   allocate(mol%atoms(mol%natom))

   read (unit, '(a)', end=99) buffer
   mol%title = trim(buffer)

   do i = 1, mol%natom
      read (unit, *, end=99) buffer, mol%atoms(i)%coords(:)
      mol%atoms(i)%label = buffer
   end do

   allocate(mol%adjmat(mol%natom, mol%natom))

   mol%adjmat(:, :) = .false.

   return

   99 continue
   write (error_unit, '(a)') 'Unexpected end of file!'
   stop

end subroutine

subroutine readmol2(unit, mol)
   integer, intent(in) :: unit
   type(Molecule), intent(out) :: mol
   character(ll) :: buffer
   integer :: i, id

   do
      read (unit, '(a)', end=99) buffer
      if (buffer == '@<TRIPOS>MOLECULE') exit
   end do

   read (unit, '(a)', end=99) buffer
   mol%title = trim(buffer)
   read (unit, *, end=99) mol%natom, mol%nbond

   allocate(mol%atoms(mol%natom))
   allocate(mol%bonds(mol%nbond))

   do
      read (unit, '(a)', end=99) buffer
      if (buffer == '@<TRIPOS>ATOM') exit
   end do

   do i = 1, mol%natom
      read (unit, *, end=99) id, buffer, mol%atoms(i)%coords(:)
      mol%atoms(i)%label = buffer
   end do

   do
      read (unit, '(a)', end=99) buffer
      if (buffer == '@<TRIPOS>BOND') exit
   end do

   do i = 1, mol%nbond
      read (unit, *, end=99) id, mol%bonds(i)%atom1, mol%bonds(i)%atom2, mol%bonds(i)%order
   end do

   allocate(mol%adjmat(mol%natom, mol%natom))

   mol%adjmat(:, :) = .false.

   do i = 1, mol%nbond
      mol%adjmat(mol%bonds(i)%atom1, mol%bonds(i)%atom2) = .true.
      mol%adjmat(mol%bonds(i)%atom2, mol%bonds(i)%atom1) = .true.
   end do

   return

   99 continue
   write (error_unit, '(a)') 'Unexpected end of file!'
   stop

end subroutine

end module
