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
use parameters
use globals
use strutils
use chemdata
use chemutils
use molecule

implicit none

contains

subroutine readxyz(unit, mol)
   integer, intent(in) :: unit
   type(mol_type), intent(out) :: mol
   ! Local variables
   character(ll) :: buffer
   character(wl) :: tag
   integer :: i, num_atoms, elnum, label
   real(rk) :: coords(3)

   read (unit, *, end=99) num_atoms
   allocate (mol%atoms(num_atoms))
   read (unit, '(a)', end=99) buffer
   mol%title = trim(buffer)

   do i = 1, num_atoms
      read (unit, *, end=99) tag, coords
      call split_tag(tag, elnum, label)
      mol%atoms(i)%elnum = elnum
      mol%atoms(i)%weight = atomic_weights(elnum)
      mol%atoms(i)%label = label
      mol%atoms(i)%coords = coords
   end do

   return

   99 continue
   write (stderr, '(a)') 'Unexpected end of file!'
   stop

end subroutine

subroutine readmol2(unit, mol)
   integer, intent(in) :: unit
   type(mol_type), intent(out) :: mol
   character(ll) :: buffer
   integer :: i, id, num_atoms, elnum, label
   integer :: nbond, atom1, atom2, bondorder
   character(wl) :: tag
   real(rk) :: coords(3)
   integer, allocatable :: nadjs(:), adjlists(:, :)

   do
      read (unit, '(a)', end=99) buffer
      if (buffer == '@<TRIPOS>MOLECULE') exit
   end do

   read (unit, '(a)', end=99) buffer
   mol%title = trim(buffer)
   read (unit, *, end=99) num_atoms, nbond

   allocate (mol%atoms(num_atoms))
   allocate (nadjs(num_atoms))
   allocate (adjlists(num_atoms, num_atoms))

   do
      read (unit, '(a)', end=99) buffer
      if (buffer == '@<TRIPOS>ATOM') exit
   end do

   do i = 1, num_atoms
      read (unit, *, end=99) id, tag, coords
      call split_tag(tag, elnum, label)
      mol%atoms(i)%elnum = elnum
      mol%atoms%weight = atomic_weights(elnum)
      mol%atoms(i)%label = label
      mol%atoms(i)%coords = coords
   end do

   do
      read (unit, '(a)', end=99) buffer
      if (buffer == '@<TRIPOS>BOND') exit
   end do

   ! Bond initialization
   nadjs(:) = 0
   adjlists(:, :) = 0

   do i = 1, nbond
      read (unit, *, end=99) id, atom1, atom2, bondorder
      nadjs(atom1) = nadjs(atom1) + 1
      nadjs(atom2) = nadjs(atom2) + 1
      adjlists(nadjs(atom1), atom1) = atom2
      adjlists(nadjs(atom2), atom2) = atom1
   end do

   call set_adjlists(mol, nadjs, adjlists)

   return

   99 continue
   write (stderr, '(a)') 'Unexpected end of file!'
   stop

end subroutine

subroutine set_bonds(mol)
   type(mol_type), target, intent(inout) :: mol
   ! Local variables
   integer :: i, j, num_atoms
   integer, allocatable :: nadjs(:), adjlists(:, :)
   type(atom_type), pointer :: atoms(:)
   real(rk), allocatable :: adjrads(:)
   real(rk) :: atomdist

   atoms => mol%atoms
   num_atoms = size(atoms)
   allocate (nadjs(num_atoms))
   allocate (adjlists(max_coord, num_atoms))

   ! Bond initialization
   nadjs(:) = 0

   ! Set adjacency radii
   adjrads = 0.75*covalent_radii(atoms%elnum) + 0.25*vdw_radii(atoms%elnum)

   ! Register adjacency matrix i,j if atoms i and j are closer
   ! than the sum of their adjacency radius
   do i = 1, num_atoms
      do j = i + 1, num_atoms
         atomdist = sqrt(sum((atoms(i)%coords - atoms(j)%coords)**2))
         if (atomdist < adjrads(i) + adjrads(j)) then
            nadjs(i) = nadjs(i) + 1
            nadjs(j) = nadjs(j) + 1
            if (nadjs(i) > max_coord .or. nadjs(j) > max_coord) then
               write (stderr, '("Maximum coordination number exceeded!")')
               stop
            end if
            adjlists(nadjs(i), i) = j
            adjlists(nadjs(j), j) = i
         end if
      end do
   end do

   call set_adjlists(mol, nadjs, adjlists)

end subroutine

end module
