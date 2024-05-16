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
use types
use flags
use bounds
use strutils
use chemdata

implicit none

contains

subroutine readxyz(unit, mol)
   integer, intent(in) :: unit
   type(Molecule), intent(out) :: mol
   character(ll) :: buffer

   integer :: i

   read (unit, *, end=99) mol%natom

   allocate(mol%atoms(mol%natom))

   read (unit, '(a)', end=99) buffer
   mol%title = trim(buffer)

   do i = 1, mol%natom
      read (unit, *, end=99) buffer, mol%atoms(i)%coords(:)
      mol%atoms(i)%label = buffer
   end do

   return

   99 continue
   write (stderr, '(a)') 'Unexpected end of file!'
   stop

end subroutine

subroutine readmol2(unit, mol)
   integer, intent(in) :: unit
   type(Molecule), intent(out) :: mol
   character(ll) :: buffer
   integer :: i, id
   integer :: atom1, atom2, bondorder, nbond
   integer, allocatable :: nadj(:), adjlist(:, :)

   do
      read (unit, '(a)', end=99) buffer
      if (buffer == '@<TRIPOS>MOLECULE') exit
   end do

   read (unit, '(a)', end=99) buffer
   mol%title = trim(buffer)
   read (unit, *, end=99) mol%natom, nbond

   allocate(mol%atoms(mol%natom))
   allocate(nadj(mol%natom))
   allocate(adjlist(mol%natom, mol%natom))

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

   ! Bond initialization
   nadj(:) = 0
   adjlist(:, :) = 0

   ! Return if bonds are not requested
   if (.not. bond_flag) then
      call mol%set_adjlist(nadj, adjlist)
      return
   end if

   do i = 1, nbond
      read (unit, *, end=99) id, atom1, atom2, bondorder
      nadj(atom1) = nadj(atom1) + 1
      nadj(atom2) = nadj(atom2) + 1
      adjlist(nadj(atom1), atom1) = atom2
      adjlist(nadj(atom2), atom2) = atom1
   end do

   call mol%set_adjlist(nadj, adjlist)

   return

   99 continue
   write (stderr, '(a)') 'Unexpected end of file!'
   stop

end subroutine

subroutine set_bonds(mol)
   type(Molecule), intent(inout) :: mol

   integer :: i, j
   integer, allocatable :: znums(:)
   integer :: nadj(mol%natom), adjlist(maxcoord, mol%natom)
   real(wp) :: atomdist
   real(wp), allocatable :: adjrad(:)
   real(wp), allocatable :: coords(:, :)

   ! Bond initialization
   nadj(:) = 0

   ! Return if bonds are not requested
   if (.not. bond_flag) then
      call mol%set_adjlist(nadj, adjlist)
      return
   end if

   znums = mol%get_znums()
   coords = mol%get_coords()

   ! Set adjacency radii
   adjrad = covrad(znums) + 0.25*(vdwrad(znums) - covrad(znums))

   ! Register adjacency matrix i,j if atoms i and j are closer
   ! than the sum of their adjacency radius
   do i = 1, mol%natom
      do j = i + 1, mol%natom
         atomdist = sqrt(sum((coords(:, i) - coords(:, j))**2))
         if (atomdist < adjrad(i) + adjrad(j)) then
            nadj(i) = nadj(i) + 1
            nadj(j) = nadj(j) + 1
            if (nadj(i) > maxcoord .or. nadj(j) > maxcoord) then
               write (stderr, '("Maximum coordination number exceeded!")')
               stop
            end if
            adjlist(nadj(i), i) = j
            adjlist(nadj(j), j) = i
         end if
      end do
   end do

   call mol%set_adjlist(nadj, adjlist)

end subroutine set_bonds

end module
