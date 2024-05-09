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
   allocate(adjlist(mol%natom,mol%natom))

   nadj(:) = 0
   adjlist(:, :) = 0

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

   ! Return if bonds are not required
   if (.not. bond_flag) then
      mol%atoms(:)%nadj = 0
      return
   end if

   do i = 1, nbond
      read (unit, *, end=99) id, atom1, atom2, bondorder
      nadj(atom1) = nadj(atom1) + 1
      nadj(atom2) = nadj(atom2) + 1
      adjlist(nadj(atom1),atom1) = atom2
      adjlist(nadj(atom2),atom2) = atom1
   end do

   do i = 1, mol%natom
      mol%atoms(i)%nadj = nadj(i)
      allocate(mol%atoms(i)%adjlist(maxcoord))
!      allocate(mol%atoms(i)%adjlist(nadj(i)))
      mol%atoms(i)%adjlist(1:nadj(i)) = adjlist(1:nadj(i),i)
   end do

   return

   99 continue
   write (stderr, '(a)') 'Unexpected end of file!'
   stop

end subroutine

subroutine set_bonds (mol)
   type(Molecule), intent(inout) :: mol
   integer :: i, j
   real(wp) :: atomdist
   real(wp), dimension(:), allocatable :: adjrad
   real(wp), dimension(:, :), allocatable :: coords
   integer, dimension(:), allocatable :: znums
   integer, allocatable :: nadj(:), adjlist(:, :)

   ! Quick return if bonds are not required
   if (.not. bond_flag) then
      mol%atoms(:)%nadj = 0
      return
   end if

   allocate(adjrad(mol%natom))
   allocate(coords(3,mol%natom))
   allocate(znums(mol%natom))
   allocate(nadj(mol%natom),adjlist(maxcoord,mol%natom))
   znums = mol%get_znums()
   coords = mol%get_coords()

   ! initialization
   nadj(:) = 0

   ! Set adjacency radii
   adjrad(:) = covrad(znums) + 0.25*(vdwrad(znums) - covrad(znums))

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

   ! update adjecency lists in mol structure from adjecency matrix
   do i = 1, mol%natom
      mol%atoms(i)%nadj = nadj(i)
!      allocate(mol%atoms(i)%adjlist(maxcoord))   ! fixed array size
      allocate(mol%atoms(i)%adjlist(nadj(i)))   ! exact array size
      mol%atoms(i)%adjlist(1:nadj(i)) = adjlist(1:nadj(i),i)
   end do

end subroutine set_bonds

end module
