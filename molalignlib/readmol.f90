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
use pointers
use strutils
use chemdata
use chemutils

implicit none

contains

subroutine readxyz(unit, mol)
   integer, intent(in) :: unit
   type(Molecule), intent(out) :: mol
   character(ll) :: buffer
   character(wl) :: element
   integer, allocatable :: atomelnums(:), atomlabels(:)
   real(wp), allocatable :: weights(:)
   real(wp), allocatable :: atomcoords(:, :)

   integer :: i

   read (unit, *, end=99) mol%natom

   allocate(atomelnums(mol%natom))
   allocate(atomlabels(mol%natom))
   allocate(weights(mol%natom))
   allocate(atomcoords(3, mol%natom))
   allocate(mol%atoms(mol%natom))

   read (unit, '(a)', end=99) buffer
   mol%title = trim(buffer)

   do i = 1, mol%natom
      read (unit, *, end=99) element, atomcoords(:, i)
      call readlabel(element, atomelnums(i), atomlabels(i))
      weights(i) = weight_func(atomelnums(i))
   end do

   call mol%set_elnums(atomelnums)
   call mol%set_labels(atomlabels)
   call mol%set_weights(weights)
   call mol%set_coords(atomcoords)

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
   character(wl) :: element
   integer, allocatable :: atomelnums(:), atomlabels(:)
   real(wp), allocatable :: weights(:)
   real(wp), allocatable :: atomcoords(:, :)
   integer, allocatable :: nadjs(:), adjlists(:, :)

   do
      read (unit, '(a)', end=99) buffer
      if (buffer == '@<TRIPOS>MOLECULE') exit
   end do

   read (unit, '(a)', end=99) buffer
   mol%title = trim(buffer)
   read (unit, *, end=99) mol%natom, nbond

   allocate(atomelnums(mol%natom))
   allocate(atomlabels(mol%natom))
   allocate(weights(mol%natom))
   allocate(atomcoords(3, mol%natom))
   allocate(nadjs(mol%natom))
   allocate(adjlists(mol%natom, mol%natom))
   allocate(mol%atoms(mol%natom))

   do
      read (unit, '(a)', end=99) buffer
      if (buffer == '@<TRIPOS>ATOM') exit
   end do

   do i = 1, mol%natom
      read (unit, *, end=99) id, element, atomcoords(:, i)
      call readlabel(element, atomelnums(i), atomlabels(i))
      weights(i) = weight_func(atomelnums(i))
   end do

   do
      read (unit, '(a)', end=99) buffer
      if (buffer == '@<TRIPOS>BOND') exit
   end do

   ! Bond initialization
   nadjs(:) = 0
   adjlists(:, :) = 0

   ! Return if bonds are not requested
   if (.not. bond_flag) then
      call mol%set_adjlists(nadjs, adjlists)
      return
   end if

   do i = 1, nbond
      read (unit, *, end=99) id, atom1, atom2, bondorder
      nadjs(atom1) = nadjs(atom1) + 1
      nadjs(atom2) = nadjs(atom2) + 1
      adjlists(nadjs(atom1), atom1) = atom2
      adjlists(nadjs(atom2), atom2) = atom1
   end do

   call mol%set_elnums(atomelnums)
   call mol%set_labels(atomlabels)
   call mol%set_weights(weights)
   call mol%set_coords(atomcoords)
   call mol%set_adjlists(nadjs, adjlists)

   return

   99 continue
   write (stderr, '(a)') 'Unexpected end of file!'
   stop

end subroutine

subroutine set_bonds(mol)
   type(Molecule), intent(inout) :: mol

   integer :: i, j
   integer :: nadjs(mol%natom)
   integer :: adjlists(maxcoord, mol%natom)
   integer, allocatable :: atomelnums(:)
   real(wp) :: atomdist
   real(wp), allocatable :: adjrads(:)
   real(wp), allocatable :: atomcoords(:, :)

   ! Bond initialization
   nadjs(:) = 0

   ! Return if bonds are not requested
   if (.not. bond_flag) then
      call mol%set_adjlists(nadjs, adjlists)
      return
   end if

   atomelnums = mol%get_elnums()
   atomcoords = mol%get_coords()

   ! Set adjacency radii
   adjrads = covrad(atomelnums) + 0.25*(vdwrad(atomelnums) - covrad(atomelnums))

   ! Register adjacency matrix i,j if atoms i and j are closer
   ! than the sum of their adjacency radius
   do i = 1, mol%natom
      do j = i + 1, mol%natom
         atomdist = sqrt(sum((atomcoords(:, i) - atomcoords(:, j))**2))
         if (atomdist < adjrads(i) + adjrads(j)) then
            nadjs(i) = nadjs(i) + 1
            nadjs(j) = nadjs(j) + 1
            if (nadjs(i) > maxcoord .or. nadjs(j) > maxcoord) then
               write (stderr, '("Maximum coordination number exceeded!")')
               stop
            end if
            adjlists(nadjs(i), i) = j
            adjlists(nadjs(j), j) = i
         end if
      end do
   end do

   call mol%set_adjlists(nadjs, adjlists)

end subroutine

end module
