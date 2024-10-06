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

module molalignlib
use stdio
use kinds
use molecule
use bounds
use random
use linalg
use sorting
use permutation
use rotation
use translation
use assorting
use adjacency
use alignment
use remapping
use reactivity
use assignment
use writemol
use biasing
use tracking

implicit none

contains

! Assign atoms0 and atoms1
subroutine molecule_remap( &
   mol0, &
   mol1, &
   nrec, &
   maplist, &
   countlist)

   type(molecule_type), intent(inout) :: mol0, mol1
   integer, intent(out) :: nrec
   integer, dimension(:, :), intent(inout) :: maplist
   integer, dimension(:), intent(inout) :: countlist

   real(rk) :: travec0(3), travec1(3)
   real(rk), allocatable, dimension(:, :) :: coords0, coords1

!   integer, dimension(mol0%natom) :: mapping
!   integer :: nbond0, bonds0(2, maxcoord*mol0%natom)
!   integer :: nbond1, bonds1(2, maxcoord*mol1%natom)

   call mol0%set_atompartition(mol0%mnatypes)
   call mol1%set_atompartition(mol1%mnatypes)

   ! Abort if molecules have different number of atoms

   if (mol0%get_natom() /= mol1%get_natom()) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if molecules are not isomers

   if (any(mol0%gather_elnums() /= mol1%gather_elnums())) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if there are conflicting atomic types

   if (any(mol0%gather_atomeltypes() /= mol1%gather_atomeltypes())) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      stop
   end if

   ! Backup coordinates

   coords0 = mol0%get_coords()
   coords1 = mol1%get_coords()

   ! Mirror coordinates

   if (mirror_flag) then
      call mol1%mirror_coords()
   end if

   ! Calculate centroids

   travec0 = -centroid(mol0)
   travec1 = -centroid(mol1)

   ! Center coordinates at the centroids

   call mol0%translate_coords(travec0)
   call mol1%translate_coords(travec1)

   ! Initialize random number generator

   call random_initialize()

   ! Optimize assignment to minimize the AdjD and RMSD

   call remap_atoms(mol0, mol1, maplist, countlist, nrec)

   ! Remove bonds from reactive sites and reoptimize assignment

   if (reac_flag) then
!      call mol0%print_bonds()
!      call mol1%print_bonds()
      call remove_reactive_bonds(mol0, mol1, maplist(:, 1))
      call find_molfrags(mol0)
      call find_molfrags(mol1)
      call assort_neighbors(mol0)
      call assort_neighbors(mol1)
      call remap_atoms(mol0, mol1, maplist, countlist, nrec)
   end if

   ! Restore coordinates

   call mol0%set_coords(coords0)
   call mol1%set_coords(coords1)

!   ! Print coordinates with internal order

!   call rotate(natom1, workcoords1, leastrotquat(natom0, workweights0, workcoords0, workcoords1, mapping))
!   open(unit=99, file='aligned_debug.mol2', action='write', status='replace')
!   call adjmat2bonds(natom0, workadjmat0, nbond0, bonds0)
!   call adjmat2bonds(natom1, workadjmat1, nbond1, bonds1)
!!   call adjmat2bonds(natom1, workadjmat1(mapping, mapping), nbond1, bonds1)
!   call writemol2(99, 'coords0', natom0, workznums0, workcoords0, nbond0, bonds0)
!   call writemol2(99, 'coords1', natom1, workznums1, workcoords1, nbond1, bonds1)
!!   call writemol2(99, 'coords1', natom1, workznums1(mapping), workcoords1(:, mapping), nbond1, bonds1)

end subroutine

! Align atoms
subroutine remapped_molecule_align( &
   mol0, &
   mol1, &
   mapping, &
   travec0, &
   travec1, &
   rotquat)

   type(molecule_type), intent(in) :: mol0, mol1
   integer, intent(in) :: mapping(:)
   real(rk), intent(out) :: travec0(3), travec1(3), rotquat(4)
   ! Local variables
   real(rk), allocatable :: weights(:)
   real(rk), allocatable, dimension(:, :) :: coords0, coords1

   weights = mol0%get_weights()
   coords0 = mol0%get_coords()
   coords1 = mol1%get_coords()

   ! Calculate centroids

   travec0 = -centroid(mol0)
   travec1 = -centroid(mol1)

   ! Calculate optimal rotation matrix

   rotquat = leastrotquat( &
      mol0%natom, &
      weights, &
      translated(mol0%natom, coords0, travec0), &
      translated(mol1%natom, coords1, travec1), &
      mapping &
   )

end subroutine

subroutine molecule_align( &
! Purpose: Align atoms0 and atoms1
   mol0, &
   mol1, &
   travec0, &
   travec1, &
   rotquat)

   type(molecule_type), intent(in) :: mol0, mol1
   real(rk), intent(out) :: travec0(3), travec1(3), rotquat(4)

   ! Abort if molecules have different number of atoms

   if (mol0%get_natom() /= mol1%get_natom()) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if molecules are not isomers

   if (any(mol0%gather_elnums() /= mol1%gather_elnums())) then
      write (stderr, '(a)') '*Error: These molecules are not isomers'
      stop
   end if

   ! Abort if there are conflicting atomic types

   if (any(mol0%gather_atomeltypes() /= mol1%gather_atomeltypes())) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      stop
   end if

   ! Abort if atoms are not ordered

   if (any(mol0%get_elnums() /= mol1%get_elnums())) then
      write (stderr, '(a)') 'Error: The atoms are not in the same order'
      stop
   end if

   ! Abort if atomic types are not ordered

   if (any(mol0%get_atomeltypes() /= mol1%get_atomeltypes())) then
      write (stderr, '(a)') 'Error: Atomic types are not in the same order'
      stop
   end if

   ! Calculate centroids

   travec0 = -centroid(mol0)
   travec1 = -centroid(mol1)

   ! Calculate optimal rotation matrix

   rotquat = leastrotquat( &
      mol0%natom, &
      mol0%get_weights(), &
      translated(mol0%natom, mol0%get_coords(), travec0), &
      translated(mol1%natom, mol1%get_coords(), travec1), &
      identity_permutation(mol0%natom) &
   )

end subroutine

function get_rmsd(mol0, mol1, mapping) result(rmsd)
   type(molecule_type), intent(in) :: mol0, mol1
   integer :: mapping(:)
   real(rk) :: rmsd

   rmsd = sqrt(squaredist(mol0%natom, mol0%get_weights(), mol0%get_coords(), &
         mol1%get_coords(), mapping) / sum(mol0%get_weights()))

end function

function get_adjd(mol0, mol1, mapping) result(adjd)
   type(molecule_type), intent(in) :: mol0, mol1
   integer :: mapping(:)
   integer :: adjd

   adjd = adjacencydiff(mol0%natom, mol0%get_adjmatrix(), mol1%get_adjmatrix(), mapping)

end function

function centroid(mol)
! Purpose: Get the centroid coordinates
   type(molecule_type), intent(in) :: mol
   ! Local variables
   integer :: i
   real(rk) :: centroid(3)
   real(rk), allocatable :: coords(:, :)
   real(rk), allocatable :: weights(:)

   coords = mol%get_coords()
   weights = mol%get_weights()

! Calculate the coordinates of the center of mass

   centroid(:) = 0

   do i = 1, size(mol%atoms)
      centroid(:) = centroid(:) + weights(i)*coords(:, i)
   end do

   centroid(:) = centroid(:)/sum(weights)

end function

end module
