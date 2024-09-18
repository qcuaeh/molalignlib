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
use types
use bounds
use random
use sorting
use discrete
use rotation
use translation
use assorting
use adjacency
use alignment
use remapping
use assignment
use writemol
use biasing
use tracking

implicit none

contains

! Assign atoms0 and atoms1
subroutine remap_atoms( &
   mol0, &
   mol1, &
   nrec, &
   maplist, &
   countlist)

   type(cMol), intent(inout) :: mol0, mol1
   integer, intent(out) :: nrec
   integer, dimension(:, :), intent(inout) :: maplist
   integer, dimension(:), intent(inout) :: countlist

   integer :: i
   integer, allocatable, dimension(:) :: mnaord0, mnaord1
   real(wp) :: travec0(3), travec1(3)
   real(wp), allocatable, dimension(:, :) :: atomcoords0, atomcoords1
   type(cAtom), allocatable, dimension(:) :: atoms0, atoms1

!   integer, dimension(mol0%natom) :: mapping
!   integer :: nbond0, bonds0(2, maxcoord*mol0%natom)
!   integer :: nbond1, bonds1(2, maxcoord*mol1%natom)

   !  Select assignment algorithm

   if (bond_flag) then
      mapatoms => mapatoms_bonded
      mapsetcrossbias => setcrossbias_mna
   else
      mapatoms => mapatoms_free
      mapsetcrossbias => setcrossbias_rd
   end if

   ! Select bias function

   if (bias_flag) then
      setcrossbias => mapsetcrossbias
   else
      setcrossbias => setcrossbias_none
   end if

   ! Set MNA equivalence partitions

   mnaord0 = sorted_order(mol0%get_atommnatypes())
   mnaord1 = sorted_order(mol1%get_atommnatypes())

   ! Get sorted atoms lists

   atoms0 = mol0%get_sorted_atoms()
   atoms1 = mol1%get_sorted_atoms()

   ! Abort if molecules have different number of atoms

   if (size(atoms0) /= size(atoms1)) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if molecules are not isomers

   if (any(atoms0%elnum /= atoms1%elnum)) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if there are conflicting atomic types

   if (any(atoms0%eltype /= atoms1%eltype)) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      stop
   end if

   ! Abort if there are conflicting weights

   if (any(abs(atoms0%weight - atoms1%weight) > 1.E-6)) then
      write (stderr, '(a)') 'Error: There are conflicting weights'
      stop
   end if

   ! Backup coordinates

   atomcoords0 = mol0%get_atomcoords()
   atomcoords1 = mol1%get_atomcoords()

   ! Mirror coordinates

   if (mirror_flag) then
      call mol1%mirror_atomcoords()
   end if

   ! Calculate centroids

   travec0 = -centroid_coords(mol0)
   travec1 = -centroid_coords(mol1)

   ! Center coordinates at the centroids

   call mol0%translate_atomcoords(travec0)
   call mol1%translate_atomcoords(travec1)

   ! Initialize random number generator

   call initialize_random()

   ! Optimize assignment to minimize the AdjD and RMSD

   call optimize_mapping(mol0, mol1, mnaord0, mnaord1, maplist, countlist, nrec)

   ! Remove bonds from reactive sites and reoptimize assignment

   if (reac_flag) then
!      call mol0%print_bonds()
!      call mol1%print_bonds()
      call remove_reactive_bonds(mol0, mol1, maplist(:, 1))
      call find_molfrags(mol0)
      call find_molfrags(mol1)
      call assort_neighbors(mol0)
      call assort_neighbors(mol1)
      call optimize_mapping(mol0, mol1, mnaord0, mnaord1, maplist, countlist, nrec)
   end if

   ! Restore coordinates

   call mol0%set_atomcoords(atomcoords0)
   call mol1%set_atomcoords(atomcoords1)

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
subroutine align_remapped_atoms( &
   mol0, &
   mol1, &
   mapping, &
   travec, &
   rotmat)

   type(cMol), intent(inout) :: mol0, mol1
   integer :: mapping(:)
   real(wp), intent(out) :: travec(3), rotmat(3, 3)
   ! Local variables
   real(wp) :: travec0(3), travec1(3), rotquat(4)
   real(wp), allocatable :: atomweights(:)
   real(wp), allocatable, dimension(:, :) :: atomcoords0, atomcoords1

   atomweights = mol0%get_atomweights()
   atomcoords0 = mol0%get_atomcoords()
   atomcoords1 = mol1%get_atomcoords()

   ! Calculate centroids

   travec0 = -centroid_coords(mol0)
   travec1 = -centroid_coords(mol1)

   ! Calculate optimal rotation matrix

   rotquat = leastrotquat( &
      mol0%natom, &
      atomweights, &
      translated(mol0%natom, atomcoords0, travec0), &
      translated(mol1%natom, atomcoords1, travec1), &
      mapping &
   )

   rotmat = quat2rotmat(rotquat)

   ! Calculate optimal translation vector

   travec = matmul(rotmat, travec1) - travec0

end subroutine

subroutine align_atoms( &
! Purpose: Align atoms0 and atoms1
   mol0, &
   mol1, &
   travec, &
   rotmat)

   type(cMol), intent(inout) :: mol0, mol1
   real(wp), intent(out) :: travec(3), rotmat(3, 3)
   real(wp) :: travec0(3), travec1(3), rotquat(4)
   type(cAtom), allocatable, dimension(:) :: atoms0, atoms1

   ! Get sorted atoms lists

   atoms0 = mol0%get_sorted_atoms()
   atoms1 = mol1%get_sorted_atoms()

   ! Abort if molecules have different number of atoms

   if (size(atoms0) /= size(atoms1)) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if molecules are not isomers

   if (any(atoms0%elnum /= atoms1%elnum)) then
      write (stderr, '(a)') '*Error: These molecules are not isomers'
      stop
   end if

   ! Abort if there are conflicting atomic types

   if (any(atoms0%eltype /= atoms1%eltype)) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      stop
   end if

   ! Get unsorted atoms lists

   atoms0 = mol0%get_atoms()
   atoms1 = mol1%get_atoms()

   ! Abort if atoms are not ordered

   if (any(atoms0%elnum /= atoms1%elnum)) then
      write (stderr, '(a)') 'Error: The atoms are not in the same order'
      stop
   end if

   ! Abort if atomic types are not ordered

   if (any(atoms0%eltype /= atoms1%eltype)) then
      write (stderr, '(a)') 'Error: Atomic types are not in the same order'
      stop
   end if

   ! Calculate centroids

   travec0 = -centroid_coords(mol0)
   travec1 = -centroid_coords(mol1)

   ! Calculate optimal rotation matrix

   rotquat = leastrotquat( &
      mol0%natom, &
      mol0%get_atomweights(), &
      translated(mol0%natom, mol0%get_atomcoords(), travec0), &
      translated(mol1%natom, mol1%get_atomcoords(), travec1), &
      identitymap(mol0%natom) &
   )

   rotmat = quat2rotmat(rotquat)

   ! Calculate optimal translation vector

   travec = matmul(rotmat, travec1) - travec0

end subroutine

function get_rmsd(mol0, mol1) result(rmsd)
   type(cMol), intent(in) :: mol0, mol1
   real(wp) :: rmsd

   rmsd = sqrt(squaredist(mol0%natom, mol0%get_atomweights(), mol0%get_atomcoords(), &
         mol1%get_atomcoords(), identitymap(mol0%natom)) / sum(mol0%get_atomweights()))

end function

function get_adjd(mol0, mol1) result(adjd)
   type(cMol), intent(in) :: mol0, mol1
   integer :: adjd

   adjd = adjacencydiff(mol0%natom, mol0%get_adjmat(), mol1%get_adjmat(), identitymap(mol0%natom))

end function

function centroid_coords(mol) result(coords)
! Purpose: Get the centroid coordinates
   type(cMol), intent(in) :: mol
   ! Local variables
   integer :: i
   real(wp) :: coords(3)
   real(wp), allocatable :: atomcoords(:, :)
   real(wp), allocatable :: atomweights(:)

   atomcoords = mol%get_atomcoords()
   atomweights = mol%get_atomweights()

! Calculate the coordinates of the center of mass

   coords(:) = 0

   do i = 1, size(mol%atoms)
      coords(:) = coords(:) + atomweights(i)*atomcoords(:, i)
   end do

   coords(:) = coords(:)/sum(atomweights)

end function

end module
