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
use parameters
use globals
use linalg
use sorting
use molecule
use rotation
use rigid_body
use permutation
use adjacency
use alignment
use atom_mapping_reac
use atom_mapping_conf
use writemol
use lcrs_tree
use partitioning

implicit none

contains

! Assign atoms0 and atoms1
subroutine molecule_remap( mol1, mol2, results)
!
! Find best atom mapping
!

   ! Arguments
   type(mol_type), intent(inout) :: mol1, mol2
   type(registry_type), intent(out) :: results

   ! Local variables
   type(tree_node), pointer :: eltypes_tree
   type(bipartition_container) :: eltypes

   ! Abort if molecules have different number of atoms
   if (size(mol1%atoms) /= size(mol2%atoms)) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if molecules are not isomers
   if (any(sorted(mol1%atoms%elnum) /= sorted(mol2%atoms%elnum))) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Compute atomic types
   call compute_eltypes(mol1, mol2, eltypes_tree)
   call partition_from_tree(eltypes_tree, eltypes)

   ! Abort if there are conflicting atomic types
   if (any(sorted(eltypes%itemdir1) /= sorted(eltypes%itemdir2))) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      stop
   end if

   ! Optimize assignment to minimize the AdjD and RMSD
   call remap_reactive_bonds( mol1, mol2, eltypes, results)

   call writemol2(stdout, mol1)
   call writemol2(stdout, mol2)

   ! Optimize assignment to minimize the AdjD and RMSD
   call remap_conformations( mol1, mol2, eltypes, results)

end subroutine

subroutine molecule_align( &
! Purpose: Align atoms0 and atoms1
   mol1, &
   mol2, &
   travec1, &
   travec2, &
   rotquat)

   type(mol_type), intent(in) :: mol1, mol2
   real(rk), intent(out) :: travec1(3), travec2(3), rotquat(4)
   type(tree_node), pointer :: eltypes_tree
   type(bipartition_container) :: eltypes
   ! Local variables
   integer :: num_atoms1
   real(rk), allocatable, dimension(:, :) :: coords1, coords2

   ! Abort if molecules have different number of atoms
   if (size(mol1%atoms) /= size(mol2%atoms)) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if molecules are not isomers
   if (any(sorted(mol1%atoms%elnum) /= sorted(mol2%atoms%elnum))) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Compute atomic types
   call compute_eltypes(mol1, mol2, eltypes_tree)
   call partition_from_tree(eltypes_tree, eltypes)
   call print_partition(eltypes)

   ! Abort if there are conflicting atomic types
   if (any(sorted(eltypes%itemdir1) /= sorted(eltypes%itemdir2))) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      stop
   end if

   ! Abort if atoms are not ordered
   if (any(mol1%atoms%elnum /= mol2%atoms%elnum)) then
      write (stderr, '(a)') 'Error: The atoms are not in the same order'
      stop
   end if

   ! Abort if atomic types are not ordered
   if (any(eltypes%itemdir1 /= eltypes%itemdir2)) then
      write (stderr, '(a)') 'Error: Atomic types are not in the same order'
      stop
   end if

   num_atoms1 = size(mol1%atoms)
   coords1 = mol1%get_weightcoords()
   coords2 = mol2%get_weightcoords()

   ! Calculate centroids
   travec1 = -centroid(coords1)
   travec2 = -centroid(coords2)

   ! Calculate optimal rotation matrix
   rotquat = leasteigquat( &
      identity_perm(num_atoms1), &
      translated_coords(coords1, travec1), &
      translated_coords(coords2, travec2) &
   )

end subroutine

end module
