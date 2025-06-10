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
use sorting
use molecule
use spatial
use permutation
use adjacency
use lcrs_tree
use eltype_compute
!use remapping_bonded
!use writemol

implicit none

contains

subroutine molecule_align( mol1, mol2, coords2)

   type(mol_type), intent(in) :: mol1, mol2
   real(rk), dimension(:,:), allocatable, intent(out) :: coords2
   ! Local variables
   type(partition_t) :: eltypes
   real(rk) :: center1(3), center2(3), rotquat(4)
   real(rk), dimension(:,:), allocatable :: coords1
   real(rk), dimension(:), allocatable :: weights1, weights2
   integer :: num_atoms1

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
   call set_eltypes( mol1%atoms, mol2%atoms, eltypes)

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
   coords1 = get_coords( mol1)
   weights1 = atomic_weights(mol1%atoms%elnum)
   call weight_coords( coords1, weights1)
   coords2 = get_coords( mol2)
   weights2 = atomic_weights(mol2%atoms%elnum)
   call weight_coords( coords2, weights2)

   ! Calculate centroids
   center1 = centroid( coords1)
   center2 = centroid( coords2)

   call translate_coords( coords2, -center2)
   call translate_coords( coords2, center1)

   ! Calculate optimal rotation matrix
   rotquat = optimal_rotation( &
      identity_perm(num_atoms1), &
      coords1, &
      coords2, &
      center1 &
   )

   call rotate_coords( coords2, rotquat, center1)

end subroutine

end module
