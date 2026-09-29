! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

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

module alignment_atoms
use parameters
use random
use chemdata
use permutation
use euclidean
use assignment_atoms
use linked_list_types
use refinement
use pruning_atoms
use recording
use flags
use error_codes
implicit none

contains

subroutine optimize_mapping_atoms(atomtypes, prunes, coords1, &
      coords2, conv_freq, max_trials, registry, error_code)
! coords1 and coords2 hold the included atoms only, in the numbering of
! atomtypes; the atom permutations stored in registry use that numbering.
   type(partition_t), intent(in) :: atomtypes
   type(bool_matrix), dimension(:), intent(in) :: prunes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer(ik), intent(in) :: conv_freq, max_trials
   type(registry_t), intent(inout) :: registry
   integer(ik), intent(out) :: error_code

   ! Local variables
   integer(ik), dimension(:), allocatable :: mapping1, new_mapping
   real(rk), dimension(:,:), allocatable :: coords2r
   real(rk) :: mapdist, steps, rotation(4), total_rotation(4)

   ! Allocations
   allocate (coords2r, mold=coords2)
   allocate (mapping1(size(atomtypes%itemdir1)))
   allocate (new_mapping(size(atomtypes%itemdir1)))

   error_code = MOLALIGN_SUCCESS

   ! Initialize local minima registry
   call reset_registry( registry)

   ! Initialize random number generator
   call random_initialize()

   ! Optimize atom permutation
   do while (registry%records(1)%freq < conv_freq .and. registry%n_trials < max_trials)

      ! Apply random rotation to coords2 copy
      coords2r = coords2
      total_rotation = randrotquat()
      call rotate_coords( coords2r, total_rotation)

      ! Assign atoms with current orientation
      call assign_atoms_pruned( atomtypes, coords1, coords2r, prunes, mapping1, error_code)
      if (error_code /= MOLALIGN_SUCCESS) return
      rotation = least_rotquat( mapping1, coords1, coords2r)
      call rotate_coords( coords2r, rotation)
      total_rotation = quatmul( total_rotation, rotation)
      steps = 1

      do
         call assign_atoms_pruned( atomtypes, coords1, coords2r, prunes, new_mapping, error_code)
         if (error_code /= MOLALIGN_SUCCESS) return
         if (all(new_mapping == mapping1)) exit
         mapping1 = new_mapping
         rotation = least_rotquat( mapping1, coords1, coords2r)
         call rotate_coords( coords2r, rotation)
         total_rotation = quatmul( total_rotation, rotation)
         steps = steps + 1
      end do

      ! Push local minimum to registry
      mapdist = sqrt( sqdistsum( mapping1, coords1, coords2r))
      call insert_record_mapping( registry, mapping1, steps, total_rotation, 0, mapdist)

   end do

end subroutine

end module
