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

module alignment_conformer
use parameters
use common_types
use str_utils
use random
use chemdata
use permutation
use euclidean
use adjacency
use pruning_atoms
use linked_list_types
use indexed_list_types
use refinement
use assignment_atoms
use assignment_conformer
use recording
use flags
use error_codes
implicit none

contains

subroutine optimize_mapping_conformer( adjcs1, adjcs2, atomtypes, &
      coords1, coords2, conv_freq, max_trials, registry, error_code)
! adjcs, atomtypes and coords hold the included atoms only, in a common
! numbering; the atom permutations stored in registry use that numbering.
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(partition_t), intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer(ik), intent(in) :: conv_freq, max_trials
   type(registry_t), target, intent(inout) :: registry
   integer(ik), intent(out) :: error_code

   ! Local variables
   integer(ik), dimension(:), allocatable :: mapping1, new_mapping
   real(rk), dimension(:,:), allocatable :: coords2r
   real(rk), dimension(4) :: rotation, total_rotation
   real(rk) :: steps, mapdist, new_mapdist
   type(chaintree_node_t), pointer :: hna_chain
   type(array_trees_t) :: cache_arrays

   ! Allocations
   allocate (coords2r, mold=coords2)

   ! Pre-compute assignment tree for decision making
   call compute_sc_hna_chain( adjcs1, adjcs2, atomtypes, hna_chain, error_code)
   if (error_code /= MOLALIGN_SUCCESS) return
   call build_assignment_tree( adjcs1, adjcs2, hna_chain%last_link, cache_arrays, error_code)
   if (error_code /= MOLALIGN_SUCCESS) return

   if (print_assigntree) then
      call print_chain_tree_array( atomtypes, cache_arrays)
   end if

   ! Reset registry for new conformer
   call reset_registry( registry)

   ! Choose the search strategy (see FORCE_EXHAUSTIVE and FORCE_STOCHASTIC)
   if (.not. FORCE_EXHAUSTIVE .and. (FORCE_STOCHASTIC .or. &
         cache_arrays%total_combinations > conv_freq*cache_arrays%partial_combinations)) then

      ! Initialize random number generator
      call random_initialize()

      ! Optimize atom permutation
      do

         ! Apply random rotation to coords2 copy
         coords2r = coords2
         total_rotation = randrotquat()
         call rotate_coords( coords2r, total_rotation)

         ! Assign atoms with current orientation
         if (PRUNE_ASSIGNMENTS) then
            call assign_atoms_greedy( coords1, coords2r, cache_arrays, mapping1, mapdist)
            call assign_atoms_local_pruned( coords1, coords2r, cache_arrays, mapping1, mapdist)
         else
            call assign_atoms_local( coords1, coords2r, cache_arrays, mapping1, mapdist)
         end if

         rotation = least_rotquat( mapping1, coords1, coords2r)
         total_rotation = quatmul( total_rotation, rotation)
         call rotate_coords( coords2r, rotation)
         mapdist = sqdistsum( mapping1, coords1, coords2r)
         steps = 1

         do
            if (PRUNE_ASSIGNMENTS) then
               new_mapdist = mapdist
               call assign_atoms_local_pruned( coords1, coords2r, cache_arrays, new_mapping, new_mapdist)
            else
               call assign_atoms_local( coords1, coords2r, cache_arrays, new_mapping, new_mapdist)
            end if
!            write (stdout,*) mapdist, new_mapdist
            if (all(mapping1 == new_mapping)) exit
            mapping1 = new_mapping
            rotation = least_rotquat( mapping1, coords1, coords2r)
            total_rotation = quatmul( total_rotation, rotation)
            call rotate_coords( coords2r, rotation)
            mapdist = sqdistsum( mapping1, coords1, coords2r)
            steps = steps + 1
         end do

         ! Update results
         call insert_record_mapping( registry, mapping1, steps, total_rotation, 0, mapdist)

         if (registry%records(1)%freq > conv_freq) then
            exit
         end if

         if (MAX_TRIALS_EXIT) then
            if (registry%n_trials > max_trials) then
               exit
            end if
         end if
      end do

   else

      ! Assign atoms using global assignment
      call assign_atoms_global( coords1, coords2, cache_arrays, mapping1)
      
      ! Calculate optimal rotation
      rotation = least_rotquat( mapping1, coords1, coords2)
      
      ! Rotate coords2 and calculate mapdist
      coords2r = coords2
      call rotate_coords( coords2r, rotation)
      mapdist = sqdistsum( mapping1, coords1, coords2r)
      
      ! Initialize registry with a single record
      registry%occ_records = 1
      registry%records(1)%mapping1 = mapping1
      registry%records(1)%mapdiff = 0
      registry%records(1)%mapdist = mapdist
      registry%records(1)%freq = 1
      registry%records(1)%steps = 1
      registry%records(1)%rotation = rotation

   end if
end subroutine

subroutine assign_mapping_conformer( adjcs1, adjcs2, atomtypes, coords1, coords2, mapping1, error_code)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(partition_t), intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer(ik), dimension(:), allocatable, intent(out) :: mapping1
   integer(ik), intent(out) :: error_code

   ! Local variables
   type(chaintree_node_t), pointer :: hna_chain
   type(array_trees_t) :: cache_arrays
   real(rk) :: mapdist

   ! Pre-compute assignment tree
   call compute_sc_hna_chain( adjcs1, adjcs2, atomtypes, hna_chain, error_code)
   if (error_code /= MOLALIGN_SUCCESS) return
   call build_assignment_tree( adjcs1, adjcs2, hna_chain%last_link, cache_arrays, error_code)
   if (error_code /= MOLALIGN_SUCCESS) return

   if (print_assigntree) then
      call print_chain_tree_array( atomtypes, cache_arrays)
   end if

   call assign_atoms_local( coords1, coords2, cache_arrays, mapping1, mapdist)
!   call assign_atoms_greedy( coords1, coords2, cache_arrays, mapping1, mapdist)
!   call assign_atoms_local_pruned( coords1, coords2, cache_arrays, mapping1, mapdist)
end subroutine

end module
