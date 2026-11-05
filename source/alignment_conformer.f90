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
      coords1, coords2, min_freq, max_trials, registry, error_code)
! Symmetry-corrected alignment of two conformers (Figure S1 of the paper).
! The assignment tree is built once and the strategy is chosen from its
! shape: when the total number of assignments exceeds min_freq times the
! sum of partial combinations, the stochastic search is used; otherwise
! every assignment is enumerated and scored after optimal superposition.
!
! Stochastic search: coords2 is given a random orientation, the optimal
! assignment at that orientation is found and the pair is iterated
! (superposition, reassignment) until the assignment no longer changes.
! Each local minimum is recorded, and the search stops when the best one
! has been found min_freq times or after max_trials trials (the same
! stopping rules as alignment_atoms::optimize_mapping_atoms).
!
! In both strategies registry holds up to size(registry%records) distinct
! ranked assignments (the -maxrecs option): the distinct local minima of the
! stochastic search, or the best complete assignments of the exhaustive one.
!
! adjcs, atomtypes and coords hold the included atoms only, in a common
! numbering; the atom permutations stored in registry use that numbering.
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(partition_t), intent(in) :: atomtypes
   real(RK), dimension(:,:), intent(in) :: coords1, coords2
   integer(IK), intent(in) :: min_freq, max_trials
   type(registry_t), target, intent(inout) :: registry
   integer(IK), intent(out) :: error_code

   ! Local variables
   integer(IK) :: i
   integer(IK), dimension(:), allocatable :: mapping1, new_mapping
   real(RK), dimension(:,:), allocatable :: coords2r
   real(RK), dimension(4) :: rotation, total_rotation
   real(RK) :: steps, mapdist, new_mapdist
   type(chaintree_node_t), pointer :: hna_chain
   type(array_trees_t) :: cache_arrays

   allocate (coords2r(size(coords2, 1), size(coords2, 2)))

   ! Self-consistent HNA partition and assignment tree
   call compute_sc_hna_chain( adjcs1, adjcs2, atomtypes, hna_chain, error_code)
   if (error_code /= MOLALIGN_SUCCESS) return
   call build_assignment_tree( adjcs1, adjcs2, hna_chain%last_link, cache_arrays, error_code)
   if (error_code /= MOLALIGN_SUCCESS) return

   if (printassigntree_flag) then
      call print_chain_tree_array( atomtypes, cache_arrays)
   end if

   call reset_registry( registry)

   ! Adaptive choice of strategy (see FORCE_EXHAUSTIVE and FORCE_STOCHASTIC)
   if (.not. FORCE_EXHAUSTIVE .and. (FORCE_STOCHASTIC .or. &
         cache_arrays%total_combinations > min_freq*cache_arrays%partial_combinations)) then

      call random_initialize()

      ! Search until the best local minimum has been found min_freq times
      do while (registry%records(1)%freq < min_freq)

         ! Random orientation of molecule 2
         coords2r = coords2
         total_rotation = randrotquat()
         call rotate_coords( coords2r, total_rotation)

         ! Optimal assignment at this orientation. With pruning, the greedy
         ! assignment provides the initial distance budget.
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

         ! Alternate superposition and assignment until the assignment is
         ! stable. The distance of the current assignment after
         ! superposition bounds the pruned search.
         do
            if (PRUNE_ASSIGNMENTS) then
               new_mapdist = mapdist
               call assign_atoms_local_pruned( coords1, coords2r, cache_arrays, new_mapping, new_mapdist)
            else
               call assign_atoms_local( coords1, coords2r, cache_arrays, new_mapping, new_mapdist)
            end if
            if (all(mapping1 == new_mapping)) exit
            mapping1 = new_mapping
            rotation = least_rotquat( mapping1, coords1, coords2r)
            total_rotation = quatmul( total_rotation, rotation)
            call rotate_coords( coords2r, rotation)
            mapdist = sqdistsum( mapping1, coords1, coords2r)
            steps = steps + 1
         end do

         call insert_record_mapping( registry, mapping1, steps, total_rotation, 0, mapdist)

         ! Stop after max_trials trials
         if (registry%n_trials >= max_trials) exit

      end do

   else

      ! Exhaustive search: the best assignments after superposition, up to
      ! the size of the registry, ranked by distance
      call assign_atoms_global( coords1, coords2, cache_arrays, registry)

      ! Optimal rotation and distance of each ranked assignment
      do i = 1, registry%n_records
         mapping1 = registry%records(i)%mapping1
         rotation = least_rotquat( mapping1, coords1, coords2)
         coords2r = coords2
         call rotate_coords( coords2r, rotation)
         registry%records(i)%rotation = rotation
         registry%records(i)%mapdist = sqdistsum( mapping1, coords1, coords2r)
      end do

   end if
end subroutine

subroutine assign_mapping_conformer( adjcs1, adjcs2, atomtypes, coords1, coords2, mapping1, error_code)
! Symmetry-corrected assignment of two conformers at their current relative
! orientation (no alignment)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(partition_t), intent(in) :: atomtypes
   real(RK), dimension(:,:), intent(in) :: coords1, coords2
   integer(IK), dimension(:), allocatable, intent(out) :: mapping1
   integer(IK), intent(out) :: error_code

   ! Local variables
   type(chaintree_node_t), pointer :: hna_chain
   type(array_trees_t) :: cache_arrays
   real(RK) :: mapdist

   ! Self-consistent HNA partition and assignment tree
   call compute_sc_hna_chain( adjcs1, adjcs2, atomtypes, hna_chain, error_code)
   if (error_code /= MOLALIGN_SUCCESS) return
   call build_assignment_tree( adjcs1, adjcs2, hna_chain%last_link, cache_arrays, error_code)
   if (error_code /= MOLALIGN_SUCCESS) return

   if (printassigntree_flag) then
      call print_chain_tree_array( atomtypes, cache_arrays)
   end if

   call assign_atoms_local( coords1, coords2, cache_arrays, mapping1, mapdist)
end subroutine

end module
