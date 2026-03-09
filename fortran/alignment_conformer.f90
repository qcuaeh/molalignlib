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
use types_basic
use str_utils
use random
use chemdata
use permutation
use euclidean
use adjacency
use pruning_atoms
use types_linked
use types_indexed
use refinement
use assignment_atoms
use assignment_conformer
use recording
use options
implicit none

contains

subroutine optimize_atomperm_conformer( atomset1, atomset2, adjcs1, adjcs2, atomtypes, &
      coords1, coords2, registry)
   integer(ik), dimension(:), intent(in) :: atomset1, atomset2
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(partition_t), intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(registry_t), target, intent(inout) :: registry

   ! Local variables
   integer(ik), dimension(:), allocatable :: atomperm1, new_atomperm
   real(rk), dimension(:,:), allocatable :: coords2r
   real(rk), dimension(4) :: rotation, total_rotation
   real(rk) :: steps, permdist, new_permdist
   type(chaintree_node_t), pointer :: hna_chain
   type(array_trees_t) :: cache_arrays

   ! Allocations
   allocate (coords2r, mold=coords2)

   ! Pre-compute assignment tree for decision making
   call compute_sc_hna_chain( adjcs1, adjcs2, atomtypes, hna_chain)
   call build_assignment_tree( adjcs1, adjcs2, hna_chain%last_link, cache_arrays)

   if (print_tree_flag) then
      call print_chain_tree_array( atomtypes, cache_arrays)
   end if

   ! Reset registry for new conformer
   call reset_registry( registry)

   if ((stoch_flag .and. .not. adaptive_flag) .or. (stoch_flag .and. adaptive_flag .and. &
         cache_arrays%total_combinations > confo_thres*cache_arrays%partial_combinations)) then

      ! Initialize random number generator
      call random_initialize()

      ! Optimize atom permutation
      do

         ! Apply random rotation to coords2 copy
         coords2r = coords2
         total_rotation = randrotquat()
         call rotate_coords( atomset2, coords2r, total_rotation)

         ! Assign atoms with current orientation
         if (PRUNE_ASSIGNMENTS) then
            call assign_atoms_greedy( coords1, coords2r, cache_arrays, atomperm1, permdist)
            call assign_atoms_local_pruned( coords1, coords2r, cache_arrays, atomperm1, permdist)
         else
            call assign_atoms_local( coords1, coords2r, cache_arrays, atomperm1, permdist)
         end if

         rotation = least_rotquat( atomset1, atomperm1, coords1, coords2r)
         total_rotation = quatmul( total_rotation, rotation)
         call rotate_coords( atomset2, coords2r, rotation)
         permdist = sqdistsum( atomset1, atomperm1, coords1, coords2r)
         steps = 1

         do
            if (PRUNE_ASSIGNMENTS) then
               new_permdist = permdist
               call assign_atoms_local_pruned( coords1, coords2r, cache_arrays, new_atomperm, new_permdist)
            else
               call assign_atoms_local( coords1, coords2r, cache_arrays, new_atomperm, new_permdist)
            end if
!            write (stdout,*) permdist, new_permdist
            if (all(atomperm1 == new_atomperm)) exit
            atomperm1 = new_atomperm
            rotation = least_rotquat( atomset1, atomperm1, coords1, coords2r)
            total_rotation = quatmul( total_rotation, rotation)
            call rotate_coords( atomset2, coords2r, rotation)
            permdist = sqdistsum( atomset1, atomperm1, coords1, coords2r)
            steps = steps + 1
         end do

         ! Update results
         call insert_record_atomperm( registry, atomperm1, steps, total_rotation, 0, permdist)

         if (registry%records(1)%freq > confo_thres) then
            exit
         end if

         if (MAX_TRIALS_EXIT) then
            if (registry%num_trials > max_trials) then
               exit
            end if
         end if
      end do

   else

      ! Assign atoms using global assignment
      call assign_atoms_global( coords1, coords2, cache_arrays, atomperm1)
      
      ! Calculate optimal rotation
      rotation = least_rotquat( atomset1, atomperm1, coords1, coords2)
      
      ! Rotate coords2 and calculate permdist
      coords2r = coords2
      call rotate_coords( atomset2, coords2r, rotation)
      permdist = sqdistsum( atomset1, atomperm1, coords1, coords2r)
      
      ! Initialize registry with a single record
      registry%occ_records = 1
      registry%records(1)%atomperm1 = atomperm1
      registry%records(1)%permdiff = 0
      registry%records(1)%permdist = permdist
      registry%records(1)%freq = 1
      registry%records(1)%steps = 1
      registry%records(1)%rotation = rotation

   end if
end subroutine

subroutine assign_atomperm_conformer( adjcs1, adjcs2, atomtypes, coords1, coords2, atomperm1)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(partition_t), intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer(ik), dimension(:), allocatable, intent(out) :: atomperm1

   ! Local variables
   type(chaintree_node_t), pointer :: hna_chain
   type(array_trees_t) :: cache_arrays
   real(rk) :: permdist

   ! Pre-compute assignment tree
   call compute_sc_hna_chain( adjcs1, adjcs2, atomtypes, hna_chain)
   call build_assignment_tree( adjcs1, adjcs2, hna_chain%last_link, cache_arrays)

   if (print_tree_flag) then
      call print_chain_tree_array( atomtypes, cache_arrays)
   end if

   call assign_atoms_local( coords1, coords2, cache_arrays, atomperm1, permdist)
!   call assign_atoms_greedy( coords1, coords2, cache_arrays, atomperm1, permdist)
!   call assign_atoms_local_pruned( coords1, coords2, cache_arrays, atomperm1, permdist)
end subroutine

end module
