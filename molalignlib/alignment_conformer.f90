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

module alignment_conformer
use parameters
use derived_types
use utils
use random
use chemistry
use permutation
use euclidean
use assignment_atoms
use adjacency
use biasing_isomer
use pruning_atoms
use lcrs_trees
use lcrs_arrays
use assignment_tree
use assignment_conformer
use recording
use options
implicit none

contains

subroutine optimize_atomperm_conform( atomset1, atomset2, adjcs1, adjcs2, atomtypes, coords1, coords2, registry)
   integer, dimension(:), intent(in) :: atomset1, atomset2
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(partition_t), intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(registry_t), target, intent(out) :: registry

   ! Local variables
   integer, dimension(:), allocatable :: atomperm1, new_atomperm
   real(rk), dimension(:,:), allocatable :: coords2r
   real(rk) :: permdist, new_permdist
   real(rk), dimension(4) :: rotation, total_rotation
   integer, pointer :: num_trials, lead_count
   integer :: steps
   type(assigntree_node_t), pointer :: hnachain
   type(array_trees_t) :: assign_arrays
   logical :: stoch_flag, count_flag

   stoch_flag = .true.
   count_flag = .true.

   ! Pre-compute assignment tree for decision making
   call compute_consistent_hna_partition( adjcs1, adjcs2, atomtypes, hnachain)
   call build_assignment_tree( adjcs1, adjcs2, hnachain%last_link, assign_arrays)

   if (tree_flag) then
      call print_chain_tree_array( assign_arrays)
   end if

   if ((stoch_flag .and. .not. count_flag) .or. (stoch_flag .and. count_flag .and. &
         assign_arrays%global_combinations > count_thres*assign_arrays%local_combinations)) then

      ! Initialize random number generator
      call random_initialize()

      ! Initialize local minima registry
      call init_registry( registry, num_records)
      num_trials => registry%num_trials
      lead_count => registry%records(1)%count

      ! Optimize atom permutation
      do while (lead_count < count_thres .and. num_trials < max_trials)

         ! Get randomly rotated coords2
         total_rotation = randrotquat()
         coords2r = rotated_coords( coords2, total_rotation)

         ! Assign atoms with current orientation
         if (PRUNE_ASSIGNMENT_TREE) then
            call assign_atoms_greedy( coords1, coords2r, assign_arrays, atomperm1, permdist)
            call assign_atoms_local_pruned( coords1, coords2r, assign_arrays, atomperm1, permdist)
         else
            call assign_atoms_local( coords1, coords2r, assign_arrays, atomperm1, permdist)
         end if
         rotation = least_rotquat( atomset1, atomperm1, coords1, coords2r)
         total_rotation = quatmul( total_rotation, rotation)
         call rotate_coords( atomset2, coords2r, rotation)
         permdist = sqdistsum( atomset1, atomperm1, coords1, coords2r)
         steps = 1

         if (iterate_flag) then
            do
               if (PRUNE_ASSIGNMENT_TREE) then
                  new_permdist = permdist
                  call assign_atoms_local_pruned( coords1, coords2r, assign_arrays, new_atomperm, new_permdist)
               else
                  call assign_atoms_local( coords1, coords2r, assign_arrays, new_atomperm, new_permdist)
               end if
!               write (stdout,*) permdist, new_permdist
               if (all(atomperm1 == new_atomperm)) exit
               atomperm1 = new_atomperm
               rotation = least_rotquat( atomset1, atomperm1, coords1, coords2r)
               total_rotation = quatmul( total_rotation, rotation)
               call rotate_coords( atomset2, coords2r, rotation)
               permdist = sqdistsum( atomset1, atomperm1, coords1, coords2r)
               steps = steps + 1
            end do
         end if

         ! Update results
         call insert_record_homo( registry, atomperm1, permdist, steps, total_rotation)

      end do

   else

      ! Assign atoms using global assignment
      call assign_atoms_global( coords1, coords2, assign_arrays, atomperm1)
      
      ! Calculate optimal rotation
      rotation = least_rotquat( atomset1, atomperm1, coords1, coords2)
      
      ! Rotate coords2 and calculate permdist
      coords2r = rotated_coords( coords2, rotation)
      permdist = sqdistsum( atomset1, atomperm1, coords1, coords2r)
      
      ! Initialize registry with a single record
      call init_registry( registry, 1)
      registry%occ_records = 1
      registry%records(1)%atomperm1 = atomperm1
      registry%records(1)%permdist = permdist
      registry%records(1)%count = 1
      registry%records(1)%steps = 1
      registry%records(1)%rotation = rotation

   end if
end subroutine

end module
