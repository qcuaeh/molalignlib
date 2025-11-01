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
use pruning_atoms
use lcrs_trees
use lcrs_arrays
use partitioning
use assignment_tree
use assignment_conformer
use recording
use options
implicit none

contains

subroutine optimize_atomperm_conformer( atomset1, atomset2, adjcs1, adjcs2, atomtypes, &
      coords1, coords2, registry)
   integer, dimension(:), intent(in) :: atomset1, atomset2
   type(adjcs_t), intent(in) :: adjcs1, adjcs2
   type(partition_t), intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(registry_t), target, intent(inout) :: registry

   ! Local variables
   integer, dimension(:), allocatable :: atomperm1, new_atomperm
   real(rk), dimension(:,:), allocatable :: coords2r
   real(rk), dimension(4) :: rotation, total_rotation
   real(rk) :: steps, permdist, new_permdist
   type(assigntree_node_t), pointer :: mlnachain
   type(array_trees_t) :: assign_arrays

   ! Pre-compute assignment tree for decision making
   call compute_scna_partition( adjcs1, adjcs2, atomtypes, mlnachain)
   call build_assignment_tree( adjcs1, adjcs2, mlnachain%last_link, assign_arrays)

   if (tree_flag) then
      call print_chain_tree_array( assign_arrays)
   end if

   ! Reset registry for new conformer
   call reset_registry( registry)

   if ((stoch_flag .and. .not. adaptive_flag) .or. (stoch_flag .and. adaptive_flag .and. &
         assign_arrays%global_combinations > confo_thres*assign_arrays%local_combinations)) then

      ! Initialize random number generator
      call random_initialize()

      ! Optimize atom permutation
      do while (registry%records(1)%count < confo_thres .and. registry%num_trials < max_trials)

         ! Get randomly rotated coords2
         total_rotation = randrotquat()
         coords2r = rotated_coords( coords2, total_rotation)

         ! Assign atoms with current orientation
         if (full_flag) then
            call assign_atoms_local_full( adjcs1, adjcs2, coords1, coords2r, assign_arrays, &
                  atomperm1, permdist)
         else
            call assign_atoms_greedy( adjcs1, adjcs2, coords1, coords2r, assign_arrays, &
                  atomperm1, permdist)
            call assign_atoms_local_pruned( adjcs1, adjcs2, coords1, coords2r, assign_arrays, &
                  atomperm1, permdist)
         end if
         rotation = least_rotquat( atomset1, atomperm1, coords1, coords2r)
         total_rotation = quatmul( total_rotation, rotation)
         call rotate_coords( atomset2, coords2r, rotation)
         permdist = sqdistsum( atomset1, atomperm1, coords1, coords2r)
         steps = 1

         do
            if (full_flag) then
               call assign_atoms_local_full( adjcs1, adjcs2, coords1, coords2r, assign_arrays, &
                     new_atomperm, new_permdist)
            else
               new_permdist = permdist
               call assign_atoms_local_pruned( adjcs1, adjcs2, coords1, coords2r, assign_arrays, &
                     new_atomperm, new_permdist)
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

      end do

   else

      ! Assign atoms using global assignment
      call assign_atoms_global( adjcs1, adjcs2, coords1, coords2, assign_arrays, atomperm1)
      
      ! Calculate optimal rotation
      rotation = least_rotquat( atomset1, atomperm1, coords1, coords2)
      
      ! Rotate coords2 and calculate permdist
      coords2r = rotated_coords( coords2, rotation)
      permdist = sqdistsum( atomset1, atomperm1, coords1, coords2r)
      
      ! Initialize registry with a single record
      registry%occ_records = 1
      registry%records(1)%atomperm1 = atomperm1
      registry%records(1)%permdiff = 0
      registry%records(1)%permdist = permdist
      registry%records(1)%count = 1
      registry%records(1)%steps = 1
      registry%records(1)%rotation = rotation

   end if
end subroutine

subroutine assign_atomperm_conformer( adjcs1, adjcs2, atomtypes, coords1, coords2, atomperm1)
   type(adjcs_t), intent(in) :: adjcs1, adjcs2
   type(partition_t), intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, dimension(:), allocatable, intent(out) :: atomperm1

   ! Local variables
   type(assigntree_node_t), pointer :: mlnachain
   type(array_trees_t) :: assign_arrays
   real(rk) :: permdist

   ! Pre-compute assignment tree
   call compute_scna_partition( adjcs1, adjcs2, atomtypes, mlnachain)
   call build_assignment_tree( adjcs1, adjcs2, mlnachain%last_link, assign_arrays)

   if (tree_flag) then
      call print_chain_tree_array( assign_arrays)
   end if

   ! Assign atoms using greedy and local pruned methods
   if (full_flag) then
      call assign_atoms_local_full( adjcs1, adjcs2, coords1, coords2, assign_arrays, &
            atomperm1, permdist)
   else
      call assign_atoms_greedy( adjcs1, adjcs2, coords1, coords2, assign_arrays, &
            atomperm1, permdist)
      call assign_atoms_local_pruned( adjcs1, adjcs2, coords1, coords2, assign_arrays, &
            atomperm1, permdist)
   end if
end subroutine

end module
