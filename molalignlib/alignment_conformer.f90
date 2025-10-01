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
use options
use random
use molecule
use utils
use chemistry
use permutation
use euclidean
use assignment_cluster
use adjacency
use biasing
use pruning
use lcrs_trees
use lcrs_arrays
use assignment_tree
use assignment_conformer
use registration
implicit none

contains

subroutine optimize_atomperm_conform( coords1, coords2, assign_arrays, registry)
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(array_trees_t), intent(inout) :: assign_arrays
   type(registry_t), target, intent(out) :: registry

   ! Local variables
   type(subperm_t) :: atomperm, new_atomperm
   real(rk), dimension(:,:), allocatable :: coords2r
   real(rk) :: permdist, prune_thres, rotation_step(4), rotation(4)
   integer, pointer :: num_trials, lead_count
   integer :: num_steps

   ! Initialize random number generator
   call random_initialize()

   ! Initialize local minima registry
   call init_rmsd_registry( registry, max_records)
   num_trials => registry%num_trials
   lead_count => registry%records(1)%count

   ! Optimize atom permutation
   do while (lead_count < count_thres .and. num_trials < max_trials)

      ! Get randomly rotated coords2
      rotation = randrotquat()
      coords2r = rotated_coords( coords2, rotation)

      ! Assign atoms with current orientation
!      call assign_atoms_local( coords1, coords2r, assign_arrays, atomperm)
      call assign_atoms_greedy( coords1, coords2r, assign_arrays, atomperm, prune_thres)
      call assign_atoms_local_pruned( coords1, coords2r, assign_arrays, prune_thres, atomperm)
      rotation_step = least_rotquat( atomperm, coords1, coords2r)
      call rotate_coords( coords2r, rotation_step)
      rotation = quatmul( rotation, rotation_step)
      num_steps = 1

      if (iterate_flag) then
         do
!            call assign_atoms_local( coords1, coords2r, assign_arrays, new_atomperm)
            prune_thres = sqdistsum( atomperm, coords1, coords2r)
            call assign_atoms_local_pruned( coords1, coords2r, assign_arrays, prune_thres, new_atomperm)
            if (new_atomperm == atomperm) exit
            atomperm = new_atomperm
            rotation_step = least_rotquat( atomperm, coords1, coords2r)
            call rotate_coords( coords2r, rotation_step)
            rotation = quatmul( rotation, rotation_step)
            num_steps = num_steps + 1
         end do
      end if

      ! Update results
      permdist = sqdistsum( atomperm, coords1, coords2r)
      call push_record( registry, atomperm, num_steps, permdist=permdist, rotation=rotation)

   end do
end subroutine

end module
