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

module alignment_cluster
use parameters
use options
use random
use molecule
use chemistry
use permutation
use euclidean
use assignment_cluster
use lcrs_trees
use hna
use partitioning
use pruning
use registration
implicit none

contains

subroutine optimize_atomperm_cluster(coords1, coords2, atomtypes, prunes, registry)
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(partition_t), intent(in) :: atomtypes
   type(registry_t), intent(out) :: registry
   type(bool_matrix), dimension(:), intent(in) :: prunes

   ! Local variables
   integer :: num_steps
   type(subperm_t) :: atomperm, new_atomperm
   real(rk), dimension(:,:), allocatable :: coords2r
   real(rk) :: permdist, rotation_step(4), rotation(4)
   integer :: insert_pos

   ! Initialize random number generator
   call random_initialize()

   ! Initialize local minima registry
   call init_rmsd_registry( registry, max_records)

   ! Optimize atom permutation
   do while (registry%records(1)%count < count_thres .and. registry%num_trials < max_trials)

      ! Aply a random rotation to coords2
      rotation = randrotquat()
      coords2r = rotated_coords( coords2, rotation)

      ! Assign atoms with current orientation
      call assign_atoms_pruned( atomtypes, coords1, coords2r, prunes, atomperm)
      rotation_step = least_rotquat( atomperm, coords1, coords2r)
      call rotate_coords( coords2r, rotation_step)
      rotation = quatmul( rotation, rotation_step)
      num_steps = 1

      if (iterate_flag) then
         do
            call assign_atoms_pruned( atomtypes, coords1, coords2r, prunes, new_atomperm)
            if (new_atomperm == atomperm) exit
            atomperm = new_atomperm
            rotation_step = least_rotquat( atomperm, coords1, coords2r)
            call rotate_coords( coords2r, rotation_step)
            rotation = quatmul( rotation, rotation_step)
            num_steps = num_steps + 1
         end do
      end if

      ! Push local minimum to registry
      permdist = sqrt( sqdistsum( atomperm, coords1, coords2r))
      insert_pos = insert_record( registry, atomperm, num_steps, permdist=permdist, rotation=rotation)

   end do

end subroutine

end module
