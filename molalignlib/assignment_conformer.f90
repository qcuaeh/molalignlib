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

module assignment_conformer
use parameters
use derived_types
use options
use random
use molecule
use utils
use chemistry
use permutation
use spatial_transforms
use assignment_atoms
use adjacency
use biasing
use pruning
use lcrs_tree
use lcrs_tree_arrays
use assigntree_build
use assigntree_distribute_linked
use assigntree_distribute
use registration

implicit none
logical, parameter :: iter_flag = .true.

contains

subroutine optimize_atomperm_conform( coords1, coords2, assign_arrays, registry)
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(array_trees_t), intent(inout) :: assign_arrays
   type(registry_t), target, intent(out) :: registry

   ! Local variables
   type(subperm_t) :: atomperm, auxperm
   real(rk), dimension(:,:), allocatable :: coords2r
   real(rk) :: rmsd, rotation_step(4), rotation(4)
   integer, pointer :: num_trials, lead_count
   integer :: num_steps

   ! Initialize random number generator
   call random_initialize()

   ! Initialize local minima registry
   call init_rmsd_registry( registry, max_records)
   num_trials => registry%num_trials
   lead_count => registry%records(1)%count

   ! Optimize atom permutation
   do while (lead_count < max_count .and. num_trials < max_trials)

      ! Get randomly rotated coords2
      rotation = randrotquat()
      coords2r = rotated_coords( coords2, rotation)

      ! Assign atoms with current orientation
      call distribute_items_parallel( coords1, coords2r, assign_arrays, atomperm)
      rotation_step = least_rotquat( atomperm, coords1, coords2r)
      call rotate_coords( coords2r, rotation_step)
      rotation = quatmul( rotation, rotation_step)
      num_steps = 1

      do while (iter_flag)
         call distribute_items_parallel( coords1, coords2r, assign_arrays, auxperm)
!         write (stderr,'(F8.4)') sqrt( total_sqdist( auxperm, coords1, coords2r))
         if (auxperm == atomperm) exit
         atomperm = auxperm
         rotation_step = least_rotquat( atomperm, coords1, coords2r)
         call rotate_coords( coords2r, rotation_step)
         rotation = quatmul( rotation, rotation_step)
         num_steps = num_steps + 1
      end do

      ! Update results
      rmsd = sqrt( total_sqdist( atomperm, coords1, coords2r))
      call push_record( registry, atomperm, num_steps, rmsd=rmsd, rotation=rotation)

   end do
end subroutine

end module
