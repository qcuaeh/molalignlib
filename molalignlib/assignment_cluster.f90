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

module assignment_cluster
use parameters
use options
use random
use molecule
use chemistry
use permutation
use spatial_transforms
use assignment_atoms
use lcrs_tree
use atom_mnas
use atom_types
use pruning
use registration

implicit none
logical, parameter :: iter_flag = .true.

contains

subroutine optimize_atomperm_cluster(coords1, coords2, atomtypes, prunes, registry)
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(partition_t), intent(in) :: atomtypes
   type(registry_t), intent(out) :: registry
   type(bool_matrix), dimension(:), intent(in) :: prunes

   ! Local variables
   integer :: num_steps
   type(subperm_t) :: atomperm, auxperm
   real(rk), dimension(:,:), allocatable :: coords2r
   real(rk) :: rmsd, rotation_step(4), rotation(4)

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

      do while (iter_flag)
         call assign_atoms_pruned( atomtypes, coords1, coords2r, prunes, auxperm)
         if (auxperm == atomperm) exit
         atomperm = auxperm
         rotation_step = least_rotquat( atomperm, coords1, coords2r)
         call rotate_coords( coords2r, rotation_step)
         rotation = quatmul( rotation, rotation_step)
         num_steps = num_steps + 1
      end do

      ! Push local minimum to registry
      rmsd = sqrt( total_sqdist( atomperm, coords1, coords2r))
      call push_record( registry, atomperm, num_steps, rmsd=rmsd, rotation=rotation)

   end do

end subroutine

end module
