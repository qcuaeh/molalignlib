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
use types_linked
use refinement
use pruning_atoms
use recording
use options
implicit none

contains

subroutine optimize_atomperm_atoms(atomset1, atomset2, atomtypes, prunes, coords1, coords2, registry)
   integer(ik), dimension(:), intent(in) :: atomset1, atomset2
   type(partition_t), intent(in) :: atomtypes
   type(bool_matrix), dimension(:), intent(in) :: prunes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(registry_t), intent(inout) :: registry

   ! Local variables
   integer(ik), dimension(:), allocatable :: atomperm1, new_atomperm
   real(rk), dimension(:,:), allocatable :: coords2r
   real(rk) :: permdist, steps, rotation(4), total_rotation(4)

   ! Allocations
   allocate (coords2r, mold=coords2)

   ! Initialize local minima registry
   call reset_registry( registry)

   ! Initialize random number generator
   call random_initialize()

   ! Optimize atom permutation
   do while (registry%records(1)%freq < ato_thres .and. registry%num_trials < max_trials)

      ! Apply random rotation to coords2 copy
      coords2r = coords2
      total_rotation = randrotquat()
      call rotate_coords( atomset2, coords2r, total_rotation)

      ! Assign atoms with current orientation
      call assign_atoms_pruned( atomtypes, coords1, coords2r, prunes, atomperm1)
      rotation = least_rotquat( atomset1, atomperm1, coords1, coords2r)
      call rotate_coords( atomset1, coords2r, rotation)
      total_rotation = quatmul( total_rotation, rotation)
      steps = 1

      do
         call assign_atoms_pruned( atomtypes, coords1, coords2r, prunes, new_atomperm)
         if (all(new_atomperm == atomperm1)) exit
         atomperm1 = new_atomperm
         rotation = least_rotquat( atomset1, atomperm1, coords1, coords2r)
         call rotate_coords( atomset1, coords2r, rotation)
         total_rotation = quatmul( total_rotation, rotation)
         steps = steps + 1
      end do

      ! Push local minimum to registry
      permdist = sqrt( sqdistsum( atomset1, atomperm1, coords1, coords2r))
      call insert_record_atomperm( registry, atomperm1, steps, total_rotation, 0, permdist)

   end do

end subroutine

end module
