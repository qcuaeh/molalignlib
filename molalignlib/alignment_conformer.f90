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
use molecule
use chemistry
use permutation
use euclidean
use assignment_atoms
use adjacency
use biasing
use pruning
use lcrs_trees
use lcrs_arrays
use assignment_tree
use assignment_conformer
use registration
use options
implicit none

contains

subroutine optimize_atomperm_conform( atomset1, atomset2, assign_arrays, coords1, coords2, registry)
   integer, dimension(:), intent(in) :: atomset1, atomset2
   type(array_trees_t), intent(inout) :: assign_arrays
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(registry_t), target, intent(out) :: registry

   ! Local variables
   integer, dimension(:), allocatable :: atomperm, new_atomperm
   real(rk), dimension(:,:), allocatable :: coords2r
   real(rk) :: permdist, new_permdist
   real(rk), dimension(4) :: rotation, total_rotation
   integer, pointer :: num_trials, lead_count
   integer :: num_steps

   ! Initialize random number generator
   call random_initialize()

   ! Initialize local minima registry
   call init_rmsd_registry( registry, num_records)
   num_trials => registry%num_trials
   lead_count => registry%records(1)%count

   ! Optimize atom permutation
   do while (lead_count < count_thres .and. num_trials < max_trials)

      ! Get randomly rotated coords2
      total_rotation = randrotquat()
      coords2r = rotated_coords( coords2, total_rotation)

      ! Assign atoms with current orientation
!      call assign_atoms_local( coords1, coords2r, assign_arrays, atomperm, permdist)
      call assign_atoms_greedy( coords1, coords2r, assign_arrays, atomperm, permdist)
      call assign_atoms_local_pruned( coords1, coords2r, assign_arrays, atomperm, permdist)
      rotation = least_rotquat( atomset1, atomperm, coords1, coords2r)
      total_rotation = quatmul( total_rotation, rotation)
      call rotate_coords( atomset2, coords2r, rotation)
      permdist = sqdistsum( atomset1, atomperm, coords1, coords2r)
      num_steps = 1

      if (iterate_flag) then
         do
!            call assign_atoms_local( coords1, coords2r, assign_arrays, new_atomperm, new_permdist)
            new_permdist = permdist
            call assign_atoms_local_pruned( coords1, coords2r, assign_arrays, new_atomperm, new_permdist)
!            write (stdout,*) permdist, new_permdist
            if (all(atomperm == new_atomperm)) exit
            atomperm = new_atomperm
            rotation = least_rotquat( atomset1, atomperm, coords1, coords2r)
            total_rotation = quatmul( total_rotation, rotation)
            call rotate_coords( atomset2, coords2r, rotation)
            permdist = sqdistsum( atomset1, atomperm, coords1, coords2r)
            num_steps = num_steps + 1
         end do
      end if

      ! Update results
      call insert_record( registry, atomperm, num_steps, permdist=permdist, rotation=total_rotation)

   end do
end subroutine

end module
