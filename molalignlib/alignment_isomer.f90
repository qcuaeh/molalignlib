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

module alignment_isomer
use parameters
use types_basic
use random
use euclidean
use adjacency
use biasing_isomer
use partitioning
use recording
use assignment_atoms
use assignment_bonds
use options
implicit none

contains

subroutine optimize_atomperm_isomer( atomset1, atomset2, atomtypes, adjcs1, adjcs2, &
                                     coords1, coords2, registry)
   integer, dimension(:), intent(in) :: atomset1, atomset2
   type(partition_t), intent(in) :: atomtypes
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(registry_t), target, intent(inout) :: registry

   ! Local variables
   integer, dimension(:), allocatable :: atomperm1, new_atomperm
   real(rk), dimension(:,:), allocatable :: coords2r
   real(rk), dimension(4) :: rotation, total_rotation
   real(rk) :: euclidean_scale, permdist
   integer :: permdiff
   type(real_matrix), dimension(:), allocatable :: costs, fixed_costs
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   integer, dimension(:,:), allocatable :: moldiffs
   type(partition_t) :: scnatypes

   ! Convert adjacency lists to adjacency matrix
   adjmat1 = adjcs_to_adjmat( adjcs1)
   adjmat2 = adjcs_to_adjmat( adjcs2)

   ! Initialize atomperm1 as an identity permutation
   allocate (atomperm1(size(coords1, 2)))
   allocate (new_atomperm(size(coords1, 2)))
   call init_identity_permutation( atomperm1)

   ! Compute constant costs
   call init_costs( atomtypes, fixed_costs)
   call add_hna_costs( atomtypes, adjcs1, adjcs2, scnatypes, fixed_costs)
   euclidean_scale = 0.99_rk/longest_distance( atomtypes, coords1, coords2)**2

   ! Initialize local minima registry
   call reset_registry( registry)

   ! Initialize random number generator
   call random_initialize()

   ! Optimize atom permutation
   do while (registry%records(1)%freq < iso_thres .and. registry%num_trials < max_trials)

      ! Get randomly rotated coords2
      total_rotation = randrotquat()
      coords2r = rotated_coords( coords2, total_rotation)

      ! Assign atoms with current orientation using costs
      costs = fixed_costs
      call add_euclidean_costs( atomtypes, coords1, coords2r, euclidean_scale, costs)
!      call add_random_costs( atomtypes, coords1, coords2r, costs)
      call assign_atoms( atomtypes, costs, atomperm1)
      call minimize_adjdiff( atomset1, atomtypes, scnatypes, adjcs1, adjcs2, adjmat2, &
            coords1, coords2, atomperm1)
      rotation = least_rotquat( atomset1, atomperm1, coords1, coords2r)
      total_rotation = quatmul( total_rotation, rotation)
      call rotate_coords( atomset2, coords2r, rotation)
      permdiff = adjacencydiff( atomset1, atomperm1, adjcs1, adjcs2)

      ! Insert atom permutation into registry
      permdist = sqdistsum( atomset1, atomperm1, coords1, coords2r)
      call compute_differing_bonds( atomset1, atomperm1, adjmat1, adjmat2, moldiffs)
      call insert_record_moldiff( registry, moldiffs, atomperm1, 1, total_rotation, permdist)

   end do
end subroutine

end module
