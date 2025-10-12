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
use derived_types
use random
use molecule
use euclidean
use adjacency
!use tracking
use biasing
use partitioning
use recording
use assignment_atoms
use options
implicit none

contains

subroutine optimize_atomperm_isomer(atomset1, atomset2, atomtypes, adjcs1, adjcs2, &
                                     coords1, coords2, registry)
   integer, dimension(:), intent(in) :: atomset1, atomset2
   type(partition_t), intent(in) :: atomtypes
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(registry_t), target, intent(out) :: registry

   ! Local variables
   integer, dimension(:), allocatable :: atomperm, new_atomperm
   real(rk), dimension(:,:), allocatable :: coords2r
   real(rk) :: permdist
   real(rk), dimension(4) :: rotation, total_rotation
   integer, pointer :: num_trials, lead_count
   integer :: steps
   integer :: permdiff
!   type(int_list), dimension(:), allocatable :: molfrags1
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   type(int_matrix), dimension(:), allocatable :: biases
   type(partition_t) :: scnatypes

   ! Initialize random number generator
   call random_initialize()

   ! Compute adjacency matrices from adjacency lists
   adjmat1 = adjcs_to_adjmat(adjcs1)
   adjmat2 = adjcs_to_adjmat(adjcs2)

   ! Find molecular fragments
!   call find_molfrags(adjcs1, first_partition(atomtypes), molfrags1)

   ! Compute HNA types and biases
   call compute_hna_biases(adjcs1, adjcs2, atomtypes, biases, scnatypes)

   ! Initialize local minima registry (dual mode: adjacency + position)
   call init_permutation_grouped_registry(registry, num_records)
   num_trials => registry%num_trials
   lead_count => registry%records(1)%count

   ! Optimize atom permutation
   do while (lead_count < count_thres .and. num_trials < max_trials)

      ! Get randomly rotated coords2
      total_rotation = randrotquat()
      coords2r = rotated_coords(coords2, total_rotation)

      ! Assign atoms with current orientation using biases
      call assign_atoms_biased(atomtypes, biases, coords1, coords2r, atomperm)
!      call minimize_adjdiff(atomtypes, scnatypes, molfrags1, adjcs1, adjcs2, &
!            coords1, coords2, atomperm)

      ! Optimize rotation
      rotation = least_rotquat(atomset1, atomperm, coords1, coords2r)
      total_rotation = quatmul(total_rotation, rotation)
      call rotate_coords(atomset2, coords2r, rotation)

      permdiff = adjacencydiff(atomperm, adjmat1, adjmat2)
      permdist = sqdistsum(atomset1, atomperm, coords1, coords2r)
      steps = 1

      if (iterate_flag) then
         do
            ! Try to improve assignment
            call assign_atoms_biased(atomtypes, biases, coords1, coords2r, new_atomperm)
!            call minimize_adjdiff(atomtypes, scnatypes, molfrags1, adjcs1, adjcs2, &
!                  coords1, coords2, atomperm)

            if (all(atomperm == new_atomperm)) exit
            atomperm = new_atomperm

            rotation = least_rotquat(atomset1, atomperm, coords1, coords2r)
            total_rotation = quatmul(total_rotation, rotation)
            call rotate_coords(atomset2, coords2r, rotation)

            permdiff = adjacencydiff(atomperm, adjmat1, adjmat2)
            permdist = sqdistsum(atomset1, atomperm, coords1, coords2r)
            steps = steps + 1
         end do
      end if

      ! Update results
      call insert_record(registry, atomperm, permdiff, permdist, steps, total_rotation)

   end do
end subroutine

function adjcs_to_adjmat(adjcs) result(adjmat)
   type(adjc_t), dimension(:), intent(in) :: adjcs
   logical, dimension(:,:), allocatable :: adjmat
   integer :: i, j, k

   allocate(adjmat(size(adjcs), size(adjcs)))
   adjmat = .false.

   do i = 1, size(adjcs)
      do j = 1, size(adjcs(i)%adjlist)
         k = adjcs(i)%adjlist(j)
         adjmat(i, k) = .true.
      end do
   end do
end function

end module
