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

module biasing_isomer
use parameters
use types_basic
use sorting
use utils
use permutation
use adjacency
use types_linked
use refining
use euclidean
use options
use randlib
implicit none

contains

subroutine init_costs(atomtypes, costs)
   type(partition_t), intent(in) :: atomtypes
   type(real_matrix), dimension(:), allocatable, intent(out) :: costs
   ! Local variables
   integer :: h

   allocate(costs(atomtypes%num_parts))

   do h = 1, atomtypes%num_parts
      allocate(costs(h)%a(atomtypes%parts(h)%num_items1, atomtypes%parts(h)%num_items2))
      costs(h)%a = 0
   end do
end subroutine

real(rk) function longest_distance(atomtypes, coords1, coords2)
   type(partition_t), target, intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   ! Local variables
   type(partition_part_t), pointer :: part
   real(rk) :: length, maxlength1, maxlength2
   integer :: h, i

   maxlength1 = 0
   maxlength2 = 0
   do h = 1, atomtypes%num_parts
      part => atomtypes%parts(h)
      do i = 1, part%num_items1
         length = sqrt(sum(coords1(:,part%items1(i))**2))
         if (length > maxlength1) then
            maxlength1 = length
         end if
      end do
      do i = 1, part%num_items2
         length = sqrt(sum(coords2(:,part%items2(i))**2))
         if (length > maxlength2) then
            maxlength2 = length
         end if
      end do
   end do

   longest_distance = maxlength1 + maxlength2
end function

subroutine add_euclidean_costs(atomtypes, coords1, coords2, euclidean_scale, costs)
   type(partition_t), target, intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   real(rk), intent(in) :: euclidean_scale
   type(real_matrix), dimension(:), allocatable, intent(inout) :: costs
   ! Local variables
   type(partition_part_t), pointer :: part
   integer :: h, i, j

   do h = 1, atomtypes%num_parts
      part => atomtypes%parts(h)
      do j = 1, part%num_items2
         do i = 1, part%num_items1
            costs(h)%a(i,j) = costs(h)%a(i,j) + euclidean_scale * &
                  sum((coords1(:,part%items1(i)) - coords2(:,part%items2(j)))**2)
         end do
      end do
   end do
end subroutine

subroutine add_random_costs(atomtypes, coords1, coords2, costs)
   type(partition_t), target, intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(real_matrix), dimension(:), allocatable, intent(inout) :: costs
   ! Local variables
   type(partition_part_t), pointer :: part
   integer :: h, i, j

   do h = 1, atomtypes%num_parts
      part => atomtypes%parts(h)
      do j = 1, part%num_items2
         do i = 1, part%num_items1
            costs(h)%a(i,j) = costs(h)%a(i,j) + random_uniform_integer(0, 1)
         end do
      end do
   end do
end subroutine

subroutine add_hna_costs(atomtypes, adjcs1, adjcs2, costs)
! Iteratively compute HNAs
   type(partition_t), target, intent(in) :: atomtypes
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(real_matrix), dimension(:), allocatable, intent(inout) :: costs
   ! Local variables
   type(partition_part_t), pointer :: part
   type(chaintree_node_t), pointer :: hna_chain
   integer :: h, i, j, iatom, jatom, num_splits
!   integer :: link_idx

   ! Initialize HNA chain with element types
   hna_chain => chain_from_partition(atomtypes)

!   link_idx = 0
   do

!      write(stderr, *)
!      write(stderr, '(a)') repeat('-- link_idx '//str(link_idx)//' --', 6)

      ! Call refine_hna_partition and get the number of splits
      call refine_hna_partition(adjcs1, adjcs2, hna_chain, num_splits)

      ! Exit loop if no splits occurred in the last iteration
      if (num_splits == 0) exit

      ! Update costs with HNAs at current level
      do h = 1, atomtypes%num_parts
         part => atomtypes%parts(h)
         do j = 1, part%num_items2
            jatom = part%items2(j)
            do i = 1, part%num_items1
               iatom = part%items1(i)
               if (.not. associated( &
                  hna_chain%last_link%itemdir1(iatom)%ptr, &
                  hna_chain%last_link%itemdir2(jatom)%ptr) &
               ) then
                  costs(h)%a(i, j) = costs(h)%a(i, j) + 1
               end if
            end do
         end do
      end do

!      link_idx = link_idx + 1
   end do

!   do h = 1, atomtypes%num_parts
!      write(stderr, *)
!      do j = 1, atomtypes%parts(h)%num_items2
!         write(stderr, '(*(i2))') costs(h)%a(:atomtypes%parts(h)%num_items1, j)
!      end do
!   end do

   call delete_chain(hna_chain)  ! Cleanup
end subroutine

end module
