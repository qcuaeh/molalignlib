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

module biasing
use parameters
use options
use derived_types
use sorting
use utils
use permutation
use molecule
use lcrs_trees
use partitioning
use hna
use euclidean
implicit none
contains

subroutine compute_hna_biases(adjcs1, adjcs2, atomtypes, biases)
! Iteratively compute HNA types
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(partition_t), intent(in) :: atomtypes
   type(int_matrix), dimension(:), allocatable, intent(out) :: biases
   ! Local variables
   type(assigntree_node_t), pointer :: hnachain
   integer :: h, i, j, iatom, jatom
   integer :: num_splits
!   integer :: link_idx

   allocate(biases(atomtypes%num_parts))

   do h = 1, atomtypes%num_parts
      allocate(biases(h)%ee(atomtypes%parts(h)%num_items1, atomtypes%parts(h)%num_items2))
      biases(h)%ee = 0
   end do

   ! Initialize HNA chain with element types
   hnachain => chain_from_partition(atomtypes)

!   link_idx = 0
   do

!      write(stderr, *)
!      write(stderr, '(a)') repeat('-- link_idx '//str(link_idx)//' --', 6)

      ! Call refine_hna_partition and get the number of splits
      call refine_hna_partition(adjcs1, adjcs2, hnachain, num_splits)

      ! Exit loop if no splits occurred in the last iteration
      if (num_splits == 0) exit

      ! Update biases with HNAs at current level
      do h = 1, atomtypes%num_parts
         do j = 1, atomtypes%parts(h)%num_items2
            jatom = atomtypes%parts(h)%items2(j)
            do i = 1, atomtypes%parts(h)%num_items1
               iatom = atomtypes%parts(h)%items1(i)
               if (associated( &
                  hnachain%last_link%itemdir1(iatom)%ptr, &
                  hnachain%last_link%itemdir2(jatom)%ptr) &
               ) then
                  biases(h)%ee(i, j) = biases(h)%ee(i, j) + 1
               end if
            end do
         end do
      end do

!      link_idx = link_idx + 1
   end do

!   do h = 1, atomtypes%num_parts
!      write(stderr, *)
!      do j = 1, atomtypes%parts(h)%num_items2
!         write(stderr, '(*(i2))') biases(h)%ee(:atomtypes%parts(h)%num_items1, j)
!      end do
!   end do

   call delete_chain(hnachain)  ! Cleanup
end subroutine

end module
