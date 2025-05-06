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
use globals
use basetypes
use sorting
use strutils
use molecule
use lcrs_tree
use mna_compute
implicit none
contains

subroutine compute_mna_biases(mol1, mol2, eltypes, biases)
! Iteratively compute MNA types
   type(mol_type), intent(in) :: mol1, mol2
   type(partitionarray_t), intent(in) :: eltypes
   type(int_matrix), dimension(:), allocatable, intent(out) :: biases
   ! Local variables
   type(tree_node_t), pointer :: mnachain
   integer :: link_idx, prev_num_parts
   integer :: h, i, j, iatom, jatom

   allocate(biases(eltypes%num_parts))

   do h = 1, eltypes%num_parts
      allocate(biases(h)%ee(eltypes%parts(h)%num_items1, eltypes%parts(h)%num_items2))
      biases(h)%ee = 0
   end do

   ! Initialize mna chain with element types
   mnachain => tree_from_partitionarray(eltypes)

   link_idx = 0
   do

!      write(stderr, *)
!      write(stderr, '(a)') repeat('-- link_idx '//str(link_idx)//' --', 6)
!      call print_link(mnachain%last_link)

      ! Save the current number of parts before computation
      prev_num_parts = mnachain%last_link%num_parts

      ! Compute next level and update current branch
      call compute_nextlevel_mnas(mol1, mol2, mnachain)

      ! Exit the loop if no change
      if (mnachain%last_link%num_parts == prev_num_parts) exit

      link_idx = link_idx + 1

      ! Update biases based on the current branch
      do h = 1, eltypes%num_parts
         do j = 1, eltypes%parts(h)%num_items2
            jatom = eltypes%parts(h)%items2(j)
            do i = 1, eltypes%parts(h)%num_items1
               iatom = eltypes%parts(h)%items1(i)
               if (.not. associated( &
                  mnachain%last_link%itemdir1(iatom)%ptr, &
                  mnachain%last_link%itemdir2(jatom)%ptr) &
               ) then
                  biases(h)%ee(i, j) = biases(h)%ee(i, j) + 1
               end if
            end do
         end do
      end do

   end do

!   do h = 1, eltypes%num_parts
!      write(stderr, *)
!      do j = 1, eltypes%parts(h)%num_items2
!         write(stderr, '(*(i2))') biases(h)%ee(:eltypes%parts(h)%num_items1, j)
!      end do
!   end do

   call delete_tree(mnachain)  ! Cleanup
end subroutine

end module
