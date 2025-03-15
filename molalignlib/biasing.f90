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
use partitioning
implicit none
contains

! Iteratively compute MNA types
subroutine compute_mna_biases( mol1, mol2, eltypes, biases)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_container), intent(in) :: eltypes
   type(int_matrix), dimension(:), allocatable, intent(out) :: biases
   ! Local variables
   type(tree_node), pointer :: mnatypes
   type(tree_node_ptr), dimension(:), allocatable :: itemdir1, itemdir2
   integer :: h, i, j, iatom, jatom, level

   allocate (biases(eltypes%num_parts))

   do h = 1, eltypes%num_parts
      allocate (biases(h)%ee(eltypes%parts(h)%num_items1, eltypes%parts(h)%num_items2))
      biases(h)%ee = 0
   end do

   call tree_from_partition(eltypes, mnatypes)
   level = 0

   do

!      write (stderr, *)
!      write (stderr, '(a)') repeat('-- level '//str(level)//' --', 6)
!      call print_tree(mnatypes)

      itemdir1 = mnatypes%itemdir1
      itemdir2 = mnatypes%itemdir2
      ! Compute next level MNA types
      call compute_nextlevelmnatypes(mol1, mol2, itemdir1, itemdir2, mnatypes)
      ! Exit loop if types did not change
      if (all(mnatypes%itemdir1 == itemdir1) .and. &
          all(mnatypes%itemdir2 == itemdir2)) exit

      do h = 1, eltypes%num_parts
         do j = 1, eltypes%parts(h)%num_items2
            jatom = eltypes%parts(h)%indices2(j)
            do i = 1, eltypes%parts(h)%num_items1
               iatom = eltypes%parts(h)%indices1(i)
               if (.not. associated(mnatypes%itemdir1(iatom)%ptr, mnatypes%itemdir2(jatom)%ptr)) then
                  biases(h)%ee(i, j) = biases(h)%ee(i, j) + 1
               end if
            end do
         end do
      end do

      level = level + 1

   end do

!   do h = 1, eltypes%num_parts
!      write (stderr, *)
!      do j = 1, eltypes%parts(h)%num_items2
!         write (stderr, '(*(i2))') biases(h)%ee(:eltypes%parts(h)%num_items1, j)
!      end do
!   end do

end subroutine

end module
