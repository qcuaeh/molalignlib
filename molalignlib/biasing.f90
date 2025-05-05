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

subroutine compute_mna_biases(mol1, mol2, eltypes, biases)
! Iteratively compute MNA types
   type(mol_type), intent(in) :: mol1, mol2
   type(item_partition), intent(in) :: eltypes
   type(int_matrix), dimension(:), allocatable, intent(out) :: biases
   ! Local variables
   type(poly_node), pointer :: mnapolytree
   integer :: level, prev_num_leaves
   integer :: h, i, j, iatom, jatom

   allocate(biases(eltypes%num_parts))

   do h = 1, eltypes%num_parts
      allocate(biases(h)%ee(eltypes%parts(h)%num_items1, eltypes%parts(h)%num_items2))
      biases(h)%ee = 0
   end do

   ! Initialize mnapolytree with first tree from partition
   mnapolytree => make_new_poly( size(eltypes%itemdir1), size(eltypes%itemdir2))
   call add_root( mnapolytree, tree_from_partition( eltypes))
   level = 0

   do
!      write(stderr, *)
!      write(stderr, '(a)') repeat('-- level '//str(level)//' --', 6)
!      call print_tree(mnapolytree%last_root)

      ! Compute MNA upper level types
      prev_num_leaves = mnapolytree%last_root%num_leaves
      call compute_nextlevel_mnas(mol1, mol2, mnapolytree)

      ! Exit loop if types did not change
      if (mnapolytree%last_root%num_leaves == prev_num_leaves) exit

      level = level + 1

      do h = 1, eltypes%num_parts
         do j = 1, eltypes%parts(h)%num_items2
            jatom = eltypes%parts(h)%items2(j)
            do i = 1, eltypes%parts(h)%num_items1
               iatom = eltypes%parts(h)%items1(i)
               if (.not. associated( &
                  mnapolytree%last_root%itemdir1(iatom)%ptr, &
                  mnapolytree%last_root%itemdir2(jatom)%ptr) &
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

   call delete_polytree(mnapolytree)  ! Added cleanup
end subroutine

end module
