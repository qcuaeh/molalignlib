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

module partitioning
use parameters
use sorting
use chemdata
use molecule
use lcrs_tree

implicit none

type :: typehood_item
   type(tree_node), pointer :: node
   type(tree_node_ptr), allocatable :: typehood(:)
end type

type :: typehood_table
   integer :: num_items
   type(typehood_item), allocatable :: items(:)
end type

type :: atomtype_item
   integer :: elnum
   integer :: label
   type(tree_node), pointer :: node
end type

type :: atomtype_table
   integer :: num_items
   type(atomtype_item), allocatable :: items(:)
end type

interface operator (.equiv.)
   module procedure typehood_equivalence
end interface

contains

function typehood_equivalence( array1, array2) result(equiv)
   type(tree_node_ptr), dimension(:), intent(in) :: array1, array2
   logical :: equiv
   integer :: i, j, matches

   ! Check if sizes are equal
   if (size(array1) /= size(array2)) then
      equiv = .false.
      return
   end if

   ! Early exit: Check if any value appears more times in one array
   ! than it does in the other
   do i = 1, size(array1)
       matches = 0
       do j = 1, size(array1)
           if (associated(array1(i)%ptr, array2(j)%ptr)) matches = matches + 1
           if (associated(array1(i)%ptr, array1(j)%ptr)) matches = matches - 1
       end do
       if (matches /= 0) then
           equiv = .false.
           return
       end if
   end do

   equiv = .true.

end function

subroutine add_atomtype(atomtypetable, elnum, label, node)
   type(atomtype_table), intent(inout) :: atomtypetable
   integer, intent(in) :: elnum
   integer, intent(in) :: label
   type(tree_node), pointer, intent(in) :: node

   atomtypetable%num_items = atomtypetable%num_items + 1
   atomtypetable%items(atomtypetable%num_items)%elnum = elnum
   atomtypetable%items(atomtypetable%num_items)%label = label
   atomtypetable%items(atomtypetable%num_items)%node => node

end subroutine

function find_atomtype(atomtypetable, elnum, label) result(node)
   type(atomtype_table), intent(in) :: atomtypetable
   integer, intent(in) :: elnum
   integer, intent(in) :: label
   type(tree_node), pointer :: node
   integer :: i

   do i = 1, atomtypetable%num_items
      if (atomtypetable%items(i)%elnum == elnum .and. &
          atomtypetable%items(i)%label == label) then
         node => atomtypetable%items(i)%node
         return
      end if
   end do

   node => null()

end function

subroutine add_typehood(typehoodtable, typehood, node)
   type(typehood_table), intent(inout) :: typehoodtable
   type(tree_node_ptr), intent(in) :: typehood(:)
   type(tree_node), pointer, intent(in) :: node

   typehoodtable%num_items = typehoodtable%num_items + 1
   typehoodtable%items(typehoodtable%num_items)%typehood = typehood
   typehoodtable%items(typehoodtable%num_items)%node => node

end subroutine

function find_typehood(typehoodtable, typehood) result(node)
   type(typehood_table), intent(in) :: typehoodtable
   type(tree_node_ptr), intent(in) :: typehood(:)
   type(tree_node), pointer :: node
   integer :: i

   do i = 1, typehoodtable%num_items
      if (typehoodtable%items(i)%typehood .equiv. typehood) then
         node => typehoodtable%items(i)%node
         return
      end if
   end do

   node => null()

end function

subroutine assign_single_item(leaf)
   type(tree_node), intent(inout) :: leaf
   type(tree_node), pointer :: child
   type(item_node), pointer :: item1, item2, next_item1, next_item2

   item1 => leaf%first_item1
   item2 => leaf%first_item2
   child => add_new_child(leaf)
   next_item1 => item1%next
   next_item2 => item2%next
   call add_item1(child, item1)
   call add_item2(child, item2)
   item1 => next_item1
   item2 => next_item2
   child => add_new_child(leaf)
   do while (associated(item1))
      next_item1 => item1%next
      next_item2 => item2%next
      call add_item1(child, item1)
      call add_item2(child, item2)
      item1 => next_item1
      item2 => next_item2
   end do

   leaf%num_items1 = 0
   leaf%num_items2 = 0
   leaf%first_item1 => null()
   leaf%first_item2 => null()
end subroutine

recursive subroutine reduce_partial_matches(node, assigned)
   type(tree_node), intent(inout) :: node
   logical, intent(inout) :: assigned
   type(tree_node), pointer :: child

   if (associated(node%first_child)) then
      ! Internal node - process children
      child => node%first_child
      do
         call reduce_partial_matches(child, assigned)
         if (assigned) return
         if (.not. associated(child%next_sibling)) exit
         child => child%next_sibling
      end do
   else
      ! Leaf node - assign single item
      if (node%num_items1 == node%num_items2 .and. &
          node%num_items1 > 1) then
         call assign_single_item(node)
         assigned = .true.
      end if
   end if
end subroutine

! Partition atoms by atomic number and label
subroutine compute_eltypes(mol1, mol2, eltree)
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_node), pointer, intent(out) :: eltree
   ! Local variables
   type(tree_node), pointer :: inode
   type(atomtype_table) :: atomtypetable
   integer :: i, elnum, label, num_atoms1, num_atoms2

   num_atoms1 = size(mol1%atoms) 
   num_atoms2 = size(mol2%atoms) 

   eltree => make_new_root(num_atoms1, num_atoms2)
   allocate (atomtypetable%items(num_atoms1 + num_atoms2))
   atomtypetable%num_items = 0

   ! First molecule
   do i = 1, num_atoms1
      elnum = mol1%atoms(i)%elnum
      label = mol1%atoms(i)%label
      inode => find_atomtype(atomtypetable, elnum, label)
      if (.not. associated(inode)) then
         inode => add_new_child(eltree)
         call add_atomtype(atomtypetable, elnum, label, inode)
      end if
      call add_new_item1(inode, i)
   end do

   ! Second molecule
   do i = 1, num_atoms2
      elnum = mol2%atoms(i)%elnum
      label = mol2%atoms(i)%label
      inode => find_atomtype(atomtypetable, elnum, label)
      if (.not. associated(inode)) then
         inode => add_new_child(eltree)
         call add_atomtype(atomtypetable, elnum, label, inode)
      end if
      call add_new_item2(inode, i)
   end do

end subroutine

! Compute Next Level MNA Types
recursive subroutine compute_nextlevelmnatypes(mol1, mol2, itemdir1, itemdir2, inode)
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_node_ptr), dimension(:), intent(in) :: itemdir1, itemdir2
   type(tree_node), intent(inout) :: inode
   ! Local variables
   type(tree_node), pointer :: child

   if (associated(inode%first_child)) then
      ! Internal node - process children
      child => inode%first_child
      do
         call compute_nextlevelmnatypes(mol1, mol2, itemdir1, itemdir2, child)
         if (.not. associated(child%next_sibling)) return
         child => child%next_sibling
      end do
   else
      ! Leaf node
      if (inode%num_items1 + inode%num_items2 > 1) then
         ! Multiple items - update MNA type
         call compute_mnatype(mol1, mol2, itemdir1, itemdir2, inode)
      end if
   end if

end subroutine

subroutine compute_mnatype(mol1, mol2, itemdir1, itemdir2, inode)
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_node_ptr), dimension(:), intent(in) :: itemdir1, itemdir2
   type(tree_node), intent(inout) :: inode
   ! Local variables
   type(tree_node), pointer :: child
   type(item_node), pointer :: item, next_item
   type(tree_node_ptr), allocatable :: typehood(:)
   type(typehood_table) :: typehoodtable

   allocate (typehoodtable%items(inode%num_items1 + inode%num_items2))
   typehoodtable%num_items = 0

   ! First molecule
   item => inode%first_item1
   do while (associated(item))
      typehood = itemdir1(mol1%atoms(item%index)%adjlist)
      child => find_typehood(typehoodtable, typehood)
      if (.not. associated(child)) then
         child => add_new_child(inode)
         call add_typehood(typehoodtable, typehood, child)
      end if
      next_item => item%next
      call add_item1(child, item)
      item => next_item
   end do

   ! Second molecule
   item => inode%first_item2
   do while (associated(item))
      typehood = itemdir2(mol2%atoms(item%index)%adjlist)
      child => find_typehood(typehoodtable, typehood)
      if (.not. associated(child)) then
         child => add_new_child(inode)
         call add_typehood(typehoodtable, typehood, child)
      end if
      next_item => item%next
      call add_item2(child, item)
      item => next_item
   end do

   inode%num_items1 = 0
   inode%num_items2 = 0
   inode%first_item1 => null()
   inode%first_item2 => null()

   ! Revert changes if only child
   if (inode%num_childs == 1) then
      call move_node_items(inode%first_child, inode)
      deallocate (inode%first_child)
      inode%num_childs = 0
   end if

end subroutine

! Iteratively compute MNA types
subroutine compute_consistent_mnatypes(mol1, mol2, mnatree)
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_node), intent(inout) :: mnatree
   ! Local variables
   type(tree_node_ptr), dimension(:), allocatable :: itemdir1, itemdir2

   do

      itemdir1 = mnatree%itemdir1
      itemdir2 = mnatree%itemdir2

      ! Compute MNA upper level types
      call compute_nextlevelmnatypes(mol1, mol2, itemdir1, itemdir2, mnatree)
!      call print_tree(mnatree)

      ! Exit loop if types did not change
      if (all(mnatree%itemdir1 == itemdir1) .and. &
          all(mnatree%itemdir2 == itemdir2)) exit

   end do

end subroutine

end module
