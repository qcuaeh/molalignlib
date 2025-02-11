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
use tupledict
use partition
use lcrs_tree
use metapartition
use partitiondict
use permutation

implicit none

type :: typehood_item
   type(treenode_ptr), allocatable :: typehood(:)
   type(tree_node), pointer :: node
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

interface find_node
   module procedure atomtypetable_find_item
   module procedure typehoodtable_find_item
end interface

interface add_node
   module procedure atomtypetable_append_item
   module procedure typehoodtable_append_item
end interface

interface operator (.equiv.)
   module procedure typehood_equivalence
end interface

interface operator (.equal.)
   module procedure pointer_equality
end interface

contains

function pointer_equality( array1, array2) result(equal)
   type(treenode_ptr), dimension(:), intent(in) :: array1, array2
   logical :: equal
   integer :: i

   ! Check if sizes are equal
   if (size(array1) /= size(array2)) then
      equal = .false.
      return
   end if

   ! Early exit: Check if any value appears more times in one array
   ! than it does in the other
   do i = 1, size(array1)
      if (.not. associated(array1(i)%ptr, array2(i)%ptr)) then
         equal = .false.
         return
      end if
   end do

   equal = .true.

end function

function typehood_equivalence( array1, array2) result(equiv)
   type(treenode_ptr), dimension(:), intent(in) :: array1, array2
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

subroutine atomtypetable_append_item(atomtypetable, elnum, label, node)
   type(atomtype_table), intent(inout) :: atomtypetable
   integer, intent(in) :: elnum
   integer, intent(in) :: label
   type(tree_node), pointer, intent(in) :: node

   atomtypetable%num_items = atomtypetable%num_items + 1
   atomtypetable%items(atomtypetable%num_items)%elnum = elnum
   atomtypetable%items(atomtypetable%num_items)%label = label
   atomtypetable%items(atomtypetable%num_items)%node => node

end subroutine

function atomtypetable_find_item(atomtypetable, elnum, label) result(node)
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

subroutine typehoodtable_append_item(typehoodtable, typehood, node)
   type(typehood_table), intent(inout) :: typehoodtable
   type(treenode_ptr), intent(in) :: typehood(:)
   type(tree_node), pointer, intent(in) :: node

   typehoodtable%num_items = typehoodtable%num_items + 1
   typehoodtable%items(typehoodtable%num_items)%typehood = typehood
   typehoodtable%items(typehoodtable%num_items)%node => node

end subroutine

function typehoodtable_find_item(typehoodtable, typehood) result(node)
   type(typehood_table), intent(in) :: typehoodtable
   type(treenode_ptr), intent(in) :: typehood(:)
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

subroutine unfold_leaves(partition)
   type(tree_type), intent(inout) :: partition
   type(tree_node), pointer :: current, child_node
   type(item_node), pointer :: curr_item1, curr_item2, next_item1, next_item2

   if (.not. associated(partition%tree%first_child)) return

   current => partition%tree%first_child
   do while (associated(current))
      if (associated(current%first_child) .or. current%item_count1 <= 1) then
         current => current%next_sibling
         cycle
      end if

      ! Process leaf node with multiple items
      curr_item1 => current%first_item1
      curr_item2 => current%first_item2
      current%first_item1 => null()
      current%first_item2 => null()
      current%item_count1 = 0
      current%item_count2 = 0

      do while (associated(curr_item1))
         next_item1 => curr_item1%next_item
         next_item2 => curr_item2%next_item
         curr_item1%next_item => null()
         curr_item2%next_item => null()

         child_node => new_child(current)
         call add_item1(child_node, curr_item1)
         call add_item2(child_node, curr_item2)
         partition%itemdir1(curr_item1%idx)%ptr => child_node
         partition%itemdir2(curr_item2%idx)%ptr => child_node

         curr_item1 => next_item1
         curr_item2 => next_item2
      end do

      return
   end do
end subroutine

! Partition atoms by atomic number and label
subroutine compute_eltypes_tree(mol1, mol2, eltypes)
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_type), intent(out) :: eltypes
   ! Local variables
   type(tree_node), pointer :: inode
   type(atomtype_table) :: atomtypetable
   integer :: i, elnum, label, num_atoms1, num_atoms2

   num_atoms1 = size(mol1%atoms) 
   num_atoms2 = size(mol2%atoms) 
   allocate (eltypes%itemdir1(num_atoms1))
   allocate (eltypes%itemdir2(num_atoms2))
   allocate (atomtypetable%items(num_atoms1 + num_atoms2))
   atomtypetable%num_items = 0
   eltypes%tree => new_tree()

   ! First molecule
   do i = 1, num_atoms1
      elnum = mol1%atoms(i)%elnum
      label = mol1%atoms(i)%label
      inode => find_node(atomtypetable, elnum, label)
      if (.not. associated(inode)) then
         inode => new_child(eltypes%tree)
         call add_node(atomtypetable, elnum, label, inode)
      end if
      call add_new_item1(inode, i)
      eltypes%itemdir1(i)%ptr => inode
   end do

   ! Second molecule
   do i = 1, num_atoms2
      elnum = mol2%atoms(i)%elnum
      label = mol2%atoms(i)%label
      inode => find_node(atomtypetable, elnum, label)
      if (.not. associated(inode)) then
         inode => new_child(eltypes%tree)
         call add_node(atomtypetable, elnum, label, inode)
      end if
      call add_new_item2(inode, i)
      eltypes%itemdir2(i)%ptr => inode
   end do

end subroutine

! Compute Next Level MNA Types
subroutine nextlevel_mnatypes(mol1, mol2, mnatypes)
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_type), intent(inout) :: mnatypes
   ! Local variables
   type(tree_node), pointer :: subtree
   type(treenode_ptr), allocatable :: typehood(:)

   call traverse_tree(mnatypes%tree)
   call delete_tree(subtree)

   contains

   recursive subroutine traverse_tree(inode)
      type(tree_node), intent(inout) :: inode
      type(tree_node), pointer :: child

      ! Internal node, process its children
      if (associated(inode%first_child)) then
         child => inode%first_child
         do
            call traverse_tree(child)
            if (.not. associated(child%next_sibling)) return
            child => child%next_sibling
         end do
      end if

      ! Unique or empty leaf node, do nothing
      if (inode%item_count1 <= 1 .and. inode%item_count2 <= 1) then
         return
      end if

      call next_level_leaf(mol1%atoms, mol2%atoms, mnatypes%itemdir1, mnatypes%itemdir2, inode)

end subroutine

subroutine next_level_leaf(atoms1, atoms2, itemdir1, itemdir2, inode)
   type(atom_type), dimension(:), intent(in) :: atoms1, atoms2
   type(treenode_ptr), dimension(:), intent(inout) :: itemdir1, itemdir2
   type(tree_node), intent(inout) :: inode
   ! Local variables
   type(tree_node), pointer :: child
   type(item_node), pointer :: item, next_item
   type(typehood_table) :: typehoodtable

   subtree => new_tree()
   allocate (typehoodtable%items(inode%item_count1 + inode%item_count2))
   typehoodtable%num_items = 0

   ! First molecule
   item => inode%first_item1
   do while (associated(item))
      next_item => item%next_item
      typehood = itemdir1(atoms1(item%idx)%adjlist)
      child => find_node(typehoodtable, typehood)
      if (.not. associated(child)) then
         child => new_child(subtree)
         call add_node(typehoodtable, typehood, child)
      end if
      call add_item1(child, item)
      item => next_item
   end do

   ! Second molecule
   item => inode%first_item2
   do while (associated(item))
      next_item => item%next_item
      typehood = itemdir2(atoms2(item%idx)%adjlist)
      child => find_node(typehoodtable, typehood)
      if (.not. associated(child)) then
         child => new_child(subtree)
         call add_node(typehoodtable, typehood, child)
      end if
      call add_item2(child, item)
      item => next_item
   end do

   ! Attach subtree to tree if leaf is branched
   if (associated(subtree%first_child%next_sibling)) then
      call update_itemdir(itemdir1, itemdir2, subtree)
      inode%first_item1 => null()
      inode%first_item2 => null()
      inode%first_child => subtree%first_child
      subtree%first_child => null()
   else
      subtree%first_child%first_item1 => null()
      subtree%first_child%first_item2 => null()
   end if

   call delete_tree(subtree)

   end subroutine

end subroutine

! Iteratively compute MNA types
subroutine compute_mnatypes_tree(mol1, mol2, mnatypes)
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_type), intent(inout) :: mnatypes
   ! Local variables
   type(treenode_ptr), dimension(:), allocatable :: itemdir1, itemdir2

   allocate (itemdir1(size(mol1%atoms)))
   allocate (itemdir2(size(mol2%atoms)))

   do

!      call mnatypes%print_tree()
      itemdir1 = mnatypes%itemdir1
      itemdir2 = mnatypes%itemdir2

      ! Compute MNA upper level types
      call nextlevel_mnatypes(mol1, mol2, mnatypes)

      ! Exit loop if types did not change
      if ((mnatypes%itemdir1 .equal. itemdir1) .and. &
          (mnatypes%itemdir2 .equal. itemdir2)) exit

   end do

end subroutine

! Level up MNA types
subroutine levelup_mnatypes(mol, mnatypes, subtypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(in) :: mnatypes
   type(partition_type), intent(out) :: subtypes
   ! Local variables
   integer :: h, i, iatom
   type(tupledict_type) :: typedict
   type(partpointer_type), allocatable :: typelist(:)
   integer, allocatable :: typehood(:)

   call subtypes%initialize(mnatypes%num_items)
   call typedict%initialize(mnatypes%largest_part_size, 'unordered')
   allocate (typelist(typedict%num_slots))

   do h = 1, mnatypes%num_parts

      do i = 1, mnatypes%parts(h)%part_size
         iatom = mnatypes%parts(h)%items(i)
         typehood = mnatypes%idcs(mol%atoms(iatom)%adjlist)
         if (.not. (typehood .in. typedict)) then
            typelist(typedict%new_index(typehood))%ptr => &
               subtypes%new_part(mnatypes%parts(h)%part_size)
         end if
         call typelist(typedict%get_index(typehood))%ptr%add(iatom)
      end do

      call typedict%reset()

   end do

end subroutine

! Iteratively compute MNA types
subroutine compute_mnatypes(mol, mnatypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(inout) :: mnatypes
   ! Local variables
   type(partition_type) :: subtypes

   do

!      write (stderr, *)
!      call mnatypes%print_parts()

      ! Compute MNA upper level types
      call levelup_mnatypes(mol, mnatypes, subtypes)

      ! Exit loop if types did not change
      if (subtypes == mnatypes) then
         mnatypes = subtypes
         exit
      end if

      ! Update mnatypes
      mnatypes = subtypes

   end do

end subroutine

end module
