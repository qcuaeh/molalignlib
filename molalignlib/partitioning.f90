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

subroutine unfold_leaves(partition)
   type(type_tree), intent(inout) :: partition
   type(tree_node), pointer :: current, child_node
   type(item_node), pointer :: curr_item1, curr_item2

   if (.not. associated(partition%tree_root%first_child)) return

   current => partition%tree_root%first_child
   do while (associated(current))
      if (associated(current%first_child) .or. current%item_count < 3) then
         current => current%next_sibling
         cycle
      end if

      ! Process leaf node with multiple items
      curr_item1 => current%first_item1
      curr_item2 => current%first_item2
      current%first_item1 => null()
      current%first_item2 => null()
      current%item_count = 0

      do while (associated(curr_item1))
         child_node => add_new_child(current)
         partition%itemdir1(curr_item1%index)%ptr => child_node
         partition%itemdir2(curr_item2%index)%ptr => child_node
         call move_item1(child_node, curr_item1)
         call move_item2(child_node, curr_item2)
      end do

      return
   end do
end subroutine

! Partition atoms by atomic number and label
subroutine compute_eltypes(mol1, mol2, eltypes)
   type(mol_type), intent(in) :: mol1, mol2
   type(type_tree), intent(out) :: eltypes
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
   eltypes%tree_root => make_new_root()

   ! First molecule
   do i = 1, num_atoms1
      elnum = mol1%atoms(i)%elnum
      label = mol1%atoms(i)%label
      inode => find_atomtype(atomtypetable, elnum, label)
      if (.not. associated(inode)) then
         inode => add_new_child(eltypes%tree_root)
         call add_atomtype(atomtypetable, elnum, label, inode)
      end if
      call add_new_item1(inode, i, eltypes%itemdir1)
   end do

   ! Second molecule
   do i = 1, num_atoms2
      elnum = mol2%atoms(i)%elnum
      label = mol2%atoms(i)%label
      inode => find_atomtype(atomtypetable, elnum, label)
      if (.not. associated(inode)) then
         inode => add_new_child(eltypes%tree_root)
         call add_atomtype(atomtypetable, elnum, label, inode)
      end if
      call add_new_item2(inode, i, eltypes%itemdir2)
   end do

end subroutine

! Compute Next Level MNA Types
subroutine compute_nextlevelmnatypes(mol1, mol2, itemdir1, itemdir2, tree_root)
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_node_ptr), dimension(:), intent(in) :: itemdir1, itemdir2
   type(tree_node), intent(inout) :: tree_root

   call traverse_tree(tree_root)

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
      if (inode%item_count <= 1) then
         return
      end if

      call refine_mnatype(mol1, mol2, itemdir1, itemdir2, inode)

   end subroutine

end subroutine

subroutine refine_mnatype(mol1, mol2, itemdir1, itemdir2, inode)
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_node_ptr), dimension(:), intent(in) :: itemdir1, itemdir2
   type(tree_node), intent(inout) :: inode
   ! Local variables
   type(item_node), pointer :: item
   type(tree_node), pointer :: child
   type(tree_node_ptr), allocatable :: typehood(:)
   type(typehood_table) :: typehoodtable

   allocate (typehoodtable%items(inode%item_count))
   typehoodtable%num_items = 0

   ! First molecule
   item => inode%first_item1
   inode%first_item1 => null()
   do while (associated(item))
      typehood = itemdir1(mol1%atoms(item%index)%adjlist)
      child => find_typehood(typehoodtable, typehood)
      if (.not. associated(child)) then
         child => add_new_child(inode)
         call add_typehood(typehoodtable, typehood, child)
      end if
      call move_item1(child, item)
   end do

   ! Second molecule
   item => inode%first_item2
   inode%first_item2 => null()
   do while (associated(item))
      typehood = itemdir2(mol2%atoms(item%index)%adjlist)
      child => find_typehood(typehoodtable, typehood)
      if (.not. associated(child)) then
         child => add_new_child(inode)
         call add_typehood(typehoodtable, typehood, child)
      end if
      call move_item2(child, item)
   end do

   ! Revert changes if only child
   if (inode%child_count == 1) then
      call move_node_items(inode%first_child, inode)
      call prune_branch(inode)
   end if

end subroutine

! Iteratively compute MNA types
subroutine compute_consistent_mnatypes(mol1, mol2, mnatypes)
   type(mol_type), intent(in) :: mol1, mol2
   type(type_tree), intent(inout) :: mnatypes
   ! Local variables
   type(tree_node), pointer :: tree_root
   type(tree_node_ptr), dimension(:), allocatable :: itemdir1, itemdir2

   tree_root => mnatypes%tree_root

   do

      itemdir1 = mnatypes%itemdir1
      itemdir2 = mnatypes%itemdir2

      ! Compute MNA upper level types
      call compute_nextlevelmnatypes(mol1, mol2, itemdir1, itemdir2, tree_root)
!      call print_tree(mnatypes%tree_root)

      ! Exit loop if types did not change
      if (all(mnatypes%itemdir1 == itemdir1) .and. &
          all(mnatypes%itemdir2 == itemdir2)) exit

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
