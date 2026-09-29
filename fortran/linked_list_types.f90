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

module linked_list_types
! Linked-list (left-child right-sibling) representation of the HNA
! partition tree and the assignment tree, used while they are built, when
! parts are subdivided dynamically. See indexed_list_types for the array
! form used by the searches.
!
! Every node type keeps a global index, assigned from counters shared by
! all nodes of a tree (allocated by the root, pointed to by the others),
! which becomes its position in the array form.
use parameters
use common_types
implicit none
private
public address
public isdescendant
public new_root_part
public new_child_part
public new_bare_link
public new_chain_link
public new_root_chain
public new_child_chain
public chain_from_partition
public link_part
public update_itemdir
public add_new_item1
public add_new_item2
public move_first_item1
public move_first_item2
public copy_part_items
public move_part_items
public link_to_partition
public find_child_part
public add_branch_part
public delete_chain
public delete_part
public delete_part_tree
public print_part_tree
public print_part_items
public print_tree_items
public print_leaf_items
public print_link_itemdir
public print_part_signature
public print_tree_signatures
public print_part_indices
public print_chain_indices
public print_chain_tree
public is_partition_uneven
public print_partition_chain
public operator(==)
public operator(.equiv.)

! Atom of a part
type, public :: item_node_t
   integer(ik) :: idx
   integer(ik) :: global_idx
   type(item_node_t), pointer :: next_item
end type

! Link of a chain: one refinement level, listing its parts and the vertex
! directory (part of each atom placed at that level, null otherwise)
type, public :: chain_node_t
   integer(ik) :: n_parts
   integer(ik) :: global_idx
   integer(ik), pointer :: total_partrefs => null()
   type(chain_node_t), pointer :: next_link
   type(partref_node_t), pointer :: first_partref
   type(partref_node_t), pointer :: last_partref
   type(part_nodeptr_t), dimension(:), pointer :: itemdir1
   type(part_nodeptr_t), dimension(:), pointer :: itemdir2
end type

! Chain: assignment tree node, or the sequence of refinement levels of a
! partition. split_part is the part individualized at the node (null for
! the root).
type, public :: chaintree_node_t
   integer(ik) :: n_links
   integer(ik) :: n_children
   integer(ik) :: global_idx
   integer(ik), pointer :: n_atoms1 => null()
   integer(ik), pointer :: n_atoms2 => null()
   integer(ik), pointer :: total_chains => null()
   integer(ik), pointer :: total_links => null()
   integer(ik), pointer :: total_partrefs => null()
   type(chaintree_node_t), pointer :: parent_chain
   type(chaintree_node_t), pointer :: first_child_chain
   type(chaintree_node_t), pointer :: last_child_chain
   type(chaintree_node_t), pointer :: next_sibling_chain
   type(chain_node_t), pointer :: first_link
   type(chain_node_t), pointer :: last_link
   type(partition_node_t), pointer :: split_part
end type

! Part of the HNA partition tree
type, public :: partition_node_t
   integer(ik) :: depth
   integer(ik) :: global_idx
   integer(ik) :: n_items1
   integer(ik) :: n_items2
   integer(ik) :: n_children
   integer(ik), pointer :: total_parts => null()
   integer(ik), pointer :: total_items1 => null()
   integer(ik), pointer :: total_items2 => null()
   type(item_node_t), pointer :: first_item1
   type(item_node_t), pointer :: first_item2
   type(item_node_t), pointer :: last_item1
   type(item_node_t), pointer :: last_item2
   type(partition_node_t), pointer :: parent_part
   type(partition_node_t), pointer :: first_child_part
   type(partition_node_t), pointer :: last_child_part
   type(partition_node_t), pointer :: next_sibling_part
   ! Bond-typed signature: one edge_code(part global_idx, bond type)
   ! per neighbor with a known part, compared as a multiset
   integer(ik), dimension(:), pointer :: signature
end type

! Reference to a part in a link
type, public :: partref_node_t
   integer(ik) :: global_idx
   type(partition_node_t), pointer :: part
   type(partref_node_t), pointer :: nextref
end type

! Part node pointer
type, public :: part_nodeptr_t
   type(partition_node_t), pointer :: ptr
end type

interface address
   module procedure address_part
end interface

interface operator(==)
   module procedure part_nodeptr_equality
end interface

interface operator (.equiv.)
   module procedure signature_equivalence
end interface

contains

character(4) function address_part(nodeptr) result(address)
! Last four hex digits of the address of a part, to identify it in
! debugging output
   use iso_c_binding, only: c_loc, c_intptr_t
   type(partition_node_t), pointer, intent(in) :: nodeptr
   if (associated(nodeptr)) then
      write (address, '(Z4.4)') modulo(transfer(c_loc(nodeptr), c_intptr_t), 16**4)
   else
      address = 'NULL'
   end if
end function

elemental function part_nodeptr_equality(left, right) result(equality)
   type(part_nodeptr_t), intent(in) :: left, right
   logical(lk) :: equality
   equality = associated(left%ptr, right%ptr)
end function

function signature_equivalence(signature1, signature2) result(equiv)
! Multiset equality of two integer-encoded signatures
   integer(ik), dimension(:), intent(in) :: signature1, signature2
   logical(lk) :: equiv
   integer(ik) :: matches1, matches2
   integer(ik) :: i, j

   if (size(signature1) /= size(signature2)) then
      equiv = .FALSE.
      return
   end if

   do i = 1, size(signature1)
      matches1 = 0
      matches2 = 0
      do j = 1, size(signature1)
         if (signature1(i) == signature2(j)) matches1 = matches1 + 1
         if (signature1(i) == signature1(j)) matches2 = matches2 + 1
      end do
      if (matches1 /= matches2) then
         equiv = .FALSE.
         return
      end if
   end do

   equiv = .TRUE.
end function

logical(lk) function isdescendant(part, top_part)
! Whether part lies in the subtree of top_part (excluding top_part itself)
   type(partition_node_t), pointer, intent(in) :: part, top_part
   ! Local variables
   type(partition_node_t), pointer :: up_part

   if (part%depth <= top_part%depth) then
      isdescendant = .FALSE.
      return
   end if

   up_part => part%parent_part
   do while (up_part%depth > top_part%depth)
      up_part => up_part%parent_part
   end do

   if (associated(up_part, top_part)) then
      isdescendant = .TRUE.
   else
      isdescendant = .FALSE.
   end if
end function

function new_bare_part() result(part)
   type(partition_node_t), pointer :: part

   allocate (part)

   part%global_idx = 0  ! Set when added to a tree
   part%n_items1 = 0
   part%n_items2 = 0
   part%n_children = 0
   part%total_parts => null()
   part%total_items1 => null()
   part%total_items2 => null()
   part%parent_part => null()
   part%first_item1 => null()
   part%last_item1 => null()
   part%first_item2 => null()
   part%last_item2 => null()
   part%first_child_part => null()
   part%last_child_part => null()
   part%next_sibling_part => null()
   allocate (part%signature(0))
end function

function new_root_part() result(part)
   type(partition_node_t), pointer :: part

   part => new_bare_part()
   part%depth = 0
   part%global_idx = 1

   ! Counters shared by the whole tree
   allocate(part%total_parts)
   allocate(part%total_items1)
   allocate(part%total_items2)

   part%total_parts = 1
   part%total_items1 = 0
   part%total_items2 = 0
end function

function new_child_part(parent_part) result(child_part)
   type(partition_node_t), pointer, intent(inout) :: parent_part
   type(partition_node_t), pointer :: child_part

   child_part => new_bare_part()
   child_part%depth = parent_part%depth + 1
   child_part%parent_part => parent_part
   child_part%next_sibling_part => null()

   ! Share the tree counters
   child_part%total_parts => parent_part%total_parts
   child_part%total_items1 => parent_part%total_items1
   child_part%total_items2 => parent_part%total_items2

   child_part%total_parts = child_part%total_parts + 1
   child_part%global_idx = child_part%total_parts

   ! Append to the children of the parent (the part is not linked)
   if (.not. associated(parent_part%first_child_part)) then
      parent_part%first_child_part => child_part
   else
      parent_part%last_child_part%next_sibling_part => child_part
   end if
   parent_part%last_child_part => child_part
   parent_part%n_children = parent_part%n_children + 1
end function

subroutine link_part(link, part)
! Append a reference to part to link (the vertex directory is not updated)
   type(chain_node_t), target, intent(inout) :: link
   type(partition_node_t), target, intent(inout) :: part
   type(partref_node_t), pointer :: newref

   allocate(newref)
   newref%part => part
   newref%nextref => null()

   link%total_partrefs = link%total_partrefs + 1
   newref%global_idx = link%total_partrefs

   if (.not. associated(link%first_partref)) then
      link%first_partref => newref
   else
      link%last_partref%nextref => newref
   end if
   link%last_partref => newref
   link%n_parts = link%n_parts + 1
end subroutine

subroutine add_new_item1(part, idx)
! Append atom idx of molecule 1 to part (no vertex directory is updated)
   type(partition_node_t), target, intent(inout) :: part
   integer(ik), intent(in) :: idx
   type(item_node_t), pointer :: new_item

   allocate(new_item)
   new_item%idx = idx
   new_item%next_item => null()

   part%total_items1 = part%total_items1 + 1
   new_item%global_idx = part%total_items1

   if (.not. associated(part%first_item1)) then
      part%first_item1 => new_item
   else
      part%last_item1%next_item => new_item
   end if
   part%last_item1 => new_item
   part%n_items1 = part%n_items1 + 1
end subroutine

subroutine add_new_item2(part, idx)
! Append atom idx of molecule 2 to part (no vertex directory is updated)
   type(partition_node_t), target, intent(inout) :: part
   integer(ik), intent(in) :: idx
   type(item_node_t), pointer :: new_item

   allocate(new_item)
   new_item%idx = idx
   new_item%next_item => null()

   part%total_items2 = part%total_items2 + 1
   new_item%global_idx = part%total_items2

   if (.not. associated(part%first_item2)) then
      part%first_item2 => new_item
   else
      part%last_item2%next_item => new_item
   end if
   part%last_item2 => new_item
   part%n_items2 = part%n_items2 + 1
end subroutine

subroutine copy_part_items(orig, dest)
! Append copies of all atoms of orig to dest (no vertex directory is updated)
   type(partition_node_t), intent(inout) :: orig, dest
   type(item_node_t), pointer :: item

   item => orig%first_item1
   do while (associated(item))
      call add_new_item1(dest, item%idx)
      item => item%next_item
   end do

   item => orig%first_item2
   do while (associated(item))
      call add_new_item2(dest, item%idx)
      item => item%next_item
   end do
end subroutine

subroutine move_first_item1(orig, dest)
! Move the first atom of molecule 1 from orig to the end of dest (no vertex
! directory is updated)
   type(partition_node_t), intent(inout) :: orig
   type(partition_node_t), target, intent(inout) :: dest
   type(item_node_t), pointer :: item_to_move

   item_to_move => orig%first_item1
   if ((.not. associated(item_to_move))) error stop

   orig%first_item1 => item_to_move%next_item
   if (.not. associated(orig%first_item1)) then
      orig%last_item1 => null()
   end if
   orig%n_items1 = orig%n_items1 - 1

   item_to_move%next_item => null()

   if (.not. associated(dest%first_item1)) then
      dest%first_item1 => item_to_move
   else
      dest%last_item1%next_item => item_to_move
   end if
   dest%last_item1 => item_to_move
   dest%n_items1 = dest%n_items1 + 1
end subroutine

subroutine move_first_item2(orig, dest)
! Move the first atom of molecule 2 from orig to the end of dest (no vertex
! directory is updated)
   type(partition_node_t), intent(inout) :: orig
   type(partition_node_t), target, intent(inout) :: dest
   type(item_node_t), pointer :: item_to_move

   item_to_move => orig%first_item2
   if ((.not. associated(item_to_move))) error stop

   orig%first_item2 => item_to_move%next_item
   if (.not. associated(orig%first_item2)) then
      orig%last_item2 => null()
   end if
   orig%n_items2 = orig%n_items2 - 1

   item_to_move%next_item => null()

   if (.not. associated(dest%first_item2)) then
      dest%first_item2 => item_to_move
   else
      dest%last_item2%next_item => item_to_move
   end if
   dest%last_item2 => item_to_move
   dest%n_items2 = dest%n_items2 + 1
end subroutine

subroutine move_part_items(orig, dest)
! Move all atoms of orig to dest (no vertex directory is updated)
   type(partition_node_t), intent(inout) :: orig, dest

   do while (associated(orig%first_item1))
      call move_first_item1(orig, dest)
   end do

   do while (associated(orig%first_item2))
      call move_first_item2(orig, dest)
   end do
end subroutine

subroutine delete_chain(root_chain)
! Delete a chain and its links; the referenced parts are preserved
   type(chaintree_node_t), pointer, intent(inout) :: root_chain
   type(chain_node_t), pointer :: link, next_link

   if ((.not. associated(root_chain))) error stop

   link => root_chain%first_link
   do while (associated(link))
      next_link => link%next_link
      call delete_link(link)
      link => next_link
   end do

   deallocate(root_chain)
   root_chain => null()
end subroutine

subroutine delete_link(link)
! Delete a link and its part references, but not the parts
   type(chain_node_t), pointer, intent(inout) :: link
   type(partref_node_t), pointer :: partref, nextref

   if ((.not. associated(link))) error stop

   partref => link%first_partref
   do while (associated(partref))
      nextref => partref%nextref
      deallocate(partref)
      partref => nextref
   end do

   deallocate(link%itemdir1)
   deallocate(link%itemdir2)
   deallocate(link)
   link => null()
end subroutine

subroutine delete_part_tree(partition_tree)
! Delete a whole part tree and its counters
   type(partition_node_t), pointer, intent(inout) :: partition_tree

   if (.not. associated(partition_tree)) error stop

   call delete_part_children(partition_tree)

   deallocate(partition_tree%total_parts)
   deallocate(partition_tree%total_items1)
   deallocate(partition_tree%total_items2)

   call delete_part(partition_tree)
end subroutine

recursive subroutine delete_part_children(parent_part)
! Delete all descendants of a part
   type(partition_node_t), pointer, intent(in) :: parent_part
   type(partition_node_t), pointer :: child_part, next_child

   if (.not. associated(parent_part)) error stop

   child_part => parent_part%first_child_part
   do while (associated(child_part))
      next_child => child_part%next_sibling_part
      call delete_part_children(child_part)
      call delete_part(child_part)

      child_part => next_child
   end do
end subroutine

subroutine delete_part(part_node)
! Delete a single part and its atoms (its children must be deleted already)
   type(partition_node_t), pointer, intent(inout) :: part_node

   if ((.not. associated(part_node))) error stop

   call deallocate_items(part_node%first_item1)
   call deallocate_items(part_node%first_item2)

   if (associated(part_node%signature)) then
      deallocate(part_node%signature)
   end if

   deallocate(part_node)
   part_node => null()
end subroutine

subroutine deallocate_items(first_item)
   type(item_node_t), pointer, intent(inout) :: first_item
   type(item_node_t), pointer :: item, next_item

   item => first_item
   do while (associated(item))
      next_item => item%next_item
      deallocate(item)
      item => next_item
   end do

   first_item => null()
end subroutine

subroutine link_to_partition(link, partition)
! Flat partition (partition_t) with the parts of link
   type(chain_node_t), target, intent(in) :: link
   type(partition_t), intent(out) :: partition
   type(partref_node_t), pointer :: partref
   type(item_node_t), pointer :: item
   integer(ik) :: i, j

   partition%n_parts = link%n_parts
   allocate(partition%parts(partition%n_parts))
   allocate(partition%itemdir1(size(link%itemdir1)))
   allocate(partition%itemdir2(size(link%itemdir2)))

   partref => link%first_partref
   i = 1
   do while (associated(partref))
      if (.not. associated(partref%part)) error stop

      partition%parts(i)%n_items1 = partref%part%n_items1
      partition%parts(i)%n_items2 = partref%part%n_items2
      allocate(partition%parts(i)%items1(partref%part%n_items1))
      allocate(partition%parts(i)%items2(partref%part%n_items2))

      item => partref%part%first_item1
      do j = 1, partref%part%n_items1
         partition%parts(i)%items1(j) = item%idx
         partition%itemdir1(item%idx) = i
         item => item%next_item
      end do

      item => partref%part%first_item2
      do j = 1, partref%part%n_items2
         partition%parts(i)%items2(j) = item%idx
         partition%itemdir2(item%idx) = i
         item => item%next_item
      end do

      partref => partref%nextref
      i = i + 1
   end do
end subroutine

subroutine print_part_items(part)
   type(partition_node_t), pointer, intent(in) :: part
   type(item_node_t), pointer :: item

   item => part%first_item1
   do while (associated(item))
      write(stderr, '(1X,I0)', advance='no') item%idx
      item => item%next_item
   end do

   write(stderr, '(A)', advance='no') ' /'

   item => part%first_item2
   do while (associated(item))
      write(stderr, '(1X,I0)', advance='no') item%idx
      item => item%next_item
   end do

   write(stderr, *)
end subroutine

subroutine print_link_itemdir(link)
   type(chain_node_t), intent(in) :: link
   integer(ik) :: i

   write(stderr,'(A)') "Item Directory 1:"
   do i = 1, size(link%itemdir1)
      write(stderr,'(2X,I0,1X,A,1X,A)') i, "->", address(link%itemdir1(i)%ptr)
   end do

   write(stderr,'(A)') "Item Directory 2:"
   do i = 1, size(link%itemdir2)
      write(stderr,'(2X,I0,1X,A,1X,A)') i, "->", address(link%itemdir2(i)%ptr)
   end do

   write(stderr, *)
end subroutine

function new_bare_chain() result(chain)
   type(chaintree_node_t), pointer :: chain

   allocate(chain)

   chain%n_links = 0
   chain%n_children = 0
   chain%global_idx = 0
   chain%n_atoms1 => null()
   chain%n_atoms2 => null()
   chain%total_chains => null()
   chain%total_links => null()
   chain%total_partrefs => null()
   chain%first_child_chain => null()
   chain%last_child_chain => null()
   chain%next_sibling_chain => null()
   chain%first_link => null()
   chain%last_link => null()
end function

function new_root_chain(n_atoms1, n_atoms2) result(chain)
   integer(ik), intent(in) :: n_atoms1, n_atoms2
   type(chaintree_node_t), pointer :: chain

   chain => new_bare_chain()
   chain%parent_chain => null()
   chain%split_part => null()
   chain%global_idx = 1

   ! Counters shared by the whole tree
   allocate(chain%n_atoms1)
   allocate(chain%n_atoms2)
   allocate(chain%total_chains)
   allocate(chain%total_links)
   allocate(chain%total_partrefs)
   chain%n_atoms1 = n_atoms1
   chain%n_atoms2 = n_atoms2
   chain%total_chains = 1
   chain%total_links = 0
   chain%total_partrefs = 0
end function

function new_bare_link() result(link)
   type(chain_node_t), pointer :: link

   allocate(link)
   link%n_parts = 0
   link%global_idx = 0  ! Set when added to a chain
   link%total_partrefs => null()
   link%first_partref => null()
   link%last_partref => null()
   link%next_link => null()
end function

function new_chain_link(chain) result(link)
! Append a new link, with an empty vertex directory, to chain
   type(chaintree_node_t), target, intent(inout) :: chain
   type(chain_node_t), pointer :: link
   integer(ik) :: i

   link => new_bare_link()
   link%total_partrefs => chain%total_partrefs

   chain%total_links = chain%total_links + 1
   link%global_idx = chain%total_links

   allocate(link%itemdir1(chain%n_atoms1))
   allocate(link%itemdir2(chain%n_atoms2))
   do i = 1, chain%n_atoms1
      link%itemdir1(i)%ptr => null()
   end do
   do i = 1, chain%n_atoms2
      link%itemdir2(i)%ptr => null()
   end do

   if (.not. associated(chain%first_link)) then
      chain%first_link => link
   else
      chain%last_link%next_link => link
   end if
   chain%last_link => link
   chain%n_links = chain%n_links + 1
end function

function find_child_part(part, signature) result(child_part)
! Child of part with the given signature, or null if there is none
   type(partition_node_t), intent(in) :: part
   integer(ik), dimension(:), intent(in) :: signature
   type(partition_node_t), pointer :: child_part

   child_part => part%first_child_part
   do while (associated(child_part))
      if (child_part%signature .equiv. signature) then
         return
      end if
      child_part => child_part%next_sibling_part
   end do

   child_part => null()
end function

subroutine print_part_tree(partition_tree)
   type(partition_node_t), pointer, intent(in) :: partition_tree
   logical(lk), dimension(:), allocatable :: is_last_child

   if (.not. associated(partition_tree)) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "       Part Tree"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Allocate tracking array for tree lines (max depth 100)
   allocate(is_last_child(100))
   is_last_child = .FALSE.

   write(stderr, '(A)') 'ROOT'
   call print_part_recurse(partition_tree, 0, is_last_child)
   write(stderr, *)

   deallocate(is_last_child)
end subroutine

function chain_from_partition(partition) result(chain)
! Chain whose first link holds the parts of partition (typically the atom
! types), as children of a new part tree, with a complete vertex directory
   type(partition_t), intent(in) :: partition
   ! Local variables
   type(chaintree_node_t), pointer :: chain
   type(partition_node_t), pointer :: partition_tree
   type(chain_node_t), pointer :: first_link
   type(partition_node_t), pointer :: new_part
   integer(ik) :: i, j

   partition_tree => new_root_part()
   chain => new_root_chain(size(partition%itemdir1), size(partition%itemdir2))
   first_link => new_chain_link(chain)

   do i = 1, partition%n_parts
      new_part => new_child_part(partition_tree)

      do j = 1, partition%parts(i)%n_items1
         call add_new_item1(new_part, partition%parts(i)%items1(j))
      end do

      do j = 1, partition%parts(i)%n_items2
         call add_new_item2(new_part, partition%parts(i)%items2(j))
      end do

      call link_part(first_link, new_part)
      call update_itemdir(first_link, new_part)
   end do
end function

function new_child_chain(chain, split_part) result(new_chain)
! New assignment tree node below chain that individualizes split_part
   type(chaintree_node_t), pointer, intent(inout) :: chain
   type(partition_node_t), pointer, intent(in) :: split_part
   type(chaintree_node_t), pointer :: new_chain

   new_chain => new_bare_chain()
   new_chain%split_part => split_part
   new_chain%parent_chain => chain

   ! Share the tree counters
   new_chain%n_atoms1 => chain%n_atoms1
   new_chain%n_atoms2 => chain%n_atoms2
   new_chain%total_chains => chain%total_chains
   new_chain%total_links => chain%total_links
   new_chain%total_partrefs => chain%total_partrefs

   new_chain%total_chains = new_chain%total_chains + 1
   new_chain%global_idx = new_chain%total_chains

   ! Append to the children of the parent
   if (.not. associated(chain%first_child_chain)) then
      chain%first_child_chain => new_chain
      chain%last_child_chain => new_chain
   else
      chain%last_child_chain%next_sibling_chain => new_chain
      chain%last_child_chain => new_chain
   end if

   chain%n_children = chain%n_children + 1
end function

subroutine add_branch_part(link, part)
! Insert a reference to part in link, keeping the list sorted by increasing
! size, so that the smallest parts are individualized first. Parts with a
! single pair need no individualization and are skipped. Used for
! candidate lists only, so the global partref counter is not updated.
   type(chain_node_t), target, intent(inout) :: link
   type(partition_node_t), target, intent(in) :: part
   type(partref_node_t), pointer :: partref, prevref, newref

   if (part%n_items1 < 2) return

   allocate(newref)
   newref%part => part
   newref%nextref => null()

   if (.not. associated(link%first_partref)) then
      link%first_partref => newref
      link%last_partref => newref
      link%n_parts = link%n_parts + 1
      return
   end if

   ! Insert before the first part that is not smaller
   partref => link%first_partref
   prevref => null()

   do while (associated(partref))
      if (part%n_items1 <= partref%part%n_items1) then
         exit
      end if
      prevref => partref
      partref => partref%nextref
   end do

   newref%nextref => partref

   if (associated(prevref)) then
      prevref%nextref => newref
      if (.not. associated(partref)) then
         link%last_partref => newref
      end if
   else
      link%first_partref => newref
      if (.not. associated(partref)) then
         link%last_partref => newref
      end if
   end if

   link%n_parts = link%n_parts + 1
end subroutine

subroutine update_itemdir(link, part)
! Point the vertex directory of link to part for all atoms of part
   type(chain_node_t), target, intent(inout) :: link
   type(partition_node_t), target, intent(in) :: part
   type(item_node_t), pointer :: item

   item => part%first_item1
   do while (associated(item))
      link%itemdir1(item%idx)%ptr => part
      item => item%next_item
   end do

   item => part%first_item2
   do while (associated(item))
      link%itemdir2(item%idx)%ptr => part
      item => item%next_item
   end do
end subroutine

subroutine print_chain_tree(assignment_tree)
   type(chaintree_node_t), pointer, intent(in) :: assignment_tree
   logical(lk), dimension(:), allocatable :: is_last_child

   if (.not. associated(assignment_tree)) error stop

   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "     Assignment Tree"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Allocate tracking array for tree lines (max depth 100)
   allocate(is_last_child(100))
   is_last_child = .FALSE.

   write(stderr, '(A)') 'ROOT'
   call print_chain_recurse(assignment_tree, 0, is_last_child)
   write(stderr, *)

   deallocate(is_last_child)
end subroutine

subroutine print_part_signature(signature)
! Print each signature entry as part_index:bond_type
   integer(ik), dimension(:), intent(in) :: signature
   integer(ik) :: i

   do i = 1, size(signature)
      write(stderr,'(1X,I0,A,I0)',advance='no') signature(i)/BOND_TYPE_RADIX, ':', &
            modulo(signature(i), BOND_TYPE_RADIX)
   end do
   write(stderr,*)
end subroutine

subroutine print_tree_signatures(partition_tree)
! Print the signatures of all parts of a part tree
   type(partition_node_t), pointer, intent(in) :: partition_tree

   if (.not. associated(partition_tree)) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "    Part Signatures"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   call print_signature_recurse(partition_tree)
   write(stderr, *)
end subroutine

recursive subroutine print_signature_recurse(part)
   type(partition_node_t), pointer, intent(in) :: part
   type(partition_node_t), pointer :: child_part

   if (.not. associated(part)) return

   child_part => part%first_child_part
   do while (associated(child_part))
      write(stderr,'(A,I0,A)',advance='no') 'Part ', child_part%global_idx, ':'
      call print_part_signature(child_part%signature)

      call print_signature_recurse(child_part)

      child_part => child_part%next_sibling_part
   end do
end subroutine

subroutine print_tree_items(partition_tree)
   type(partition_node_t), pointer, intent(in) :: partition_tree

   if (.not. associated(partition_tree)) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "      Part Items"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   call print_items_recurse(partition_tree)
   write(stderr, *)
end subroutine

recursive subroutine print_items_recurse(part)
   type(partition_node_t), pointer, intent(in) :: part
   type(partition_node_t), pointer :: child_part

   if (.not. associated(part)) return

   child_part => part%first_child_part
   do while (associated(child_part))
      write(stderr, '(A)', advance='no') address(child_part) // ':'
      call print_part_items(child_part)

      call print_items_recurse(child_part)

      child_part => child_part%next_sibling_part
   end do
end subroutine

subroutine print_leaf_items(partition_tree)
   type(partition_node_t), pointer, intent(in) :: partition_tree

   if (.not. associated(partition_tree)) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "      Leaf Items"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   call print_leaf_items_recurse(partition_tree)
   write(stderr, *)
end subroutine

recursive subroutine print_leaf_items_recurse(part)
   type(partition_node_t), pointer, intent(in) :: part
   type(partition_node_t), pointer :: child_part

   if (.not. associated(part)) return

   child_part => part%first_child_part
   do while (associated(child_part))
      if (child_part%n_children == 0) then
         write(stderr, '(A)', advance='no') address(child_part) // ':'
         call print_part_items(child_part)
      end if

      call print_leaf_items_recurse(child_part)

      child_part => child_part%next_sibling_part
   end do
end subroutine

subroutine print_part_indices(partition_tree)
! Print the global indices of all parts and atoms of a part tree
   type(partition_node_t), pointer, intent(in) :: partition_tree

   if (.not. associated(partition_tree)) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, '(A)') repeat("=", 35)
   write(stderr, '(A)') "        Part Indices"
   write(stderr, '(A)') repeat("=", 35)
   write(stderr, *)

   write(stderr, '(A,I0)') "Total parts created: ", partition_tree%total_parts
   write(stderr, '(A,I0)') "Total items1 created: ", partition_tree%total_items1
   write(stderr, '(A,I0)') "Total items2 created: ", partition_tree%total_items2
   write(stderr, *)

   call print_part_indices_recurse(partition_tree)
   write(stderr, *)
end subroutine

recursive subroutine print_part_indices_recurse(part)
   type(partition_node_t), pointer, intent(in) :: part
   type(partition_node_t), pointer :: child_part
   type(item_node_t), pointer :: item

   if (.not. associated(part)) return

   child_part => part%first_child_part
   do while (associated(child_part))
      write(stderr, '(A,I0,A,A,A)') "Part ", child_part%global_idx, " (address: ", address(child_part), ")"

      item => child_part%first_item1
      do while (associated(item))
         write(stderr, '(A,I0,A,I0,A)') "  Item1 ", item%global_idx, " (idx: ", item%idx, ")"
         item => item%next_item
      end do

      item => child_part%first_item2
      do while (associated(item))
         write(stderr, '(A,I0,A,I0,A)') "  Item2 ", item%global_idx, " (idx: ", item%idx, ")"
         item => item%next_item
      end do

      call print_part_indices_recurse(child_part)

      child_part => child_part%next_sibling_part
   end do
end subroutine

subroutine print_chain_indices(assignment_tree)
! Print the global indices of all nodes, links and part references of an
! assignment tree
   type(chaintree_node_t), pointer, intent(in) :: assignment_tree

   if (.not. associated(assignment_tree)) then
      write(stderr, '(A)') "Chain tree is empty"
      return
   end if

   write(stderr, '(A)') repeat("=", 35)
   write(stderr, '(A)') "     Chain Indices"
   write(stderr, '(A)') repeat("=", 35)
   write(stderr, *)

   write(stderr, '(A,I0)') "Total chains created: ", assignment_tree%total_chains
   write(stderr, '(A,I0)') "Total links created: ", assignment_tree%total_links
   write(stderr, '(A,I0)') "Total partrefs created: ", assignment_tree%total_partrefs
   write(stderr, *)

   call print_chain_indices_recurse(assignment_tree)
   write(stderr, *)
end subroutine

recursive subroutine print_chain_indices_recurse(chain)
   type(chaintree_node_t), pointer, intent(in) :: chain
   type(chaintree_node_t), pointer :: child_chain
   type(chain_node_t), pointer :: link
   type(partref_node_t), pointer :: partref

   if (.not. associated(chain)) return

   link => chain%first_link
   do while (associated(link))
      write(stderr, '(A,I0,A,I0,A)') "  Link ", link%global_idx, " (", link%n_parts, " parts)"

      partref => link%first_partref
      do while (associated(partref))
         write(stderr, '(A,I0,A,A,A)') "    Partref ", partref%global_idx, " (part: ", address(partref%part), ")"
         partref => partref%nextref
      end do

      link => link%next_link
   end do

   child_chain => chain%first_child_chain
   do while (associated(child_chain))
      write(stderr, '(A,I0,A,A,A)') "Chain ", child_chain%global_idx, " (split part: ", address(child_chain%split_part), ")"

      call print_chain_indices_recurse(child_chain)

      child_chain => child_chain%next_sibling_chain
   end do
end subroutine

recursive subroutine print_part_recurse(part, depth, is_last_child)
   type(partition_node_t), pointer, intent(in) :: part
   integer(ik), intent(in) :: depth
   logical(lk), dimension(:), intent(inout) :: is_last_child
   type(partition_node_t), pointer :: child_part, next_child
   integer(ik) :: i

   if (.not. associated(part)) return

   child_part => part%first_child_part
   do while (associated(child_part))
      ! Tree drawing prefix
      next_child => child_part%next_sibling_part
      is_last_child(depth + 1) = .not. associated(next_child)
      do i = 1, depth
         if (is_last_child(i)) then
            write(stderr, '(A)', advance='no') "    "
         else
            write(stderr, '(A)', advance='no') "|   "
         end if
      end do

      if (is_last_child(depth + 1)) then
         write(stderr, '(A)', advance='no') "`--"
      else
         write(stderr, '(A)', advance='no') "|--"
      end if

      write(stderr, '(A,1X,A,I0,A,I0,A)') address(child_part), &
         '(', child_part%n_items1, '/', child_part%n_items2, ')'

      call print_part_recurse(child_part, depth + 1, is_last_child)

      child_part => next_child
   end do
end subroutine

recursive subroutine print_chain_recurse(chain, depth, is_last_child)
   type(chaintree_node_t), pointer, intent(in) :: chain
   integer(ik), intent(in) :: depth
   logical(lk), dimension(:), intent(inout) :: is_last_child
   type(chaintree_node_t), pointer :: child_chain, next_child
   integer(ik) :: i

   if (.not. associated(chain)) return

   child_chain => chain%first_child_chain
   do while (associated(child_chain))
      ! Tree drawing prefix
      next_child => child_chain%next_sibling_chain
      is_last_child(depth + 1) = .not. associated(next_child)
      do i = 1, depth
         if (is_last_child(i)) then
            write(stderr, '(A)', advance='no') "    "
         else
            write(stderr, '(A)', advance='no') "|   "
         end if
      end do

      if (is_last_child(depth + 1)) then
         write(stderr, '(A)', advance='no') "`--"
      else
         write(stderr, '(A)', advance='no') "|--"
      end if

      ! Split part address with item counts
      write(stderr, '(A,1X,A,I0,A,I0,A)') address(child_chain%split_part), &
         '(', child_chain%split_part%n_items1, '/', child_chain%split_part%n_items2, ')'

      call print_chain_recurse(child_chain, depth + 1, is_last_child)

      child_chain => next_child
   end do
end subroutine

function is_partition_uneven(link) result(uneven)
! Whether some part of link holds different numbers of atoms of each molecule
   type(chain_node_t), pointer, intent(in) :: link
   logical(lk) :: uneven
   type(partref_node_t), pointer :: partref

   uneven = .FALSE.
   partref => link%first_partref
   do while (associated(partref))
      if (partref%part%n_items1 /= partref%part%n_items2) then
         uneven = .TRUE.
         return
      end if
      partref => partref%nextref
   end do
end function

subroutine print_partition_details(link, link_number)
! Print the parts of a link, even parts (as many atoms of each molecule)
! first, then uneven parts
   type(chain_node_t), pointer, intent(in) :: link
   integer(ik), intent(in) :: link_number
   type(partref_node_t), pointer :: partref
   integer(ik) :: n_even, n_uneven

   write(stderr, '(A)') repeat("=", 70)
   write(stderr, '(A,I0,A,I0,A)') "PARTITION LINK ", link_number, " (", link%n_parts, " parts)"
   write(stderr, '(A)') repeat("=", 70)

   n_even = 0
   n_uneven = 0
   partref => link%first_partref
   do while (associated(partref))
      if (partref%part%n_items1 == partref%part%n_items2) then
         n_even = n_even + 1
      else
         n_uneven = n_uneven + 1
      end if
      partref => partref%nextref
   end do

   write(stderr, '(A)') "Even parts:"
   if (n_even == 0) then
      write(stderr, '(A)') "  (none)"
   else
      partref => link%first_partref
      do while (associated(partref))
         if (partref%part%n_items1 == partref%part%n_items2) then
            call print_part_line(partref%part)
         end if
         partref => partref%nextref
      end do
   end if

   if (n_uneven > 0) then
      write(stderr, *)
      write(stderr, '(A)') "Uneven parts:"
      partref => link%first_partref
      do while (associated(partref))
         if (partref%part%n_items1 /= partref%part%n_items2) then
            call print_part_line(partref%part)
         end if
         partref => partref%nextref
      end do
   end if

   write(stderr, '(A)') repeat("=", 70)
   write(stderr, *)
end subroutine

subroutine print_part_line(part)
! Print a part and its atoms on one line
   type(partition_node_t), pointer, intent(in) :: part
   type(item_node_t), pointer :: item

   write(stderr, '(A,A,I0,A,I0,A)', advance='no') &
      address(part), " (", part%n_items1, "/", part%n_items2, "): ["

   item => part%first_item1
   do while (associated(item))
      write(stderr, '(I0)', advance='no') item%idx
      item => item%next_item
      if (associated(item)) write(stderr, '(A)', advance='no') ","
   end do

   write(stderr, '(A)', advance='no') "] / ["

   item => part%first_item2
   do while (associated(item))
      write(stderr, '(I0)', advance='no') item%idx
      item => item%next_item
      if (associated(item)) write(stderr, '(A)', advance='no') ","
   end do

   write(stderr, '(A)') "]"
end subroutine

subroutine print_partition_chain(hna_chain)
! Print all links (refinement levels) of a chain
   type(chaintree_node_t), pointer, intent(in) :: hna_chain
   type(chain_node_t), pointer :: link
   integer(ik) :: link_number

   write(stderr, *)
   write(stderr, '(A)') repeat("#", 70)
   write(stderr, '(A)') "                    PARTITION CHAIN HISTORY"
   write(stderr, '(A)') repeat("#", 70)
   write(stderr, *)

   link => hna_chain%first_link
   link_number = 1

   do while (associated(link))
      call print_partition_details(link, link_number)
      link => link%next_link
      link_number = link_number + 1
   end do
end subroutine

end module
