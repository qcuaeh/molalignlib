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

module types_linked
use parameters
use types_basic
implicit none
private
! Public procedures
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

! Item node
type, public :: item_node_t
   integer :: idx
   integer :: global_idx  ! Global index
   type(item_node_t), pointer :: next_item
end type

! Chain node
type, public :: chain_node_t
   integer :: num_parts
   integer :: global_idx
   integer, pointer :: total_partrefs => null()
   type(chain_node_t), pointer :: next_link
   type(partref_node_t), pointer :: first_partref
   type(partref_node_t), pointer :: last_partref
   type(part_nodeptr_t), dimension(:), pointer :: itemdir1
   type(part_nodeptr_t), dimension(:), pointer :: itemdir2
end type

! Assignment tree node
type, public :: chaintree_node_t
   integer :: num_links
   integer :: num_children
   integer :: global_idx
   integer, pointer :: n_atoms1 => null()
   integer, pointer :: n_atoms2 => null()
   integer, pointer :: total_chains => null()
   integer, pointer :: total_links => null()
   integer, pointer :: total_partrefs => null()
   type(chaintree_node_t), pointer :: parent_chain
   type(chaintree_node_t), pointer :: first_child_chain
   type(chaintree_node_t), pointer :: last_child_chain
   type(chaintree_node_t), pointer :: next_sibling_chain
   type(chain_node_t), pointer :: first_link
   type(chain_node_t), pointer :: last_link
   type(partition_node_t), pointer :: split_part
end type

! Part tree node
type, public :: partition_node_t
   integer :: depth
   integer :: global_idx
   integer :: num_items1
   integer :: num_items2
   integer :: num_children
   integer, pointer :: total_parts => null()
   integer, pointer :: total_items1 => null()
   integer, pointer :: total_items2 => null()
   type(item_node_t), pointer :: first_item1
   type(item_node_t), pointer :: first_item2
   type(item_node_t), pointer :: last_item1
   type(item_node_t), pointer :: last_item2
   type(partition_node_t), pointer :: parent_part
   type(partition_node_t), pointer :: first_child_part
   type(partition_node_t), pointer :: last_child_part
   type(partition_node_t), pointer :: next_sibling_part
   type(part_nodeptr_t), dimension(:), pointer :: signature
end type

! Part reference node
type, public :: partref_node_t
   integer :: global_idx
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
   logical :: equality
   equality = associated(left%ptr, right%ptr)
end function

function signature_equivalence(signature1, signature2) result(equiv)
   type(part_nodeptr_t), dimension(:), intent(in) :: signature1, signature2
   logical :: equiv
   integer :: matches1, matches2
   integer :: i, j

   if (size(signature1) /= size(signature2)) then
      equiv = .FALSE.
      return
   end if

   do i = 1, size(signature1)
      matches1 = 0
      matches2 = 0
      do j = 1, size(signature1)
         if (associated(signature1(i)%ptr, signature2(j)%ptr)) then
            matches1 = matches1 + 1
         end if
         if (associated(signature1(i)%ptr, signature1(j)%ptr)) then
            matches2 = matches2 + 1
         end if
      end do
      if (matches1 /= matches2) then
         equiv = .FALSE.
         return
      end if
   end do

   equiv = .TRUE.
end function

logical function isdescendant(part, top_part)
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

   part%global_idx = 0  ! Will be set when added to tree
   part%num_items1 = 0
   part%num_items2 = 0
   part%num_children = 0
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

   ! Allocate counters for root and assign initial values
   allocate(part%total_parts)
   allocate(part%total_items1)
   allocate(part%total_items2)

   part%total_parts = 1
   part%total_items1 = 0  ! No items initially
   part%total_items2 = 0  ! No items initially
end function

function new_child_part(parent_part) result(child_part)
   type(partition_node_t), pointer, intent(inout) :: parent_part
   type(partition_node_t), pointer :: child_part

   ! Create new part with correct depth
   child_part => new_bare_part()
   child_part%depth = parent_part%depth + 1

   ! Set parent BEFORE setting up relationships
   child_part%parent_part => parent_part
   child_part%next_sibling_part => null()

   ! Point to parent's counters
   child_part%total_parts => parent_part%total_parts
   child_part%total_items1 => parent_part%total_items1
   child_part%total_items2 => parent_part%total_items2

   ! Assign global index and increment counter
   child_part%total_parts = child_part%total_parts + 1
   child_part%global_idx = child_part%total_parts

   ! Set up parent-child relationships (NO link association)
   if (.not. associated(parent_part%first_child_part)) then
      parent_part%first_child_part => child_part
   else
      parent_part%last_child_part%next_sibling_part => child_part
   end if
   parent_part%last_child_part => child_part
   parent_part%num_children = parent_part%num_children + 1
end function

subroutine link_part(link, part)
   type(chain_node_t), target, intent(inout) :: link
   type(partition_node_t), target, intent(inout) :: part
   type(partref_node_t), pointer :: newref

   allocate(newref)
   newref%part => part
   newref%nextref => null()

   ! Assign global index and increment counter
   link%total_partrefs = link%total_partrefs + 1
   newref%global_idx = link%total_partrefs

   if (.not. associated(link%first_partref)) then
      link%first_partref => newref
   else
      link%last_partref%nextref => newref
   end if
   link%last_partref => newref
   link%num_parts = link%num_parts + 1
end subroutine

subroutine add_new_item1(part, idx)
! Add item to part without updating link itemdir (for temporary children)
   type(partition_node_t), target, intent(inout) :: part
   integer, intent(in) :: idx
   type(item_node_t), pointer :: new_item

   allocate(new_item)
   new_item%idx = idx
   new_item%next_item => null()

   ! Assign global index and increment counter
   part%total_items1 = part%total_items1 + 1
   new_item%global_idx = part%total_items1

   if (.not. associated(part%first_item1)) then
      part%first_item1 => new_item
   else
      part%last_item1%next_item => new_item
   end if
   part%last_item1 => new_item
   part%num_items1 = part%num_items1 + 1
end subroutine

subroutine add_new_item2(part, idx)
! Add item to part without updating link itemdir (for temporary children)
   type(partition_node_t), target, intent(inout) :: part
   integer, intent(in) :: idx
   type(item_node_t), pointer :: new_item

   allocate(new_item)
   new_item%idx = idx
   new_item%next_item => null()

   ! Assign global index and increment counter
   part%total_items2 = part%total_items2 + 1
   new_item%global_idx = part%total_items2

   if (.not. associated(part%first_item2)) then
      part%first_item2 => new_item
   else
      part%last_item2%next_item => new_item
   end if
   part%last_item2 => new_item
   part%num_items2 = part%num_items2 + 1
end subroutine

subroutine copy_part_items(orig, dest)
! Move all items from orig to dest without updating link itemdir
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
! Move first item from orig to dest without updating link itemdir
   type(partition_node_t), intent(inout) :: orig
   type(partition_node_t), target, intent(inout) :: dest
   type(item_node_t), pointer :: item_to_move

   item_to_move => orig%first_item1
   if ((.not. associated(item_to_move))) error stop

   ! Remove from source
   orig%first_item1 => item_to_move%next_item
   if (.not. associated(orig%first_item1)) then
      orig%last_item1 => null()
   end if
   orig%num_items1 = orig%num_items1 - 1

   ! Add to destination
   item_to_move%next_item => null()

   if (.not. associated(dest%first_item1)) then
      dest%first_item1 => item_to_move
   else
      dest%last_item1%next_item => item_to_move
   end if
   dest%last_item1 => item_to_move
   dest%num_items1 = dest%num_items1 + 1
end subroutine

subroutine move_first_item2(orig, dest)
! Move first item from orig to dest without updating link itemdir
   type(partition_node_t), intent(inout) :: orig
   type(partition_node_t), target, intent(inout) :: dest
   type(item_node_t), pointer :: item_to_move

   item_to_move => orig%first_item2
   if ((.not. associated(item_to_move))) error stop

   ! Remove from source
   orig%first_item2 => item_to_move%next_item
   if (.not. associated(orig%first_item2)) then
      orig%last_item2 => null()
   end if
   orig%num_items2 = orig%num_items2 - 1

   ! Add to destination
   item_to_move%next_item => null()

   if (.not. associated(dest%first_item2)) then
      dest%first_item2 => item_to_move
   else
      dest%last_item2%next_item => item_to_move
   end if
   dest%last_item2 => item_to_move
   dest%num_items2 = dest%num_items2 + 1
end subroutine

subroutine move_part_items(orig, dest)
! Move all items from orig to dest without updating link itemdir
   type(partition_node_t), intent(inout) :: orig, dest

   do while (associated(orig%first_item1))
      call move_first_item1(orig, dest)
   end do

   do while (associated(orig%first_item2))
      call move_first_item2(orig, dest)
   end do
end subroutine

subroutine delete_chain(root_chain)
   type(chaintree_node_t), pointer, intent(inout) :: root_chain
   type(chain_node_t), pointer :: link, next_link

   if ((.not. associated(root_chain))) error stop

   ! Delete all links (partitions only - parts are preserved)
   link => root_chain%first_link
   do while (associated(link))
      next_link => link%next_link
      call delete_link(link)
      link => next_link
   end do

   ! Then deallocate the chain root
   deallocate(root_chain)
   root_chain => null()
end subroutine

subroutine delete_link(link)
   type(chain_node_t), pointer, intent(inout) :: link
   type(partref_node_t), pointer :: partref, nextref

   if ((.not. associated(link))) error stop

   ! Delete all part references WITHOUT deleting the parts themselves
   partref => link%first_partref
   do while (associated(partref))
      nextref => partref%nextref

      ! Only delete the part reference, NOT the part itself
      deallocate(partref)
      partref => nextref
   end do

   ! Deallocate directories (always allocated)
   deallocate(link%itemdir1)
   deallocate(link%itemdir2)

   ! Deallocate partition root
   deallocate(link)
   link => null()
end subroutine

subroutine delete_part_tree(partition_tree)
! Deletes an entire part tree starting from the root
   type(partition_node_t), pointer, intent(inout) :: partition_tree

   if (.not. associated(partition_tree)) error stop

   ! Recursively delete all children first
   call delete_part_children(partition_tree)

   ! Deallocate counters
   deallocate(partition_tree%total_parts)
   deallocate(partition_tree%total_items1)
   deallocate(partition_tree%total_items2)

   ! Then delete the root part itself
   call delete_part(partition_tree)
end subroutine

recursive subroutine delete_part_children(parent_part)
! Deletes all children of a part recursively
   type(partition_node_t), pointer, intent(in) :: parent_part
   type(partition_node_t), pointer :: child_part, next_child

   if (.not. associated(parent_part)) error stop

   ! Traverse and delete all children
   child_part => parent_part%first_child_part
   do while (associated(child_part))
      next_child => child_part%next_sibling_part

      ! Recursively delete this child's subtree
      call delete_part_children(child_part)
      call delete_part(child_part)

      child_part => next_child
   end do
end subroutine

subroutine delete_part(part_node)
! Deletes a single part node (assumes children are already deleted)
   type(partition_node_t), pointer, intent(inout) :: part_node

   if ((.not. associated(part_node))) error stop

   ! Deallocate items
   call deallocate_items(part_node%first_item1)
   call deallocate_items(part_node%first_item2)

   ! Deallocate signature array if allocated
   if (associated(part_node%signature)) then
      deallocate(part_node%signature)
   end if

   ! Deallocate the part itself
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
   type(chain_node_t), target, intent(in) :: link
   type(partition_t), intent(out) :: partition
   type(partref_node_t), pointer :: partref
   type(item_node_t), pointer :: item
   integer :: i, j

   partition%num_parts = link%num_parts
   allocate(partition%parts(partition%num_parts))
   allocate(partition%itemdir1(size(link%itemdir1)))
   allocate(partition%itemdir2(size(link%itemdir2)))

   partref => link%first_partref
   i = 1
   do while (associated(partref))
      if (.not. associated(partref%part)) error stop

      partition%parts(i)%num_items1 = partref%part%num_items1
      partition%parts(i)%num_items2 = partref%part%num_items2
      allocate(partition%parts(i)%items1(partref%part%num_items1))
      allocate(partition%parts(i)%items2(partref%part%num_items2))

      item => partref%part%first_item1
      do j = 1, partref%part%num_items1
         partition%parts(i)%items1(j) = item%idx
         partition%itemdir1(item%idx) = i
         item => item%next_item
      end do

      item => partref%part%first_item2
      do j = 1, partref%part%num_items2
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
   integer :: i

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

   ! Counters will be set when added to tree
   chain%num_links = 0
   chain%num_children = 0
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
   integer, intent(in) :: n_atoms1, n_atoms2
   type(chaintree_node_t), pointer :: chain

   chain => new_bare_chain()
   chain%parent_chain => null()
   chain%split_part => null()
   chain%global_idx = 1

   ! Allocate counters for root and assign initial values
   allocate(chain%n_atoms1)
   allocate(chain%n_atoms2)
   allocate(chain%total_chains)
   allocate(chain%total_links)
   allocate(chain%total_partrefs)
   chain%n_atoms1 = n_atoms1
   chain%n_atoms2 = n_atoms2
   chain%total_chains = 1
   chain%total_links = 0    ! No links initially
   chain%total_partrefs = 0 ! No partrefs initially
end function

function new_bare_link() result(link)
   type(chain_node_t), pointer :: link

   allocate(link)
   link%num_parts = 0
   link%global_idx = 0  ! Will be set when added to chain
   link%total_partrefs => null()
   link%first_partref => null()
   link%last_partref => null()
   link%next_link => null()
end function

function new_chain_link(chain) result(link)
   type(chaintree_node_t), target, intent(inout) :: chain
   type(chain_node_t), pointer :: link
   integer :: i

   ! Use the initialization function
   link => new_bare_link()

   ! Point to chain's partref counter
   link%total_partrefs => chain%total_partrefs

   ! Assign global index and increment counter
   chain%total_links = chain%total_links + 1
   link%global_idx = chain%total_links

   ! Allocate item directories
   allocate(link%itemdir1(chain%n_atoms1))
   allocate(link%itemdir2(chain%n_atoms2))

   ! Initialize all pointers to null
   do i = 1, chain%n_atoms1
      link%itemdir1(i)%ptr => null()
   end do
   do i = 1, chain%n_atoms2
      link%itemdir2(i)%ptr => null()
   end do

   ! Add link to chain
   if (.not. associated(chain%first_link)) then
      chain%first_link => link
   else
      chain%last_link%next_link => link
   end if
   chain%last_link => link
   chain%num_links = chain%num_links + 1
end function

function find_child_part(part, signature) result(child_part)
   type(partition_node_t), intent(in) :: part
   type(part_nodeptr_t), dimension(:), intent(in) :: signature
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
   logical, dimension(:), allocatable :: is_last_child

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

   ! Print root line
   write(stderr, '(A)') 'ROOT'

   ! Print children recursively
   call print_part_recurse(partition_tree, 0, is_last_child)
   write(stderr, *)

   deallocate(is_last_child)
end subroutine

function chain_from_partition(partition) result(chain)
   type(partition_t), intent(in) :: partition
   ! Local variables
   type(chaintree_node_t), pointer :: chain
   type(partition_node_t), pointer :: partition_tree
   type(chain_node_t), pointer :: first_link
   type(partition_node_t), pointer :: new_part
   integer :: i, j

   ! Create root part (decoupled from chain)
   partition_tree => new_root_part()

   ! Create new chain
   chain => new_root_chain(size(partition%itemdir1), size(partition%itemdir2))

   ! Create first link
   first_link => new_chain_link(chain)

   ! Create parts as children of partition_tree and add them to the first link
   do i = 1, partition%num_parts
      ! Create new part as child of partition_tree
      new_part => new_child_part(partition_tree)

      ! Add items to the child part using cached approach
      do j = 1, partition%parts(i)%num_items1
         call add_new_item1(new_part, partition%parts(i)%items1(j))
      end do

      do j = 1, partition%parts(i)%num_items2
         call add_new_item2(new_part, partition%parts(i)%items2(j))
      end do

      ! Link child part to first link
      call link_part(first_link, new_part)
      call update_itemdir(first_link, new_part)
   end do
end function

function new_child_chain(chain, split_part) result(new_chain)
   type(chaintree_node_t), pointer, intent(inout) :: chain
   type(partition_node_t), pointer, intent(in) :: split_part
   type(chaintree_node_t), pointer :: new_chain

   ! Create new child chain
   new_chain => new_bare_chain()
   new_chain%split_part => split_part

   ! Set up parent-child relationship
   new_chain%parent_chain => chain

   ! Point to parent's counters
   new_chain%n_atoms1 => chain%n_atoms1
   new_chain%n_atoms2 => chain%n_atoms2
   new_chain%total_chains => chain%total_chains
   new_chain%total_links => chain%total_links
   new_chain%total_partrefs => chain%total_partrefs

   ! Assign global index and increment counter
   new_chain%total_chains = new_chain%total_chains + 1
   new_chain%global_idx = new_chain%total_chains

   ! Add as child to parent (optimized with last_child_chain pointer)
   if (.not. associated(chain%first_child_chain)) then
      ! This is the first child
      chain%first_child_chain => new_chain
      chain%last_child_chain => new_chain
   else
      ! Add as next sibling to the current last child
      chain%last_child_chain%next_sibling_chain => new_chain
      chain%last_child_chain => new_chain
   end if

   chain%num_children = chain%num_children + 1
end function

subroutine add_branch_part(link, part)
! Adds parts in sorted order by num_items1
   type(chain_node_t), target, intent(inout) :: link
   type(partition_node_t), target, intent(in) :: part
   type(partref_node_t), pointer :: partref, prevref, newref

   ! Do not add assigned parts
   if (part%num_items1 < 2) return

   ! Create new part reference
   allocate(newref)
   newref%part => part
   newref%nextref => null()

   ! If link is empty, add as first element
   if (.not. associated(link%first_partref)) then
      link%first_partref => newref
      link%last_partref => newref
      link%num_parts = link%num_parts + 1
      return
   end if

   ! Find correct insertion position (sorted by increasing num_items1)
   partref => link%first_partref
   prevref => null()

   do while (associated(partref))
      ! If new part has fewer or equal items1, insert before partref
      if (part%num_items1 <= partref%part%num_items1) then
         exit
      end if
      prevref => partref
      partref => partref%nextref
   end do

   ! Insert newref at the found position
   newref%nextref => partref

   if (associated(prevref)) then
      ! Insert in middle or at end
      prevref%nextref => newref
      ! Update last_partref if we inserted at the end
      if (.not. associated(partref)) then
         link%last_partref => newref
      end if
   else
      ! Insert at beginning
      link%first_partref => newref
      ! If this was the only element, also update last_partref
      if (.not. associated(partref)) then
         link%last_partref => newref
      end if
   end if

   link%num_parts = link%num_parts + 1
end subroutine

subroutine update_itemdir(link, part)
! Updates itemdir pointers for all items in a single part
   type(chain_node_t), target, intent(inout) :: link
   type(partition_node_t), target, intent(in) :: part
   type(item_node_t), pointer :: item

   ! Register all items from molecule 1
   item => part%first_item1
   do while (associated(item))
      link%itemdir1(item%idx)%ptr => part
      item => item%next_item
   end do

   ! Register all items from molecule 2
   item => part%first_item2
   do while (associated(item))
      link%itemdir2(item%idx)%ptr => part
      item => item%next_item
   end do
end subroutine

subroutine print_chain_tree(assignment_tree)
   type(chaintree_node_t), pointer, intent(in) :: assignment_tree
   logical, dimension(:), allocatable :: is_last_child

   if (.not. associated(assignment_tree)) error stop

   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "     Assignment Tree"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Allocate tracking array for tree lines (max depth 100)
   allocate(is_last_child(100))
   is_last_child = .FALSE.

   ! Print children recursively
   write(stderr, '(A)') 'ROOT'
   call print_chain_recurse(assignment_tree, 0, is_last_child)
   write(stderr, *)

   deallocate(is_last_child)
end subroutine

subroutine print_part_signature(signature)
   type(part_nodeptr_t), dimension(:), intent(in) :: signature
   integer :: i

   do i = 1, size(signature)
      write(stderr,'(*(1X,A))',advance='no') address(signature(i)%ptr)
   end do
   write(stderr,*)
end subroutine

! New procedure to print signatures of all parts in a part tree
subroutine print_tree_signatures(partition_tree)
   type(partition_node_t), pointer, intent(in) :: partition_tree

   if (.not. associated(partition_tree)) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "    Part Signatures"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Print children recursively
   call print_signature_recurse(partition_tree)
   write(stderr, *)
end subroutine

! Helper recursive procedure for print_tree_signatures
recursive subroutine print_signature_recurse(part)
   type(partition_node_t), pointer, intent(in) :: part
   type(partition_node_t), pointer :: child_part

   if (.not. associated(part)) return

   ! Process all children in the same order as print_part_tree
   child_part => part%first_child_part
   do while (associated(child_part))
      ! Print the child signature
      write(stderr,'(A)',advance='no') address(child_part) // ':'
      call print_part_signature(child_part%signature)

      ! Recursively print this child's children
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

   ! Print children recursively
   call print_items_recurse(partition_tree)
   write(stderr, *)
end subroutine

recursive subroutine print_items_recurse(part)
   type(partition_node_t), pointer, intent(in) :: part
   type(partition_node_t), pointer :: child_part

   if (.not. associated(part)) return

   ! Process all children in the same order as print_part_tree
   child_part => part%first_child_part
   do while (associated(child_part))
      ! Print the child items with address prefix
      write(stderr, '(A)', advance='no') address(child_part) // ':'
      call print_part_items(child_part)

      ! Recursively print this child's children
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

   ! Print children recursively
   call print_leaf_items_recurse(partition_tree)
   write(stderr, *)
end subroutine

recursive subroutine print_leaf_items_recurse(part)
   type(partition_node_t), pointer, intent(in) :: part
   type(partition_node_t), pointer :: child_part

   if (.not. associated(part)) return

   ! Process all children in the same order as print_part_tree
   child_part => part%first_child_part
   do while (associated(child_part))
      ! Only print items if this is a leaf part (no children)
      if (child_part%num_children == 0) then
         write(stderr, '(A)', advance='no') address(child_part) // ':'
         call print_part_items(child_part)
      end if

      ! Recursively traverse this child's children to find more leaves
      call print_leaf_items_recurse(child_part)

      child_part => child_part%next_sibling_part
   end do
end subroutine

! New procedure to print global indices
subroutine print_part_indices(partition_tree)
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

   ! Print children recursively
   call print_part_indices_recurse(partition_tree)
   write(stderr, *)
end subroutine

! Helper recursive procedure for print_part_indices
recursive subroutine print_part_indices_recurse(part)
   type(partition_node_t), pointer, intent(in) :: part
   type(partition_node_t), pointer :: child_part
   type(item_node_t), pointer :: item

   if (.not. associated(part)) return

   ! Process all children
   child_part => part%first_child_part
   do while (associated(child_part))
      ! Print the child part index
      write(stderr, '(A,I0,A,A,A)') "Part ", child_part%global_idx, " (address: ", address(child_part), ")"

      ! Print items in this part
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

      ! Recursively print this child's children
      call print_part_indices_recurse(child_part)

      child_part => child_part%next_sibling_part
   end do
end subroutine

! New procedure to print global indices for chain tree
subroutine print_chain_indices(assignment_tree)
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

   ! Print children recursively
   call print_chain_indices_recurse(assignment_tree)
   write(stderr, *)
end subroutine

! Helper recursive procedure for print_chain_indices
recursive subroutine print_chain_indices_recurse(chain)
   type(chaintree_node_t), pointer, intent(in) :: chain
   type(chaintree_node_t), pointer :: child_chain
   type(chain_node_t), pointer :: link
   type(partref_node_t), pointer :: partref

   if (.not. associated(chain)) return

   ! Print links in this chain
   link => chain%first_link
   do while (associated(link))
      write(stderr, '(A,I0,A,I0,A)') "  Link ", link%global_idx, " (", link%num_parts, " parts)"

      ! Print partrefs in this link
      partref => link%first_partref
      do while (associated(partref))
         write(stderr, '(A,I0,A,A,A)') "    Partref ", partref%global_idx, " (part: ", address(partref%part), ")"
         partref => partref%nextref
      end do

      link => link%next_link
   end do

   ! Process all child chains
   child_chain => chain%first_child_chain
   do while (associated(child_chain))
      ! Print the child chain index
      write(stderr, '(A,I0,A,A,A)') "Chain ", child_chain%global_idx, " (split part: ", address(child_chain%split_part), ")"

      ! Recursively print this child's content and children
      call print_chain_indices_recurse(child_chain)

      child_chain => child_chain%next_sibling_chain
   end do
end subroutine

recursive subroutine print_part_recurse(part, depth, is_last_child)
   type(partition_node_t), pointer, intent(in) :: part
   integer, intent(in) :: depth
   logical, dimension(:), intent(inout) :: is_last_child
   type(partition_node_t), pointer :: child_part, next_child
   integer :: i

   if (.not. associated(part)) return

   ! Process all children
   child_part => part%first_child_part
   do while (associated(child_part))
      ! Check if this is the last child
      next_child => child_part%next_sibling_part
      is_last_child(depth + 1) = .not. associated(next_child)

      ! Print prefix components directly
      do i = 1, depth
         if (is_last_child(i)) then
            write(stderr, '(A)', advance='no') "    "
         else
            write(stderr, '(A)', advance='no') "|   "
         end if
      end do

      ! Add branch characters
      if (is_last_child(depth + 1)) then
         write(stderr, '(A)', advance='no') "`--"
      else
         write(stderr, '(A)', advance='no') "|--"
      end if

      ! Print part address with item counts
      write(stderr, '(A,1X,A,I0,A,I0,A)') address(child_part), &
         '(', child_part%num_items1, '/', child_part%num_items2, ')'

      ! Recursively print this child's children
      call print_part_recurse(child_part, depth + 1, is_last_child)

      child_part => next_child
   end do
end subroutine

recursive subroutine print_chain_recurse(chain, depth, is_last_child)
   type(chaintree_node_t), pointer, intent(in) :: chain
   integer, intent(in) :: depth
   logical, dimension(:), intent(inout) :: is_last_child
   type(chaintree_node_t), pointer :: child_chain, next_child
   integer :: i

   if (.not. associated(chain)) return

   ! Process all children
   child_chain => chain%first_child_chain
   do while (associated(child_chain))
      ! Check if this is the last child
      next_child => child_chain%next_sibling_chain
      is_last_child(depth + 1) = .not. associated(next_child)

      ! Print prefix components directly
      do i = 1, depth
         if (is_last_child(i)) then
            write(stderr, '(A)', advance='no') "    "
         else
            write(stderr, '(A)', advance='no') "|   "
         end if
      end do

      ! Add branch characters
      if (is_last_child(depth + 1)) then
         write(stderr, '(A)', advance='no') "`--"
      else
         write(stderr, '(A)', advance='no') "|--"
      end if

      ! Print the child address with item counts
      write(stderr, '(A,1X,A,I0,A,I0,A)') address(child_chain%split_part), &
         '(', child_chain%split_part%num_items1, '/', child_chain%split_part%num_items2, ')'

      ! Recursively print this child's children
      call print_chain_recurse(child_chain, depth + 1, is_last_child)

      child_chain => next_child
   end do
end subroutine

function is_partition_uneven(link) result(uneven)
! Check if all parts in a partition have equal numbers of items from both molecules
   type(chain_node_t), pointer, intent(in) :: link
   logical :: uneven
   type(partref_node_t), pointer :: partref

   uneven = .FALSE.
   partref => link%first_partref
   do while (associated(partref))
      if (partref%part%num_items1 /= partref%part%num_items2) then
         uneven = .TRUE.
         return
      end if
      partref => partref%nextref
   end do
end function

subroutine print_partition_details(link, link_number)
! Print detailed information about a partition link in compact form
! Separates even and uneven parts
   type(chain_node_t), pointer, intent(in) :: link
   integer, intent(in) :: link_number
   type(partref_node_t), pointer :: partref
   integer :: even_count, uneven_count

   write(stderr, '(A)') repeat("=", 70)
   write(stderr, '(A,I0,A,I0,A)') "PARTITION LINK ", link_number, " (", link%num_parts, " parts)"
   write(stderr, '(A)') repeat("=", 70)

   ! First pass: count even and uneven parts
   even_count = 0
   uneven_count = 0
   partref => link%first_partref
   do while (associated(partref))
      if (partref%part%num_items1 == partref%part%num_items2) then
         even_count = even_count + 1
      else
         uneven_count = uneven_count + 1
      end if
      partref => partref%nextref
   end do

   ! Print even parts
   write(stderr, '(A)') "Even parts:"
   if (even_count == 0) then
      write(stderr, '(A)') "  (none)"
   else
      partref => link%first_partref
      do while (associated(partref))
         if (partref%part%num_items1 == partref%part%num_items2) then
            call print_part_line(partref%part)
         end if
         partref => partref%nextref
      end do
   end if

   ! Print uneven parts only if they exist
   if (uneven_count > 0) then
      write(stderr, *)
      write(stderr, '(A)') "Uneven parts:"
      partref => link%first_partref
      do while (associated(partref))
         if (partref%part%num_items1 /= partref%part%num_items2) then
            call print_part_line(partref%part)
         end if
         partref => partref%nextref
      end do
   end if

   write(stderr, '(A)') repeat("=", 70)
   write(stderr, *)
end subroutine

subroutine print_part_line(part)
! Print a single part in compact form
   type(partition_node_t), pointer, intent(in) :: part
   type(item_node_t), pointer :: item

   ! Print part address and counts
   write(stderr, '(A,A,I0,A,I0,A)', advance='no') &
      address(part), " (", part%num_items1, "/", part%num_items2, "): ["

   ! Print molecule 1 atoms
   item => part%first_item1
   do while (associated(item))
      write(stderr, '(I0)', advance='no') item%idx
      item => item%next_item
      if (associated(item)) write(stderr, '(A)', advance='no') ","
   end do

   write(stderr, '(A)', advance='no') "] / ["

   ! Print molecule 2 atoms
   item => part%first_item2
   do while (associated(item))
      write(stderr, '(I0)', advance='no') item%idx
      item => item%next_item
      if (associated(item)) write(stderr, '(A)', advance='no') ","
   end do

   write(stderr, '(A)') "]"
end subroutine

subroutine print_partition_chain(hna_chain)
! Print all links in the chain history
   type(chaintree_node_t), pointer, intent(in) :: hna_chain
   type(chain_node_t), pointer :: link
   integer :: link_number

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
