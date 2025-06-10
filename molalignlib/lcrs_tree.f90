module lcrs_tree
use parameters
implicit none
private

! Item node
type, public :: item_node_t
   integer :: value
   integer :: global_index
   type(item_node_t), pointer :: next_item
end type

! Part reference node
type, public :: partref_node_t
   integer :: global_index
   type(part_node_t), pointer :: part
   type(partref_node_t), pointer :: nextref
end type

! Link node
type, public :: link_node_t
   integer :: num_parts
   integer :: global_index
   integer, pointer :: total_partrefs => null()
   type(link_node_t), pointer :: next_link
   type(partref_node_t), pointer :: first_partref
   type(partref_node_t), pointer :: last_partref
   type(part_nodeptr_t), dimension(:), allocatable :: itemdir1
   type(part_nodeptr_t), dimension(:), allocatable :: itemdir2
end type

! LCRS chain node
type, public :: chain_node_t
   integer :: tot_items1
   integer :: tot_items2
   integer :: num_links
   integer :: num_children
   integer :: global_index
   integer, pointer :: total_chains => null()
   integer, pointer :: total_links => null()
   integer, pointer :: total_partrefs => null()
   type(chain_node_t), pointer :: parent_chain
   type(chain_node_t), pointer :: first_child_chain
   type(chain_node_t), pointer :: last_child_chain
   type(chain_node_t), pointer :: next_sibling_chain
   type(link_node_t), pointer :: first_link
   type(link_node_t), pointer :: last_link
   type(part_node_t), pointer :: split_part
end type

! LCRS part node
type, public :: part_node_t
   integer :: depth
   integer :: global_index
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
   type(part_node_t), pointer :: parent_part
   type(part_node_t), pointer :: first_child_part
   type(part_node_t), pointer :: last_child_part
   type(part_node_t), pointer :: next_sibling_part
   type(part_nodeptr_t), dimension(:), allocatable :: signature
end type

! Part node pointer
type, public :: part_nodeptr_t
   type(part_node_t), pointer :: ptr
end type

! Semipart array
type, public :: semipartition_part_t
   integer :: num_items
   integer, dimension(:), allocatable :: items
end type

! Semipartition array
type, public :: semipartition_t
   integer :: num_parts
   type(semipartition_part_t), dimension(:), allocatable :: parts
   integer, dimension(:), allocatable :: itemdir
end type

! Part array
type, public :: partition_part_t
   integer :: num_items1
   integer :: num_items2
   integer :: num_children
   integer, dimension(:), allocatable :: items1
   integer, dimension(:), allocatable :: items2
   integer, dimension(:), allocatable :: signature
   integer, dimension(:), allocatable :: children
end type

! Partition array
type, public :: partition_t
   integer :: num_parts
   type(partition_part_t), dimension(:), allocatable :: parts
   integer, dimension(:), allocatable :: itemdir1
   integer, dimension(:), allocatable :: itemdir2
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

! Make types and procedures public
public address
public isdescendant
public new_root_chain
public new_child_chain
public new_bare_link
public new_chain_link
public new_child_part
public link_part
public update_itemdir
public add_new_item1
public add_new_item2
public init_chain_from_link
public init_chain_from_partition
public move_first_item1
public move_first_item2
public move_part_items
public delete_chain
public print_part_tree
public print_part_items
public print_tree_items
public print_leaf_items
public print_link_itemdir
public print_part_signature
public print_tree_signatures
public print_part_indices
public print_chain_indices
public link_to_partition
public find_child_part
public first_partition
public second_partition
public add_branch_part
public print_chain_tree
public delete_part_tree
public delete_part
public operator(==)
public operator(.equiv.)


contains

character(4) function address_part(nodeptr) result(address)
   use iso_c_binding, only: c_loc, c_intptr_t
   type(part_node_t), pointer, intent(in) :: nodeptr
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

function signature_equivalence(array1, array2) result(equiv)
   type(part_nodeptr_t), dimension(:), intent(in) :: array1, array2
   logical :: equiv
   integer :: matches1, matches2
   integer :: i, j

   if (size(array1) /= size(array2)) then
      equiv = .false.
      return
   end if

   do i = 1, size(array1)
       matches1 = 0
       matches2 = 0
       do j = 1, size(array1)
           if (associated(array1(i)%ptr, array2(j)%ptr)) matches1 = matches1 + 1
           if (associated(array1(i)%ptr, array1(j)%ptr)) matches2 = matches2 + 1
       end do
       if (matches1 /= matches2) then
           equiv = .false.
           return
       end if
   end do

   equiv = .true.
end function

logical function isdescendant(part, top_part)
   type(part_node_t), pointer, intent(in) :: part, top_part
   ! Local variables
   type(part_node_t), pointer :: up_part

   if (part%depth <= top_part%depth) then
      isdescendant = .false.
      return
   end if

   up_part => part%parent_part
   do while (up_part%depth > top_part%depth)
      up_part => up_part%parent_part
   end do

   if (associated(up_part, top_part)) then
      isdescendant = .true.
   else
      isdescendant = .false.
   end if
end function

function new_root_part() result(new_part)
   type(part_node_t), pointer :: new_part

   new_part => new_bare_part()
   new_part%depth = 0
   new_part%global_index = 1
   
   ! Allocate counters for root and assign initial values
   allocate(new_part%total_parts)
   allocate(new_part%total_items1)
   allocate(new_part%total_items2)

   new_part%total_parts = 1
   new_part%total_items1 = 0  ! No items initially
   new_part%total_items2 = 0  ! No items initially
end function

function new_bare_part() result(new_part)
   type(part_node_t), pointer :: new_part

   allocate (new_part)

   new_part%global_index = 0  ! Will be set when added to tree
   new_part%num_items1 = 0
   new_part%num_items2 = 0
   new_part%num_children = 0
   new_part%total_parts => null()
   new_part%total_items1 => null()
   new_part%total_items2 => null()
   new_part%parent_part => null()
   new_part%first_item1 => null()
   new_part%last_item1 => null()
   new_part%first_item2 => null()
   new_part%last_item2 => null()
   new_part%first_child_part => null()
   new_part%last_child_part => null()
   new_part%next_sibling_part => null()
   new_part%signature = [part_nodeptr_t::]
end function

subroutine link_part(link, part)
   type(link_node_t), target, intent(inout) :: link
   type(part_node_t), target, intent(inout) :: part
   type(partref_node_t), pointer :: newref

   allocate(newref)
   newref%part => part
   newref%nextref => null()
   
   ! Assign global index and increment counter
   link%total_partrefs = link%total_partrefs + 1
   newref%global_index = link%total_partrefs

   if (.not. associated(link%first_partref)) then
      link%first_partref => newref
   else
      link%last_partref%nextref => newref
   end if
   link%last_partref => newref
   link%num_parts = link%num_parts + 1
end subroutine

function new_child_part(parent_part) result(new_part)
   type(part_node_t), pointer, intent(inout) :: parent_part
   type(part_node_t), pointer :: new_part

   ! Create new part with correct depth
   new_part => new_bare_part()
   new_part%depth = parent_part%depth + 1

   ! Set parent BEFORE setting up relationships
   new_part%parent_part => parent_part
   new_part%next_sibling_part => null()

   ! Point to parent's counters
   new_part%total_parts => parent_part%total_parts
   new_part%total_items1 => parent_part%total_items1
   new_part%total_items2 => parent_part%total_items2
   
   ! Assign global index and increment counter
   new_part%total_parts = new_part%total_parts + 1
   new_part%global_index = new_part%total_parts

   ! Set up parent-child relationships (NO link association)
   if (.not. associated(parent_part%first_child_part)) then
      parent_part%first_child_part => new_part
   else
      parent_part%last_child_part%next_sibling_part => new_part
   end if
   parent_part%last_child_part => new_part
   parent_part%num_children = parent_part%num_children + 1
end function

subroutine add_new_item1(part, value)
! Add item to part without updating link itemdir (for temporary children)
   type(part_node_t), target, intent(inout) :: part
   integer, intent(in) :: value
   type(item_node_t), pointer :: new_item

   allocate(new_item)
   new_item%value = value
   new_item%next_item => null()
   
   ! Assign global index and increment counter
   part%total_items1 = part%total_items1 + 1
   new_item%global_index = part%total_items1

   if (.not. associated(part%first_item1)) then
      part%first_item1 => new_item
   else
      part%last_item1%next_item => new_item
   end if
   part%last_item1 => new_item
   part%num_items1 = part%num_items1 + 1
end subroutine

subroutine add_new_item2(part, value)
! Add item to part without updating link itemdir (for temporary children)
   type(part_node_t), target, intent(inout) :: part
   integer, intent(in) :: value
   type(item_node_t), pointer :: new_item

   allocate(new_item)
   new_item%value = value
   new_item%next_item => null()
   
   ! Assign global index and increment counter
   part%total_items2 = part%total_items2 + 1
   new_item%global_index = part%total_items2

   if (.not. associated(part%first_item2)) then
      part%first_item2 => new_item
   else
      part%last_item2%next_item => new_item
   end if
   part%last_item2 => new_item
   part%num_items2 = part%num_items2 + 1
end subroutine

subroutine move_first_item1(src, dest)
! Move first item from src to dest without updating link itemdir
   type(part_node_t), intent(inout) :: src
   type(part_node_t), target, intent(inout) :: dest
   type(item_node_t), pointer :: item_to_move

   item_to_move => src%first_item1
   if ((.not. associated(item_to_move))) error stop

   ! Remove from source
   src%first_item1 => item_to_move%next_item
   if (.not. associated(src%first_item1)) then
      src%last_item1 => null()
   end if
   src%num_items1 = src%num_items1 - 1

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

subroutine move_first_item2(src, dest)
! Move first item from src to dest without updating link itemdir
   type(part_node_t), intent(inout) :: src
   type(part_node_t), target, intent(inout) :: dest
   type(item_node_t), pointer :: item_to_move

   item_to_move => src%first_item2
   if ((.not. associated(item_to_move))) error stop

   ! Remove from source
   src%first_item2 => item_to_move%next_item
   if (.not. associated(src%first_item2)) then
      src%last_item2 => null()
   end if
   src%num_items2 = src%num_items2 - 1

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

subroutine move_part_items(src, dest)
! Move all items from src to dest without updating link itemdir
   type(part_node_t), intent(inout) :: src, dest

   do while (associated(src%first_item1))
      call move_first_item1(src, dest)
   end do

   do while (associated(src%first_item2))
      call move_first_item2(src, dest)
   end do
end subroutine

subroutine delete_chain(chain_root)
   type(chain_node_t), pointer, intent(inout) :: chain_root
   type(link_node_t), pointer :: link, next_link

   if ((.not. associated(chain_root))) error stop

   ! Delete all links (partitions only - parts are preserved)
   link => chain_root%first_link
   do while (associated(link))
      next_link => link%next_link
      call delete_link(link)
      link => next_link
   end do

   ! Then deallocate the chain root
   deallocate(chain_root)
   chain_root => null()
end subroutine

subroutine delete_link(link)
   type(link_node_t), pointer, intent(inout) :: link
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

subroutine delete_part_tree(root_part)
! Deletes an entire part tree starting from the root
   type(part_node_t), pointer, intent(inout) :: root_part

   if (.not. associated(root_part)) error stop

   ! Recursively delete all children first
   call delete_part_children(root_part)

   ! Deallocate counters
   deallocate(root_part%total_parts)
   deallocate(root_part%total_items1)
   deallocate(root_part%total_items2)

   ! Then delete the root part itself
   call delete_part(root_part)
end subroutine

recursive subroutine delete_part_children(parent_part)
! Deletes all children of a part recursively
   type(part_node_t), pointer, intent(in) :: parent_part
   type(part_node_t), pointer :: child_part, next_child

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
   type(part_node_t), pointer, intent(inout) :: part_node

   if ((.not. associated(part_node))) error stop

   ! Deallocate items
   call deallocate_items(part_node%first_item1)
   call deallocate_items(part_node%first_item2)

   ! Deallocate signature array if allocated
   if (allocated(part_node%signature)) then
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
   type(link_node_t), target, intent(in) :: link
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
         partition%parts(i)%items1(j) = item%value
         partition%itemdir1(item%value) = i
         item => item%next_item
      end do

      item => partref%part%first_item2
      do j = 1, partref%part%num_items2
         partition%parts(i)%items2(j) = item%value
         partition%itemdir2(item%value) = i
         item => item%next_item
      end do

      partref => partref%nextref
      i = i + 1
   end do
end subroutine

subroutine print_part_items(part)
   type(part_node_t), pointer, intent(in) :: part
   type(item_node_t), pointer :: item

   item => part%first_item1
   do while (associated(item))
      write(stderr, '(1X,I0)', advance='no') item%value
      item => item%next_item
   end do

   write(stderr, '(A)', advance='no') ' /'

   item => part%first_item2
   do while (associated(item))
      write(stderr, '(1X,I0)', advance='no') item%value
      item => item%next_item
   end do

   write(stderr, *)
end subroutine

subroutine print_link_itemdir(link)
   type(link_node_t), intent(in) :: link
   integer :: i

   write(stderr, *)
   write(stderr,'(A)') "Item Directory 1:"
   do i = 1, size(link%itemdir1)
      write(stderr,'(2X,I0,1X,A,1X,A)') i, "->", address(link%itemdir1(i)%ptr)
   end do

   write(stderr, *)
   write(stderr,'(A)') "Item Directory 2:"
   do i = 1, size(link%itemdir2)
      write(stderr,'(2X,I0,1X,A,1X,A)') i, "->", address(link%itemdir2(i)%ptr)
   end do
end subroutine

function new_bare_chain(tot_items1, tot_items2) result(new_chain)
   integer, intent(in) :: tot_items1, tot_items2
   type(chain_node_t), pointer :: new_chain

   allocate(new_chain)
   new_chain%num_links = 0
   new_chain%num_children = 0
   new_chain%global_index = 0  ! Will be set when added to tree
   new_chain%tot_items1 = tot_items1
   new_chain%tot_items2 = tot_items2
   new_chain%total_chains => null()
   new_chain%total_links => null()
   new_chain%total_partrefs => null()
   new_chain%first_child_chain => null()
   new_chain%last_child_chain => null()
   new_chain%next_sibling_chain => null()
   new_chain%first_link => null()
   new_chain%last_link => null()
end function

function new_root_chain(tot_items1, tot_items2) result(new_chain)
   integer, intent(in) :: tot_items1, tot_items2
   type(chain_node_t), pointer :: new_chain

   new_chain => new_bare_chain(tot_items1, tot_items2)
   new_chain%parent_chain => null()
   new_chain%split_part => null()
   new_chain%global_index = 1
   
   ! Allocate counters for root and assign initial values
   allocate(new_chain%total_chains)
   allocate(new_chain%total_links)
   allocate(new_chain%total_partrefs)
   new_chain%total_chains = 1
   new_chain%total_links = 0    ! No links initially
   new_chain%total_partrefs = 0 ! No partrefs initially
end function

function new_bare_link() result(new_link)
   type(link_node_t), pointer :: new_link

   allocate(new_link)
   new_link%num_parts = 0
   new_link%global_index = 0  ! Will be set when added to chain
   new_link%total_partrefs => null()
   new_link%first_partref => null()
   new_link%last_partref => null()
   new_link%next_link => null()
end function

function new_chain_link(chain) result(new_link)
   type(chain_node_t), target, intent(inout) :: chain
   type(link_node_t), pointer :: new_link
   integer :: i

   ! Use the initialization function
   new_link => new_bare_link()
   
   ! Point to chain's partref counter
   new_link%total_partrefs => chain%total_partrefs
   
   ! Assign global index and increment counter
   chain%total_links = chain%total_links + 1
   new_link%global_index = chain%total_links

   ! Allocate item directories
   allocate(new_link%itemdir1(chain%tot_items1))
   allocate(new_link%itemdir2(chain%tot_items2))

   ! Initialize all pointers to null
   do i = 1, chain%tot_items1
      new_link%itemdir1(i)%ptr => null()
   end do
   do i = 1, chain%tot_items2
      new_link%itemdir2(i)%ptr => null()
   end do

   ! Add link to chain
   if (.not. associated(chain%first_link)) then
      chain%first_link => new_link
   else
      chain%last_link%next_link => new_link
   end if
   chain%last_link => new_link
   chain%num_links = chain%num_links + 1
end function

function find_child_part(part, signature) result(child_part)
   type(part_node_t), intent(in) :: part
   type(part_nodeptr_t), dimension(:), intent(in) :: signature
   type(part_node_t), pointer :: child_part

   child_part => part%first_child_part
   do while (associated(child_part))
      if (child_part%signature .equiv. signature) then
         return
      end if
      child_part => child_part%next_sibling_part
   end do

   child_part => null()
end function

subroutine print_part_tree(root_part)
   type(part_node_t), pointer, intent(in) :: root_part
   logical, dimension(:), allocatable :: is_last_child

   if (.not. associated(root_part)) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "       Part Tree"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Allocate tracking array for tree lines (max depth 100)
   allocate(is_last_child(100))
   is_last_child = .false.

   ! Print root line
   write(stderr, '(A)') 'ROOT'

   ! Print children recursively
   call print_part_recursive(root_part, 0, is_last_child)

   deallocate(is_last_child)
   write(stderr, *)
end subroutine

recursive subroutine print_part_recursive(part, depth, is_last_child)
   type(part_node_t), pointer, intent(in) :: part
   integer, intent(in) :: depth
   logical, dimension(:), intent(inout) :: is_last_child
   type(part_node_t), pointer :: child_part, next_child
   integer :: i, pos
   character(len=200) :: prefix

   if (.not. associated(part)) return

   ! Process all children
   child_part => part%first_child_part
   do while (associated(child_part))
      ! Check if this is the last child
      next_child => child_part%next_sibling_part
      is_last_child(depth + 1) = .not. associated(next_child)

      ! Build prefix for this level
      prefix = " "
      pos = 2
      do i = 1, depth
         if (is_last_child(i)) then
            prefix(pos:pos+3) = "    "
         else
            prefix(pos:pos+3) = "|   "
         end if
         pos = pos + 4
      end do

      ! Add branch characters (removed trailing space)
      if (is_last_child(depth + 1)) then
         prefix(pos:pos+2) = "`--"
      else
         prefix(pos:pos+2) = "|--"
      end if
      pos = pos + 3

      ! Print part address with item counts
      write(stderr, '(A,A,1X,A,I0,A,I0,A)') prefix(1:pos-1), address(child_part), &
         '(', child_part%num_items1, '/', child_part%num_items2, ')'

      ! Recursively print this child's children
      call print_part_recursive(child_part, depth + 1, is_last_child)

      child_part => next_child
   end do
end subroutine

subroutine init_chain_from_link(input_link, chain_root, root_part)
   type(link_node_t), intent(in) :: input_link
   type(chain_node_t), pointer, intent(out) :: chain_root
   type(part_node_t), pointer, intent(out) :: root_part
   type(link_node_t), pointer :: first_link
   type(part_node_t), pointer :: new_part
   type(partref_node_t), pointer :: partref
   type(item_node_t), pointer :: item

   ! Create root part (decoupled from chain)
   root_part => new_root_part()

   ! Create root chain using itemdir sizes from input link
   chain_root => new_root_chain(size(input_link%itemdir1), size(input_link%itemdir2))

   ! Create first link
   first_link => new_chain_link(chain_root)

   ! Process all parts referenced in the input link
   partref => input_link%first_partref
   do while (associated(partref))
      ! Create new part as child of root_part
      new_part => new_child_part(root_part)

      ! Copy items from the existing part to the new part
      item => partref%part%first_item1
      do while (associated(item))
         call add_new_item1(new_part, item%value)
         item => item%next_item
      end do

      item => partref%part%first_item2
      do while (associated(item))
         call add_new_item2(new_part, item%value)
         item => item%next_item
      end do

      ! Link new part to first link (bulk itemdir update)
      call link_part(first_link, new_part)

      partref => partref%nextref
   end do
end subroutine

subroutine init_chain_from_partition(partition, chain_root, root_part)
   type(partition_t), intent(in) :: partition
   type(chain_node_t), pointer, intent(out) :: chain_root
   type(part_node_t), pointer, intent(out) :: root_part
   type(link_node_t), pointer :: first_link
   type(part_node_t), pointer :: new_part
   integer :: i, j

   ! Create root part (decoupled from chain)
   root_part => new_root_part()

   ! Create root chain
   chain_root => new_root_chain(size(partition%itemdir1), size(partition%itemdir2))

   ! Create first link
   first_link => new_chain_link(chain_root)

   ! Create parts as children of root_part and add them to the first link
   do i = 1, partition%num_parts
      ! Create new part as child of root_part
      new_part => new_child_part(root_part)

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
end subroutine

function first_partition(partition) result(semipartition)
   type(partition_t), intent(in) :: partition
   type(semipartition_t) :: semipartition
   integer :: i

   semipartition%num_parts = partition%num_parts
   allocate(semipartition%parts(semipartition%num_parts))
   allocate(semipartition%itemdir(size(partition%itemdir1)))
   semipartition%itemdir = partition%itemdir1

   do i = 1, semipartition%num_parts
      semipartition%parts(i)%num_items = partition%parts(i)%num_items1
      allocate(semipartition%parts(i)%items(semipartition%parts(i)%num_items))
      semipartition%parts(i)%items = partition%parts(i)%items1
   end do
end function

function second_partition(partition) result(semipartition)
   type(partition_t), intent(in) :: partition
   type(semipartition_t) :: semipartition
   integer :: i

   semipartition%num_parts = partition%num_parts
   allocate(semipartition%parts(semipartition%num_parts))
   allocate(semipartition%itemdir(size(partition%itemdir2)))
   semipartition%itemdir = partition%itemdir2

   do i = 1, semipartition%num_parts
      semipartition%parts(i)%num_items = partition%parts(i)%num_items2
      allocate(semipartition%parts(i)%items(semipartition%parts(i)%num_items))
      semipartition%parts(i)%items = partition%parts(i)%items2
   end do
end function

function new_child_chain(chain, split_part) result(new_chain)
   type(chain_node_t), pointer, intent(inout) :: chain
   type(part_node_t), pointer, intent(in) :: split_part
   type(chain_node_t), pointer :: new_chain

   ! Create new child chain
   new_chain => new_bare_chain(chain%tot_items1, chain%tot_items2)
   new_chain%split_part => split_part

   ! Set up parent-child relationship
   new_chain%parent_chain => chain
   
   ! Point to parent's counters
   new_chain%total_chains => chain%total_chains
   new_chain%total_links => chain%total_links
   new_chain%total_partrefs => chain%total_partrefs
   
   ! Assign global index and increment counter
   new_chain%total_chains = new_chain%total_chains + 1
   new_chain%global_index = new_chain%total_chains

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
   type(link_node_t), target, intent(inout) :: link
   type(part_node_t), target, intent(in) :: part
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
   type(link_node_t), target, intent(inout) :: link
   type(part_node_t), target, intent(in) :: part
   type(item_node_t), pointer :: item

   ! Register all items from molecule 1
   item => part%first_item1
   do while (associated(item))
      link%itemdir1(item%value)%ptr => part
      item => item%next_item
   end do

   ! Register all items from molecule 2
   item => part%first_item2
   do while (associated(item))
      link%itemdir2(item%value)%ptr => part
      item => item%next_item
   end do
end subroutine

subroutine print_chain_tree(root_chain)
   type(chain_node_t), pointer, intent(in) :: root_chain
   logical, dimension(:), allocatable :: is_last_child

   if (.not. associated(root_chain)) error stop

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "       Split Tree"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Allocate tracking array for tree lines (max depth 100)
   allocate(is_last_child(100))
   is_last_child = .false.

   ! Print children recursively
   write(stderr, '(A)') 'ROOT'
   call print_chain_recursive(root_chain, 0, is_last_child)

   deallocate(is_last_child)
   write(stderr, *)
end subroutine

recursive subroutine print_chain_recursive(chain, depth, is_last_child)
   type(chain_node_t), pointer, intent(in) :: chain
   integer, intent(in) :: depth
   logical, dimension(:), intent(inout) :: is_last_child
   type(chain_node_t), pointer :: child_chain, next_child
   integer :: i, pos
   character(len=200) :: prefix

   if (.not. associated(chain)) return

   ! Process all children
   child_chain => chain%first_child_chain
   do while (associated(child_chain))
      ! Check if this is the last child
      next_child => child_chain%next_sibling_chain
      is_last_child(depth + 1) = .not. associated(next_child)

      ! Build prefix for this level
      prefix = " "
      pos = 2
      do i = 1, depth
         if (is_last_child(i)) then
            prefix(pos:pos+3) = "    "
         else
            prefix(pos:pos+3) = "|   "
         end if
         pos = pos + 4
      end do

      ! Add branch characters (removed trailing space)
      if (is_last_child(depth + 1)) then
         prefix(pos:pos+2) = "`--"
      else
         prefix(pos:pos+2) = "|--"
      end if
      pos = pos + 3

      ! Print the child address with item counts
      write(stderr, '(A,A,1X,A,I0,A,I0,A)') prefix(1:pos-1), address(child_chain%split_part), &
         '(', child_chain%split_part%num_items1, '/', child_chain%split_part%num_items2, ')'

      ! Recursively print this child's children
      call print_chain_recursive(child_chain, depth + 1, is_last_child)

      child_chain => next_child
   end do
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
subroutine print_tree_signatures(root_part)
   type(part_node_t), pointer, intent(in) :: root_part

   if (.not. associated(root_part)) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "    Part Signatures"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Print children recursively
   call print_signature_recursive(root_part)

   write(stderr, *)
end subroutine

! Helper recursive procedure for print_tree_signatures
recursive subroutine print_signature_recursive(part)
   type(part_node_t), pointer, intent(in) :: part
   type(part_node_t), pointer :: child_part

   if (.not. associated(part)) return

   ! Process all children in the same order as print_part_tree
   child_part => part%first_child_part
   do while (associated(child_part))
      ! Print the child signature
      write(stderr,'(A)',advance='no') address(child_part) // ':'
      call print_part_signature(child_part%signature)

      ! Recursively print this child's children
      call print_signature_recursive(child_part)

      child_part => child_part%next_sibling_part
   end do
end subroutine

subroutine print_tree_items(root_part)
   type(part_node_t), pointer, intent(in) :: root_part

   if (.not. associated(root_part)) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "      Part Items"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Print children recursively
   call print_items_recursive(root_part)

   write(stderr, *)
end subroutine

recursive subroutine print_items_recursive(part)
   type(part_node_t), pointer, intent(in) :: part
   type(part_node_t), pointer :: child_part

   if (.not. associated(part)) return

   ! Process all children in the same order as print_part_tree
   child_part => part%first_child_part
   do while (associated(child_part))
      ! Print the child items with address prefix
      write(stderr, '(A)', advance='no') address(child_part) // ':'
      call print_part_items(child_part)

      ! Recursively print this child's children
      call print_items_recursive(child_part)

      child_part => child_part%next_sibling_part
   end do
end subroutine

subroutine print_leaf_items(root_part)
   type(part_node_t), pointer, intent(in) :: root_part

   if (.not. associated(root_part)) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "      Leaf Items"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Print children recursively
   call print_leaf_items_recursive(root_part)

   write(stderr, *)
end subroutine

recursive subroutine print_leaf_items_recursive(part)
   type(part_node_t), pointer, intent(in) :: part
   type(part_node_t), pointer :: child_part

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
      call print_leaf_items_recursive(child_part)

      child_part => child_part%next_sibling_part
   end do
end subroutine

! New procedure to print global indices
subroutine print_part_indices(root_part)
   type(part_node_t), pointer, intent(in) :: root_part

   if (.not. associated(root_part)) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 35)
   write(stderr, '(A)') "        Part Indices"
   write(stderr, '(A)') repeat("=", 35)
   write(stderr, *)
   write(stderr, '(A,I0)') "Total parts created: ", root_part%total_parts
   write(stderr, '(A,I0)') "Total items1 created: ", root_part%total_items1
   write(stderr, '(A,I0)') "Total items2 created: ", root_part%total_items2
   write(stderr, *)

   ! Print children recursively
   call print_part_indices_recursive(root_part)

   write(stderr, *)
end subroutine

! Helper recursive procedure for print_part_indices
recursive subroutine print_part_indices_recursive(part)
   type(part_node_t), pointer, intent(in) :: part
   type(part_node_t), pointer :: child_part
   type(item_node_t), pointer :: item

   if (.not. associated(part)) return

   ! Process all children
   child_part => part%first_child_part
   do while (associated(child_part))
      ! Print the child part index
      write(stderr, '(A,I0,A,A,A)') "Part ", child_part%global_index, " (address: ", address(child_part), ")"
      
      ! Print items in this part
      item => child_part%first_item1
      do while (associated(item))
         write(stderr, '(A,I0,A,I0,A)') "  Item1 ", item%global_index, " (value: ", item%value, ")"
         item => item%next_item
      end do
      
      item => child_part%first_item2
      do while (associated(item))
         write(stderr, '(A,I0,A,I0,A)') "  Item2 ", item%global_index, " (value: ", item%value, ")"
         item => item%next_item
      end do

      ! Recursively print this child's children
      call print_part_indices_recursive(child_part)

      child_part => child_part%next_sibling_part
   end do
end subroutine

! New procedure to print global indices for chain tree
subroutine print_chain_indices(root_chain)
   type(chain_node_t), pointer, intent(in) :: root_chain

   if (.not. associated(root_chain)) then
      write(stderr, '(A)') "Chain tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 35)
   write(stderr, '(A)') "     Chain Indices"
   write(stderr, '(A)') repeat("=", 35)
   write(stderr, *)
   write(stderr, '(A,I0)') "Total chains created: ", root_chain%total_chains
   write(stderr, '(A,I0)') "Total links created: ", root_chain%total_links
   write(stderr, '(A,I0)') "Total partrefs created: ", root_chain%total_partrefs
   write(stderr, *)

   ! Print children recursively
   call print_chain_indices_recursive(root_chain)

   write(stderr, *)
end subroutine

! Helper recursive procedure for print_chain_indices
recursive subroutine print_chain_indices_recursive(chain)
   type(chain_node_t), pointer, intent(in) :: chain
   type(chain_node_t), pointer :: child_chain
   type(link_node_t), pointer :: link
   type(partref_node_t), pointer :: partref

   if (.not. associated(chain)) return

   ! Print links in this chain
   link => chain%first_link
   do while (associated(link))
      write(stderr, '(A,I0,A,I0,A)') "  Link ", link%global_index, " (", link%num_parts, " parts)"
      
      ! Print partrefs in this link
      partref => link%first_partref
      do while (associated(partref))
         write(stderr, '(A,I0,A,A,A)') "    Partref ", partref%global_index, " (part: ", address(partref%part), ")"
         partref => partref%nextref
      end do
      
      link => link%next_link
   end do

   ! Process all child chains
   child_chain => chain%first_child_chain
   do while (associated(child_chain))
      ! Print the child chain index
      write(stderr, '(A,I0,A,A,A)') "Chain ", child_chain%global_index, " (split part: ", address(child_chain%split_part), ")"

      ! Recursively print this child's content and children
      call print_chain_indices_recursive(child_chain)

      child_chain => child_chain%next_sibling_chain
   end do
end subroutine

end module
