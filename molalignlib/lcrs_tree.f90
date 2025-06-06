module lcrs_tree
use parameters
implicit none
private

! Item node
type, public :: item_node_t
   integer :: value
   type(item_node_t), pointer :: next_item
end type

! LCRS branch node
type, public :: split_node_t
   integer :: num_links
   integer :: num_children
   integer :: tot_items1
   integer :: tot_items2
   type(split_node_t), pointer :: parent_branch
   type(split_node_t), pointer :: first_child_branch
   type(split_node_t), pointer :: last_child_branch
   type(split_node_t), pointer :: next_sibling_branch
   type(link_node_t), pointer :: first_link
   type(link_node_t), pointer :: last_link
   type(part_node_t), pointer :: split_part
end type

! Part reference node
type, public :: partref_node_t
   type(part_node_t), pointer :: part
   type(partref_node_t), pointer :: nextref
   type(link_node_t), pointer :: parent_link
end type

! Link node
type, public :: link_node_t
   integer :: num_parts
   type(link_node_t), pointer :: next_link
   type(partref_node_t), pointer :: first_partref
   type(partref_node_t), pointer :: last_partref
   type(part_nodeptr_t), dimension(:), allocatable :: itemdir1
   type(part_nodeptr_t), dimension(:), allocatable :: itemdir2
end type

! LCRS part node
type, public :: part_node_t
   integer :: num_items1
   integer :: num_items2
   integer :: num_children
   integer :: depth
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
type, public :: semipartarray_t
   integer :: num_items
   integer, dimension(:), allocatable :: items
end type

! Semipartition array
type, public :: semipartitionarray_t
   integer :: num_parts
   type(semipartarray_t), dimension(:), allocatable :: parts
   integer, dimension(:), allocatable :: itemdir
end type

! Part array
type, public :: partarray_t
   integer :: num_children
   integer :: num_items1
   integer :: num_items2
   integer, dimension(:), allocatable :: children
   integer, dimension(:), allocatable :: items1
   integer, dimension(:), allocatable :: items2
   integer, dimension(:), allocatable :: signature
end type

! Partition array
type, public :: partitionarray_t
   integer :: num_parts
   type(partarray_t), dimension(:), allocatable :: parts
   integer, dimension(:), allocatable :: itemdir1
   integer, dimension(:), allocatable :: itemdir2
end type

! Partchain array
type, public :: chainarray_t
   integer :: num_links
   integer :: tot_items1
   integer :: tot_items2
   type(partitionarray_t), dimension(:), allocatable :: partitions
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
public new_root_branch
public new_child_branch
public new_root_link
public new_generic_link
public new_root_part
public new_child_part
public link_part
public update_itemdir
public remove_last_link
public remove_onlychild_part
public add_new_item1
public add_new_item2
public init_chain_from_link
public move_first_item1
public move_first_item2
public move_part_items
public delete_chain
public print_part_tree
public print_chain
public print_chainarray
public print_part_items
public print_tree_items
public print_leaf_items
public print_link_itemdir
public print_part_signature
public print_tree_signatures
public sort_parts_by_size
public partition_to_partitionarray
public init_chain_from_partarray
public find_child_part
public first_partition
public second_partition
public add_branch_part
public print_split_tree
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

function new_root_part(depth) result(new_part)
   integer, intent(in) :: depth
   type(part_node_t), pointer :: new_part

   allocate (new_part)

   new_part%num_children = 0
   new_part%num_items1 = 0
   new_part%num_items2 = 0
   new_part%depth = depth
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
   newref%parent_link => link

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
   new_part => new_root_part(parent_part%depth + 1)

   ! Set parent BEFORE setting up relationships
   new_part%parent_part => parent_part
   new_part%next_sibling_part => null()

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
   type(split_node_t), pointer, intent(inout) :: chain_root
   type(link_node_t), pointer :: link, next_link

   if ((.not. associated(chain_root))) error stop

   ! Delete all links (partitions only - parts are preserved)
   link => chain_root%first_link
   do while (associated(link))
      next_link => link%next_link
      call delete_partition(link)
      link => next_link
   end do

   ! Then deallocate the chain root
   deallocate(chain_root)
   chain_root => null()
end subroutine

subroutine delete_partition(parent_link)
   type(link_node_t), pointer, intent(inout) :: parent_link
   type(partref_node_t), pointer :: partref, nextref

   if ((.not. associated(parent_link))) error stop

   ! Delete all part references WITHOUT deleting the parts themselves
   partref => parent_link%first_partref
   do while (associated(partref))
      nextref => partref%nextref

      ! Only delete the part reference, NOT the part itself
      deallocate(partref)
      partref => nextref
   end do

   ! Deallocate directories (always allocated)
   deallocate(parent_link%itemdir1)
   deallocate(parent_link%itemdir2)

   ! Deallocate partition root
   deallocate(parent_link)
   parent_link => null()
end subroutine

subroutine delete_part_tree(root_part)
! Deletes an entire part tree starting from the root
   type(part_node_t), pointer, intent(inout) :: root_part

   if (.not. associated(root_part)) error stop

   ! Recursively delete all children first
   call delete_part_children(root_part)

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

subroutine partition_to_partitionarray(parent_link, partition)
   type(link_node_t), target, intent(in) :: parent_link
   type(partitionarray_t), intent(out) :: partition
   type(partref_node_t), pointer :: partref
   type(item_node_t), pointer :: item
   integer :: i, j

   partition%num_parts = parent_link%num_parts
   allocate(partition%parts(partition%num_parts))
   allocate(partition%itemdir1(size(parent_link%itemdir1)))
   allocate(partition%itemdir2(size(parent_link%itemdir2)))

   partref => parent_link%first_partref
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

function new_root_branch(tot_items1, tot_items2) result(new_branch)
   integer, intent(in) :: tot_items1, tot_items2
   type(split_node_t), pointer :: new_branch

   allocate(new_branch)
   new_branch%num_links = 0
   new_branch%num_children = 0
   new_branch%tot_items1 = tot_items1
   new_branch%tot_items2 = tot_items2
   new_branch%split_part => null()
   new_branch%parent_branch => null()
   new_branch%first_child_branch => null()
   new_branch%last_child_branch => null()
   new_branch%next_sibling_branch => null()
   new_branch%first_link => null()
   new_branch%last_link => null()
end function

function new_root_link() result(new_link)
   type(link_node_t), pointer :: new_link

   allocate(new_link)
   new_link%num_parts = 0
   new_link%first_partref => null()
   new_link%last_partref => null()
   new_link%next_link => null()
end function

function new_generic_link(branch) result(new_link)
   type(split_node_t), target, intent(inout) :: branch
   type(link_node_t), pointer :: new_link
   integer :: i

   ! Use the initialization function
   new_link => new_root_link()

   ! Allocate item directories
   allocate(new_link%itemdir1(branch%tot_items1))
   allocate(new_link%itemdir2(branch%tot_items2))

   ! Initialize all pointers to null
   do i = 1, branch%tot_items1
      new_link%itemdir1(i)%ptr => null()
   end do
   do i = 1, branch%tot_items2
      new_link%itemdir2(i)%ptr => null()
   end do

   ! Add link to branch
   if (.not. associated(branch%first_link)) then
      branch%first_link => new_link
   else
      branch%last_link%next_link => new_link
   end if
   branch%last_link => new_link
   branch%num_links = branch%num_links + 1
end function

subroutine remove_last_link(branch)
   type(split_node_t), pointer, intent(inout) :: branch
   type(link_node_t), pointer :: link_to_remove, prev_link

   if (.not. associated(branch%last_link)) then
      error stop "Cannot remove last link: no links in branch"
   end if

   link_to_remove => branch%last_link

   ! If this is the only link
   if (associated(branch%first_link, branch%last_link)) then
      branch%first_link => null()
      branch%last_link => null()
   else
      ! Find the previous link
      prev_link => branch%first_link
      do while (.not. associated(prev_link%next_link, link_to_remove))
         prev_link => prev_link%next_link
      end do

      ! Update pointers
      prev_link%next_link => null()
      branch%last_link => prev_link
   end if

   branch%num_links = branch%num_links - 1

   ! Delete the removed link
   call delete_partition(link_to_remove)
end subroutine

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

subroutine print_chain(chain_root)
   type(split_node_t), pointer, intent(in) :: chain_root
   type(link_node_t), pointer :: link
   integer :: link_idx

   if (.not. associated(chain_root)) error stop

   write(stderr,*)
   write(stderr,'(A)') repeat("=", 17)
   write(stderr,'(A)') " Partition chain"
   write(stderr,'(A)') repeat("=", 17)

   ! Print all partitions in this chain
   link => chain_root%first_link
   link_idx = 1
   do while (associated(link))
      call print_partition(link, link_idx)
      link => link%next_link
      link_idx = link_idx + 1
   end do

   write(stderr,*)  ! Final newline
end subroutine

subroutine print_partition(parent_link, link_idx)
   type(link_node_t), pointer, intent(in) :: parent_link
   integer, intent(in) :: link_idx

   if (.not. associated(parent_link)) error stop

   write(stderr,*)
   write(stderr,'(A,I0,A)') "( Link ", link_idx, " )"

   ! Print part items
   call print_partition_items(parent_link)

   ! Print part signatures
   call print_partition_signatures(parent_link)

   ! Print part children
   call print_partition_children(parent_link)
end subroutine

subroutine print_partition_items(parent_link)
   type(link_node_t), pointer, intent(in) :: parent_link
   type(partref_node_t), pointer :: partref

   if (.not. associated(parent_link)) error stop

   write(stderr, *)
   write(stderr,'(A)') "Items"
   write(stderr,'(A)') repeat("-", 6)

   partref => parent_link%first_partref
   do while (associated(partref))
      write(stderr, '(A)', advance='no') address(partref%part) // ':'
      call print_part_items(partref%part)
      partref => partref%nextref
   end do
end subroutine

subroutine print_partition_children(parent_link)
   type(link_node_t), pointer, intent(in) :: parent_link
   type(partref_node_t), pointer :: partref
   type(part_node_t), pointer :: child_part

   if (.not. associated(parent_link)) error stop

   write(stderr,*)
   write(stderr,'(A)') "Children"
   write(stderr,'(A)') repeat("-", 9)

   partref => parent_link%first_partref
   do while (associated(partref))
      ! Print part address
      write(stderr,'(A)',advance='no') address(partref%part) // ":"

      ! Print part children
      child_part => partref%part%first_child_part
      do while (associated(child_part))
         ! Print which part in next partition this child_part points to
         write(stderr,'(1X,A)',advance='no') address(child_part)
         child_part => child_part%next_sibling_part
      end do
      write(stderr,*)

      partref => partref%nextref
   end do
end subroutine

subroutine print_partition_signatures(parent_link)
   type(link_node_t), pointer, intent(in) :: parent_link
   type(partref_node_t), pointer :: partref

   if (.not. associated(parent_link)) error stop

   write(stderr,*)
   write(stderr,'(A)') "Signatures"
   write(stderr,'(A)') repeat("-", 11)

   partref => parent_link%first_partref
   do while (associated(partref))
      ! Print part address
      write(stderr,'(A)',advance='no') address(partref%part) // ":"

      ! Print part signature using the new procedure
      call print_part_signature(partref%part%signature)

      partref => partref%nextref
   end do
end subroutine

subroutine print_chainarray(chainarray)
   type(chainarray_t), intent(in) :: chainarray
   integer :: link_idx

   write(stderr,*)
   write(stderr,'(A)') repeat("=", 17)
   write(stderr,'(A)') " Partition chain"
   write(stderr,'(A)') repeat("=", 17)

   ! Print each link
   do link_idx = 1, chainarray%num_links
      call print_partitionarray(chainarray, link_idx)
   end do

   write(stderr,*)  ! Final newline
end subroutine

subroutine print_partitionarray(chainarray, link_idx)
   type(chainarray_t), intent(in) :: chainarray
   integer, intent(in) :: link_idx

   write(stderr,*)
   write(stderr,'(A,I0,A)') "( Link ", link_idx, " )"

   ! Print part items
   call print_partitionarray_items(chainarray, link_idx)

   ! Print part signatures
   call print_partitionarray_signatures(chainarray, link_idx)

   ! Print part children
   call print_partitionarray_children(chainarray, link_idx)
end subroutine

subroutine print_partitionarray_items(chainarray, link_idx)
   type(chainarray_t), intent(in) :: chainarray
   integer, intent(in) :: link_idx
   integer :: i, j

   write(stderr, *)
   write(stderr,'(A)') "Items"
   write(stderr,'(A)') repeat("-", 6)

   do i = 1, chainarray%partitions(link_idx)%num_parts
      write(stderr, '(I3,A)', advance='no') i, ':'

      do j = 1, chainarray%partitions(link_idx)%parts(i)%num_items1
         write(stderr,'(1X,I0)',advance='no') chainarray%partitions(link_idx)%parts(i)%items1(j)
      end do

      write(stderr,'(A)',advance='no') ' /'

      do j = 1, chainarray%partitions(link_idx)%parts(i)%num_items2
         write(stderr,'(1X,I0)',advance='no') chainarray%partitions(link_idx)%parts(i)%items2(j)
      end do

      write(stderr,*)
   end do
end subroutine

subroutine print_partitionarray_children(chainarray, link_idx)
   type(chainarray_t), intent(in) :: chainarray
   integer, intent(in) :: link_idx
   integer :: i, j

   write(stderr,*)
   write(stderr,'(A)') "Children"
   write(stderr,'(A)') repeat("-", 9)

   do i = 1, chainarray%partitions(link_idx)%num_parts
      write(stderr,'(I3,A)',advance='no') i, ":"
      do j = 1, chainarray%partitions(link_idx)%parts(i)%num_children
         write(stderr,'(1X,I3)',advance='no') &
            chainarray%partitions(link_idx)%parts(i)%children(j)
      end do
      write(stderr,*)
   end do
end subroutine

subroutine print_partitionarray_signatures(chainarray, link_idx)
   type(chainarray_t), intent(in) :: chainarray
   integer, intent(in) :: link_idx
   integer :: i, j

   write(stderr,*)
   write(stderr,'(A)') "Signatures"
   write(stderr,'(A)') repeat("-", 10)

   do i = 1, chainarray%partitions(link_idx)%num_parts
      write(stderr,'(I3,A)',advance='no') i, ":"
      do j = 1, size(chainarray%partitions(link_idx)%parts(i)%signature)
         write(stderr,'(1X,I3)',advance='no') &
            chainarray%partitions(link_idx)%parts(i)%signature(j)
      end do
      write(stderr,*)
   end do
end subroutine

subroutine init_chain_from_link(input_link, chain_root, root_part)
   type(link_node_t), intent(in) :: input_link
   type(split_node_t), pointer, intent(out) :: chain_root
   type(part_node_t), pointer, intent(out) :: root_part
   type(link_node_t), pointer :: first_link
   type(part_node_t), pointer :: new_part
   type(partref_node_t), pointer :: partref
   type(item_node_t), pointer :: item

   ! Create root part (decoupled from chain)
   root_part => new_root_part(0)

   ! Create root chain using itemdir sizes from input link
   chain_root => new_root_branch(size(input_link%itemdir1), size(input_link%itemdir2))

   ! Create first link
   first_link => new_generic_link(chain_root)

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

subroutine init_chain_from_partarray(partition, chain_root, root_part)
   type(partitionarray_t), intent(in) :: partition
   type(split_node_t), pointer, intent(out) :: chain_root
   type(part_node_t), pointer, intent(out) :: root_part
   type(link_node_t), pointer :: first_link
   type(part_node_t), pointer :: new_part
   integer :: i, j

   ! Create root part (decoupled from chain)
   root_part => new_root_part(0)

   ! Create root chain
   chain_root => new_root_branch(size(partition%itemdir1), size(partition%itemdir2))

   ! Create first link
   first_link => new_generic_link(chain_root)

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

function first_partition(bipartition) result(partition)
   type(partitionarray_t), intent(in) :: bipartition
   type(semipartitionarray_t) :: partition
   integer :: i

   partition%num_parts = bipartition%num_parts
   allocate(partition%parts(partition%num_parts))
   allocate(partition%itemdir(size(bipartition%itemdir1)))
   partition%itemdir = bipartition%itemdir1

   do i = 1, partition%num_parts
      partition%parts(i)%num_items = bipartition%parts(i)%num_items1
      allocate(partition%parts(i)%items(partition%parts(i)%num_items))
      partition%parts(i)%items = bipartition%parts(i)%items1
   end do
end function

function second_partition(bipartition) result(partition)
   type(partitionarray_t), intent(in) :: bipartition
   type(semipartitionarray_t) :: partition
   integer :: i

   partition%num_parts = bipartition%num_parts
   allocate(partition%parts(partition%num_parts))
   allocate(partition%itemdir(size(bipartition%itemdir2)))
   partition%itemdir = bipartition%itemdir2

   do i = 1, partition%num_parts
      partition%parts(i)%num_items = bipartition%parts(i)%num_items2
      allocate(partition%parts(i)%items(partition%parts(i)%num_items))
      partition%parts(i)%items = bipartition%parts(i)%items2
   end do
end function

subroutine sort_parts_by_size(parent_link)
   type(link_node_t), pointer, intent(inout) :: parent_link
   type(partref_node_t), pointer :: partref, nextref, prevref
   logical :: swapped

   ! Bubble sort implementation for part reference linked list
   do
      swapped = .false.
      partref => parent_link%first_partref
      prevref => null()
      do while (associated(partref) .and. associated(partref%nextref))
         nextref => partref%nextref
         ! Check if we need to swap (current has more items than next)
         if (partref%part%num_items1 > nextref%part%num_items1) then
            swapped = .true.
            ! Perform the swap
            partref%nextref => nextref%nextref
            nextref%nextref => partref
            if (associated(prevref)) then
               prevref%nextref => nextref
            else
               ! Update first_partref if we're swapping the first element
               parent_link%first_partref => nextref
            end if
            ! Update last_partref if necessary
            if (.not. associated(partref%nextref)) then
               parent_link%last_partref => partref
            end if
            ! Update prevref for next iteration
            prevref => nextref
         else
            ! No swap needed, just advance
            prevref => partref
            partref => nextref
         end if
      end do
      ! If no swaps occurred, the list is sorted
      if (.not. swapped) exit
   end do
end subroutine

function new_child_branch(branch, split_part) result(new_branch)
   type(split_node_t), pointer, intent(inout) :: branch
   type(part_node_t), pointer, intent(in) :: split_part
   type(split_node_t), pointer :: new_branch

   ! Create new child branch
   new_branch => new_root_branch(branch%tot_items1, branch%tot_items2)
   new_branch%split_part => split_part

   ! Set up parent-child relationship
   new_branch%parent_branch => branch

   ! Add as child to parent (optimized with last_child_branch pointer)
   if (.not. associated(branch%first_child_branch)) then
      ! This is the first child
      branch%first_child_branch => new_branch
      branch%last_child_branch => new_branch
   else
      ! Add as next sibling to the current last child
      branch%last_child_branch%next_sibling_branch => new_branch
      branch%last_child_branch => new_branch
   end if

   branch%num_children = branch%num_children + 1
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
   newref%parent_link => link

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

subroutine remove_onlychild_part(parent_part)
! Remove the only child part (since we know leaf parts start with 0 children)
   type(part_node_t), pointer, intent(inout) :: parent_part

   ! Assert only child
   if (parent_part%num_children /= 1) error stop

   ! Delete the child part (no itemdir cleanup needed)
   call delete_part(parent_part%first_child_part)

   ! Reset parent's child pointers (removing the only child)
   parent_part%first_child_part => null()
   parent_part%last_child_part => null()
   parent_part%num_children = 0
end subroutine

subroutine print_split_tree(root_branch)
   type(split_node_t), pointer, intent(in) :: root_branch
   logical, dimension(:), allocatable :: is_last_child

   if (.not. associated(root_branch)) error stop

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
   call print_split_recursive(root_branch, 0, is_last_child)

   deallocate(is_last_child)
   write(stderr, *)
end subroutine

recursive subroutine print_split_recursive(branch, depth, is_last_child)
   type(split_node_t), pointer, intent(in) :: branch
   integer, intent(in) :: depth
   logical, dimension(:), intent(inout) :: is_last_child
   type(split_node_t), pointer :: child_branch, next_child
   integer :: i, pos
   character(len=200) :: prefix

   if (.not. associated(branch)) return

   ! Process all children
   child_branch => branch%first_child_branch
   do while (associated(child_branch))
      ! Check if this is the last child
      next_child => child_branch%next_sibling_branch
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
      write(stderr, '(A,A,1X,A,I0,A,I0,A)') prefix(1:pos-1), address(child_branch%split_part), &
         '(', child_branch%split_part%num_items1, '/', child_branch%split_part%num_items2, ')'

      ! Recursively print this child's children
      call print_split_recursive(child_branch, depth + 1, is_last_child)

      child_branch => next_child
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

end module
