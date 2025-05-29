module lcrs_tree
use parameters
implicit none
private

! Item node
type, public :: item_node_t
   integer :: value
   type(item_node_t), pointer :: next_item
end type

! Chain root
type, public :: chain_root_t
   integer :: num_links
   integer :: tot_items1
   integer :: tot_items2
   type(link_node_t), pointer :: first_link
   type(link_node_t), pointer :: last_link
end type

! Part reference node
type, public :: partref_node_t
   type(part_node_t), pointer :: part_node
   type(partref_node_t), pointer :: next_partref
end type

! Link node
type, public :: link_node_t
   integer :: num_parts
   type(link_node_t), pointer :: next_link
   type(chain_root_t), pointer :: chain_root
   type(partref_node_t), pointer :: first_partref
   type(partref_node_t), pointer :: last_partref
   type(part_nodeptr_t), dimension(:), allocatable :: itemdir1
   type(part_nodeptr_t), dimension(:), allocatable :: itemdir2
end type

! Part node
type, public :: part_node_t
   integer :: num_items1
   integer :: num_items2
   integer :: num_children
   integer :: depth
   type(item_node_t), pointer :: first_item1
   type(item_node_t), pointer :: first_item2
   type(item_node_t), pointer :: last_item1
   type(item_node_t), pointer :: last_item2
   type(link_node_t), pointer :: parent_link
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

interface operator(==)
   module procedure part_nodeptr_equality
end interface

interface operator (.equiv.)
   module procedure signature_equivalence
end interface

! Branch node for LCRS tree with part reference list
type, public :: branch_node_t
   integer :: num_parts
   integer :: num_children
   integer :: tot_items1
   integer :: tot_items2
   type(branch_node_t), pointer :: parent_branch
   type(branch_node_t), pointer :: first_child_branch
   type(branch_node_t), pointer :: last_child_branch
   type(branch_node_t), pointer :: next_sibling_branch
   type(partref_node_t), pointer :: first_partref ! Obsolete
   type(partref_node_t), pointer :: last_partref  ! Obsolete
   type(part_node_t), pointer :: part_node
   type(link_node_t), pointer :: first_link
   type(link_node_t), pointer :: last_link
end type

! Make types and procedures public
public address
public init_chain
public add_new_link
public add_new_link_branch
public add_part
public add_new_part
public add_new_item1
public add_new_item2
public move_first_item1
public move_first_item2
public move_node_items
public delete_chain
public print_part_tree
public print_chain
public print_chainarray
public print_items
public print_itemdirs
public sort_parts_by_size
public partition_to_partitionarray
public chain_from_partitionarray
public find_child_part_node
public first_partition
public second_partition
public make_new_branch
public add_new_branch
public update_branch_parts
public remove_branch_part
public print_branch_tree
public print_branch_contents
public operator(==)
public operator(.equiv.)


contains


integer function address(nodeptr)
   use iso_c_binding, only: c_loc, c_intptr_t
   type(part_node_t), target, intent(in) :: nodeptr
   address = modulo(transfer(c_loc(nodeptr), c_intptr_t), 16**4)
end function

elemental function part_nodeptr_equality(left, right) result(equality)
   type(part_nodeptr_t), intent(in) :: left, right
   logical :: equality
   equality = associated(left%ptr, right%ptr)
end function

function signature_equivalence(array1, array2) result(equiv)
   type(part_nodeptr_t), dimension(:), intent(in) :: array1, array2
   logical :: equiv
   integer :: i, j, matches

   if (size(array1) /= size(array2)) then
      equiv = .false.
      return
   end if

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

function init_part(depth) result(new_part)
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
   new_part%parent_link => null()
   new_part%signature = [part_nodeptr_t::]
end function


subroutine add_part(link, part)
   type(link_node_t), target, intent(inout) :: link
   type(part_node_t), target, intent(inout) :: part
   type(partref_node_t), pointer :: new_partref
   type(item_node_t), pointer :: item

   part%parent_link => link

   allocate(new_partref)
   new_partref%part_node => part
   new_partref%next_partref => null()

   if (.not. associated(link%first_partref)) then
      link%first_partref => new_partref
   else
      link%last_partref%next_partref => new_partref
   end if
   link%last_partref => new_partref
   link%num_parts = link%num_parts + 1

   ! CRITICAL FIX: Update the itemdir arrays in the new link
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

function add_new_part(parent_link, part) result(new_part)
   type(link_node_t), pointer, intent(inout) :: parent_link
   type(part_node_t), pointer, intent(inout) :: part
   type(part_node_t), pointer :: new_part

   ! Create new part with correct depth
   new_part => init_part(part%depth + 1)

   ! Set parent BEFORE calling add_part
   new_part%parent_part => part
   new_part%next_sibling_part => null()

   ! Now add to partition
   call add_part(parent_link, new_part)

   ! Set up parent-child relationships
   if (.not. associated(part%first_child_part)) then
      part%first_child_part => new_part
   else
      part%last_child_part%next_sibling_part => new_part
   end if
   part%last_child_part => new_part
   part%num_children = part%num_children + 1
end function

subroutine add_new_item1(part, value)
   type(part_node_t), target, intent(inout) :: part
   integer, intent(in) :: value
   type(item_node_t), pointer :: new_item

   allocate(new_item)
   new_item%value = value
   new_item%next_item => null()
   part%parent_link%itemdir1(value)%ptr => part

   if (.not. associated(part%first_item1)) then
      part%first_item1 => new_item
   else
      part%last_item1%next_item => new_item
   end if
   part%last_item1 => new_item
   part%num_items1 = part%num_items1 + 1
end subroutine

subroutine add_new_item2(part, value)
   type(part_node_t), target, intent(inout) :: part
   integer, intent(in) :: value
   type(item_node_t), pointer :: new_item

   allocate(new_item)
   new_item%value = value
   new_item%next_item => null()
   part%parent_link%itemdir2(value)%ptr => part

   if (.not. associated(part%first_item2)) then
      part%first_item2 => new_item
   else
      part%last_item2%next_item => new_item
   end if
   part%last_item2 => new_item
   part%num_items2 = part%num_items2 + 1
end subroutine

subroutine move_first_item1(src, dest)
   type(part_node_t), intent(inout) :: src
   type(part_node_t), target, intent(inout) :: dest
   type(item_node_t), pointer :: item_to_move

   item_to_move => src%first_item1
   if ((.not. associated(item_to_move))) error stop

   src%first_item1 => item_to_move%next_item
   if (.not. associated(src%first_item1)) then
      src%last_item1 => null()
   end if
   src%num_items1 = src%num_items1 - 1

   item_to_move%next_item => null()
   dest%parent_link%itemdir1(item_to_move%value)%ptr => dest

   if (.not. associated(dest%first_item1)) then
      dest%first_item1 => item_to_move
   else
      dest%last_item1%next_item => item_to_move
   end if
   dest%last_item1 => item_to_move
   dest%num_items1 = dest%num_items1 + 1
end subroutine

subroutine move_first_item2(src, dest)
   type(part_node_t), intent(inout) :: src
   type(part_node_t), target, intent(inout) :: dest
   type(item_node_t), pointer :: item_to_move

   item_to_move => src%first_item2
   if ((.not. associated(item_to_move))) error stop

   src%first_item2 => item_to_move%next_item
   if (.not. associated(src%first_item2)) then
      src%last_item2 => null()
   end if
   src%num_items2 = src%num_items2 - 1

   item_to_move%next_item => null()
   dest%parent_link%itemdir2(item_to_move%value)%ptr => dest

   if (.not. associated(dest%first_item2)) then
      dest%first_item2 => item_to_move
   else
      dest%last_item2%next_item => item_to_move
   end if
   dest%last_item2 => item_to_move
   dest%num_items2 = dest%num_items2 + 1
end subroutine

subroutine move_node_items(src, dest)
   type(part_node_t), intent(inout) :: src, dest

   do while (associated(src%first_item1))
      call move_first_item1(src, dest)
   end do

   do while (associated(src%first_item2))
      call move_first_item2(src, dest)
   end do
end subroutine

subroutine delete_chain(chain_root)
   type(chain_root_t), pointer, intent(inout) :: chain_root
   type(link_node_t), pointer :: link, next_link

   if ((.not. associated(chain_root))) error stop

   ! First delete all links
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
   type(partref_node_t), pointer :: partref, next_partref

   if ((.not. associated(parent_link))) error stop

   ! Delete all part references and only delete parts created by this partition
   partref => parent_link%first_partref
   do while (associated(partref))
      next_partref => partref%next_partref

      ! Only delete the part if this partition is its original creator
      if (associated(partref%part_node)) then
         if (associated(partref%part_node%parent_link, parent_link)) then
            ! This partition created this part, so we can safely delete it
            call delete_part(partref%part_node)
         end if
      end if

      ! Always delete the part reference
      deallocate(partref)
      partref => next_partref
   end do

   ! Deallocate directories (always allocated)
   deallocate(parent_link%itemdir1)
   deallocate(parent_link%itemdir2)

   ! Deallocate partition root
   deallocate(parent_link)
   parent_link => null()
end subroutine

subroutine delete_part(part_node)
   type(part_node_t), pointer, intent(inout) :: part_node

   if ((.not. associated(part_node))) error stop

   ! Deallocate items
   call deallocate_items(part_node%first_item1)
   call deallocate_items(part_node%first_item2)

   ! Deallocate signature array if allocated
   deallocate(part_node%signature)

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
      if (.not. associated(partref%part_node)) error stop 'Unexpected null part'

      partition%parts(i)%num_items1 = partref%part_node%num_items1
      partition%parts(i)%num_items2 = partref%part_node%num_items2
      allocate(partition%parts(i)%items1(partref%part_node%num_items1))
      allocate(partition%parts(i)%items2(partref%part_node%num_items2))

      item => partref%part_node%first_item1
      do j = 1, partref%part_node%num_items1
         partition%parts(i)%items1(j) = item%value
         partition%itemdir1(item%value) = i
         item => item%next_item
      end do

      item => partref%part_node%first_item2
      do j = 1, partref%part_node%num_items2
         partition%parts(i)%items2(j) = item%value
         partition%itemdir2(item%value) = i
         item => item%next_item
      end do

      partref => partref%next_partref
      i = i + 1
   end do
end subroutine

subroutine print_items(part)
   type(part_node_t), intent(in) :: part
   type(item_node_t), pointer :: item

   write(stderr, '(Z4.4,A)', advance='no') address(part), ':'

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

subroutine print_itemdirs(parent_link)
   type(link_node_t), intent(in) :: parent_link
   integer :: i

   write(stderr, *)
   write(stderr,'(A)') "Item Directory 1:"
   do i = 1, size(parent_link%itemdir1)
      if (associated(parent_link%itemdir1(i)%ptr)) then
         write(stderr,'(2X,I0,A,Z4.4)') i, " -> ", address(parent_link%itemdir1(i)%ptr)
      else
         write(stderr,'(2X,I0,A)') i, " -> Not associated"
      end if
   end do

   write(stderr, *)
   write(stderr,'(A)') "Item Directory 2:"
   do i = 1, size(parent_link%itemdir2)
      if (associated(parent_link%itemdir2(i)%ptr)) then
         write(stderr,'(2X,I0,A,Z4.4)') i, " -> ", address(parent_link%itemdir2(i)%ptr)
      else
         write(stderr,'(2X,I0,A)') i, " -> Not associated"
      end if
   end do
end subroutine

function init_chain(tot_items1, tot_items2) result(chain_root)
   integer, intent(in) :: tot_items1, tot_items2
   type(chain_root_t), pointer :: chain_root

   allocate(chain_root)
   chain_root%num_links = 0
   chain_root%tot_items1 = tot_items1
   chain_root%tot_items2 = tot_items2
   chain_root%first_link => null()
   chain_root%last_link => null()
end function

function add_new_link(chain_root) result(new_link)
   type(chain_root_t), target, intent(inout) :: chain_root
   type(link_node_t), pointer :: new_link
   integer :: i

   allocate(new_link)
   new_link%num_parts = 0
   new_link%first_partref => null()
   new_link%last_partref => null()
   new_link%next_link => null()
   new_link%chain_root => chain_root

   ! Allocate item directories
   allocate(new_link%itemdir1(chain_root%tot_items1))
   allocate(new_link%itemdir2(chain_root%tot_items2))

   ! Initialize all pointers to null
   do i = 1, chain_root%tot_items1
      new_link%itemdir1(i)%ptr => null()
   end do
   do i = 1, chain_root%tot_items2
      new_link%itemdir2(i)%ptr => null()
   end do

   ! Add to link to chain
   if (.not. associated(chain_root%first_link)) then
      chain_root%first_link => new_link
   else
      chain_root%last_link%next_link => new_link
   end if
   chain_root%last_link => new_link
   chain_root%num_links = chain_root%num_links + 1
end function

function add_new_link_branch(branch) result(new_link)
   type(branch_node_t), target, intent(inout) :: branch
   type(link_node_t), pointer :: new_link
   integer :: i

   allocate(new_link)
   new_link%num_parts = 0
   new_link%first_partref => null()
   new_link%last_partref => null()
   new_link%next_link => null()
!   new_link%branch => branch

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

   ! Add to link to chain
   if (.not. associated(branch%first_link)) then
      branch%first_link => new_link
   else
      branch%last_link%next_link => new_link
   end if
   branch%last_link => new_link
!   branch%num_links = branch%num_links + 1
end function

function find_child_part_node(part, signature) result(child_part)
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

   ! Print root with item counts and leaf indicator
   if (root_part%num_children == 0) then
      write(stderr, '(Z4.4,A,I0,A,I0,A)') address(root_part), ' (', &
         root_part%num_items1, '/', root_part%num_items2, ') *'
   else
      write(stderr, '(Z4.4,A,I0,A,I0,A)') address(root_part), ' (', &
         root_part%num_items1, '/', root_part%num_items2, ')'
   end if

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
      prefix = ""
      pos = 1
      do i = 1, depth
         if (is_last_child(i)) then
            prefix(pos:pos+3) = "    "
         else
            prefix(pos:pos+3) = "|   "
         end if
         pos = pos + 4
      end do

      ! Add branch characters
      if (is_last_child(depth + 1)) then
         prefix(pos:pos+3) = "`-- "
      else
         prefix(pos:pos+3) = "|-- "
      end if
      pos = pos + 4

      ! Print the child address with item counts and leaf indicator
      if (child_part%num_children == 0) then
         write(stderr, '(A,Z4.4,A,I0,A,I0,A)') prefix(1:pos-1), address(child_part), &
            ' (', child_part%num_items1, '/', child_part%num_items2, ') *'
      else
         write(stderr, '(A,Z4.4,A,I0,A,I0,A)') prefix(1:pos-1), address(child_part), &
            ' (', child_part%num_items1, '/', child_part%num_items2, ')'
      end if

      ! Recursively print this child's children
      call print_part_recursive(child_part, depth + 1, is_last_child)

      child_part => next_child
   end do
end subroutine

subroutine print_chain(chain_root)
   type(chain_root_t), pointer, intent(in) :: chain_root
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
   write(stderr,'(A,I0,A)') "( Level ", link_idx, " )"

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
      call print_items(partref%part_node)
      partref => partref%next_partref
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
      write(stderr,'(Z4.4,A)',advance='no') address(partref%part_node), ":"

      ! Print part children
      child_part => partref%part_node%first_child_part
      do while (associated(child_part))
         ! Print which part in next partition this child_part points to
         write(stderr,'(1X,Z4.4)',advance='no') address(child_part)
         child_part => child_part%next_sibling_part
      end do
      write(stderr,*)

      partref => partref%next_partref
   end do
end subroutine

subroutine print_partition_signatures(parent_link)
   type(link_node_t), pointer, intent(in) :: parent_link
   type(partref_node_t), pointer :: partref
   integer :: i

   if (.not. associated(parent_link)) error stop

   write(stderr,*)
   write(stderr,'(A)') "Signatures"
   write(stderr,'(A)') repeat("-", 11)

   partref => parent_link%first_partref
   do while (associated(partref))
      ! Print part address
      write(stderr,'(Z4.4,A)',advance='no') address(partref%part_node), ":"

      ! Print part signature
      do i = 1, size(partref%part_node%signature)
         write(stderr,'(*(1X,Z4.4))',advance='no') address(partref%part_node%signature(i)%ptr)
      end do
      write(stderr,*)

      partref => partref%next_partref
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
   write(stderr,'(A,I0,A)') "( Level ", link_idx, " )"

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

subroutine chain_from_partitionarray(partition, chain_root, root_branch)
   type(partitionarray_t), intent(in) :: partition
   type(chain_root_t), pointer, intent(out) :: chain_root
   type(branch_node_t), pointer, intent(out) :: root_branch
   type(link_node_t), pointer :: first_link, second_link, branch_link
   type(part_node_t), pointer :: root_part, new_part
   integer :: i, j

   ! Create root chain
   chain_root => init_chain(size(partition%itemdir1), size(partition%itemdir2))

   ! Create first and second links
   first_link => add_new_link(chain_root)
   second_link => add_new_link(chain_root)

   root_part => init_part(0)
   call add_part(first_link, root_part)
   root_branch => make_new_branch(root_part, size(partition%itemdir1), size(partition%itemdir2))
   branch_link => add_new_link_branch(root_branch)

   ! Add all items to the root_part
   do i = 1, partition%num_parts
      do j = 1, partition%parts(i)%num_items1
         call add_new_item1(root_part, partition%parts(i)%items1(j))
      end do

      do j = 1, partition%parts(i)%num_items2
         call add_new_item2(root_part, partition%parts(i)%items2(j))
      end do
   end do

   ! Create parts as children of root_part and add them to the first link
   do i = 1, partition%num_parts
      ! Create new part as child of root_part
      new_part => add_new_part(second_link, root_part)
      call add_part(branch_link, new_part)

      ! Add items to the child part (items remain in root_part too)
      do j = 1, partition%parts(i)%num_items1
         call add_new_item1(new_part, partition%parts(i)%items1(j))
      end do

      do j = 1, partition%parts(i)%num_items2
         call add_new_item2(new_part, partition%parts(i)%items2(j))
      end do
   end do

   call update_branch_parts(root_branch, root_part)
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
   type(partref_node_t), pointer :: partref, next_partref, prev_partref
   logical :: swapped

   ! Bubble sort implementation for part reference linked list
   do
      swapped = .false.
      partref => parent_link%first_partref
      prev_partref => null()
      do while (associated(partref) .and. associated(partref%next_partref))
         next_partref => partref%next_partref
         ! Check if we need to swap (current has more items than next)
         if (partref%part_node%num_items1 > next_partref%part_node%num_items1) then
            swapped = .true.
            ! Perform the swap
            partref%next_partref => next_partref%next_partref
            next_partref%next_partref => partref
            if (associated(prev_partref)) then
               prev_partref%next_partref => next_partref
            else
               ! Update first_partref if we're swapping the first element
               parent_link%first_partref => next_partref
            end if
            ! Update last_partref if necessary
            if (.not. associated(partref%next_partref)) then
               parent_link%last_partref => partref
            end if
            ! Update prev_partref for next iteration
            prev_partref => next_partref
         else
            ! No swap needed, just advance
            prev_partref => partref
            partref => next_partref
         end if
      end do
      ! If no swaps occurred, the list is sorted
      if (.not. swapped) exit
   end do
end subroutine

function make_new_branch(part, tot_items1, tot_items2) result(new_branch)
   type(part_node_t), target, intent(in) :: part
   integer, intent(in) :: tot_items1, tot_items2
   type(branch_node_t), pointer :: new_branch

   allocate(new_branch)
   new_branch%first_partref => null()
   new_branch%last_partref => null()
   new_branch%first_link => null()
   new_branch%last_link => null()
   new_branch%tot_items1 = tot_items1
   new_branch%tot_items2 = tot_items2
   new_branch%num_parts = 0
   new_branch%parent_branch => null()
   new_branch%first_child_branch => null()
   new_branch%last_child_branch => null()
   new_branch%next_sibling_branch => null()
   new_branch%num_children = 0
   new_branch%part_node => part
end function

function add_new_branch(branch, part) result(new_branch)
   type(branch_node_t), pointer, intent(inout) :: branch
   type(part_node_t), target, intent(in) :: part
   type(branch_node_t), pointer :: new_branch

   ! Create new child branch
   new_branch => make_new_branch(part, branch%tot_items1, branch%tot_items2)

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

subroutine add_branch_part(branch, part)
! Adds parts in sorted order by num_items1
   type(branch_node_t), target, intent(inout) :: branch
   type(part_node_t), target, intent(in) :: part
   type(partref_node_t), pointer :: new_partref, current, prev

   ! Create new part reference
   allocate(new_partref)
   new_partref%part_node => part
   new_partref%next_partref => null()

   ! If branch is empty, add as first element
   if (.not. associated(branch%first_partref)) then
      branch%first_partref => new_partref
      branch%last_partref => new_partref
      branch%num_parts = branch%num_parts + 1
      return
   end if

   ! Find correct insertion position (sorted by increasing num_items1)
   current => branch%first_partref
   prev => null()

   do while (associated(current))
      ! If new part has fewer or equal items1, insert before current
      if (part%num_items1 <= current%part_node%num_items1) then
         exit
      end if
      prev => current
      current => current%next_partref
   end do

   ! Insert new_partref at the found position
   new_partref%next_partref => current

   if (associated(prev)) then
      ! Insert in middle or at end
      prev%next_partref => new_partref
      ! Update last_partref if we inserted at the end
      if (.not. associated(current)) then
         branch%last_partref => new_partref
      end if
   else
      ! Insert at beginning
      branch%first_partref => new_partref
      ! If this was the only element, also update last_partref
      if (.not. associated(current)) then
         branch%last_partref => new_partref
      end if
   end if

   branch%num_parts = branch%num_parts + 1
end subroutine

! New procedure to remove a specific part from branch
subroutine remove_branch_part(branch, part_to_remove)
   type(branch_node_t), target, intent(inout) :: branch
   type(part_node_t), target, intent(in) :: part_to_remove
   type(partref_node_t), pointer :: current, prev, to_delete

   if (.not. associated(branch%first_partref)) return

   ! Search for the part to remove
   current => branch%first_partref
   prev => null()

   do while (associated(current))
      if (associated(current%part_node, part_to_remove)) then
         ! Found the part to remove
         to_delete => current

         ! Update links to bypass the node to be deleted
         if (associated(prev)) then
            ! Removing from middle or end
            prev%next_partref => current%next_partref
            ! Update last_partref if we're removing the last element
            if (.not. associated(current%next_partref)) then
               branch%last_partref => prev
            end if
         else
            ! Removing first element
            branch%first_partref => current%next_partref
            ! Update last_partref if we're removing the only element
            if (.not. associated(current%next_partref)) then
               branch%last_partref => null()
            end if
         end if

         ! Deallocate the partref node
         deallocate(to_delete)
         branch%num_parts = branch%num_parts - 1
         return
      end if

      prev => current
      current => current%next_partref
   end do

   ! Part not found - no error, just return silently as requested
end subroutine

subroutine update_branch_parts(branch, part)
! Adds all children of a part to the branch
! Only adds children with more than 1 items1
   type(branch_node_t), target, intent(inout) :: branch
   type(part_node_t), target, intent(in) :: part
   type(part_node_t), pointer :: child_part

   ! Remove parent branch if it exists
   call remove_branch_part(branch, part)

   ! Traverse all children of the parent part
   child_part => part%first_child_part
   do while (associated(child_part))
      ! Only add children with more than 1 items1
      if (child_part%num_items1 > 1) then
         call add_branch_part(branch, child_part)
      end if
      child_part => child_part%next_sibling_part
   end do
end subroutine

subroutine print_branch_tree(root_branch)
   type(branch_node_t), pointer, intent(in) :: root_branch
   logical, dimension(:), allocatable :: is_last_child

   if (.not. associated(root_branch)) then
      write(stderr, '(A)') "Branch tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "       Branch Tree"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Allocate tracking array for tree lines (max depth 100)
   allocate(is_last_child(100))
   is_last_child = .false.

   ! Print root with part address, partref count and leaf indicator
   if (associated(root_branch%part_node)) then
      if (root_branch%num_children == 0) then
         write(stderr, '(Z4.4,A,I0,A)') address(root_branch%part_node), ' (', &
            root_branch%num_parts, ') *'
      else
         write(stderr, '(Z4.4,A,I0,A)') address(root_branch%part_node), ' (', &
            root_branch%num_parts, ')'
      end if
   else
      if (root_branch%num_children == 0) then
         write(stderr, '(A,I0,A)') 'NULL (', root_branch%num_parts, ') *'
      else
         write(stderr, '(A,I0,A)') 'NULL (', root_branch%num_parts, ')'
      end if
   end if

   ! Print children recursively
   call print_branch_recursive(root_branch, 0, is_last_child)

   deallocate(is_last_child)
   write(stderr, *)
end subroutine

recursive subroutine print_branch_recursive(branch, depth, is_last_child)
   type(branch_node_t), pointer, intent(in) :: branch
   integer, intent(in) :: depth
   logical, dimension(:), intent(inout) :: is_last_child
   type(branch_node_t), pointer :: child_branch, next_child
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
      prefix = ""
      pos = 1
      do i = 1, depth
         if (is_last_child(i)) then
            prefix(pos:pos+3) = "    "
         else
            prefix(pos:pos+3) = "|   "
         end if
         pos = pos + 4
      end do

      ! Add branch characters
      if (is_last_child(depth + 1)) then
         prefix(pos:pos+3) = "`-- "
      else
         prefix(pos:pos+3) = "|-- "
      end if
      pos = pos + 4

      ! Print the child address with partref count and leaf indicator
      if (associated(child_branch%part_node)) then
         if (child_branch%num_children == 0) then
            write(stderr, '(A,Z4.4,A,I0,A)') prefix(1:pos-1), address(child_branch%part_node), &
               ' (', child_branch%num_parts, ') *'
         else
            write(stderr, '(A,Z4.4,A,I0,A)') prefix(1:pos-1), address(child_branch%part_node), &
               ' (', child_branch%num_parts, ')'
         end if
      else
         if (child_branch%num_children == 0) then
            write(stderr, '(A,A,I0,A)') prefix(1:pos-1), 'NULL (', &
               child_branch%num_parts, ') *'
         else
            write(stderr, '(A,A,I0,A)') prefix(1:pos-1), 'NULL (', &
               child_branch%num_parts, ')'
         end if
      end if

      ! Recursively print this child's children
      call print_branch_recursive(child_branch, depth + 1, is_last_child)

      child_branch => next_child
   end do
end subroutine

subroutine print_branch_contents(root_branch)
   type(branch_node_t), pointer, intent(in) :: root_branch

   if (.not. associated(root_branch)) then
      write(stderr, '(A)') "Branch tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "    Branch Contents"
   write(stderr, '(A)') repeat("=", 25)

   call print_branch_contents_recursive(root_branch)
   write(stderr, *)
end subroutine

recursive subroutine print_branch_contents_recursive(branch)
   type(branch_node_t), pointer, intent(in) :: branch
   type(branch_node_t), pointer :: child_branch
   type(partref_node_t), pointer :: partref

   if (.not. associated(branch)) return

   ! Print branch node information
   write(stderr, '(A,Z4.4,A,I0,A)') "Branch ", address(branch%part_node), " (", branch%num_parts, " parts)"

   ! Print all parts in this branch
   partref => branch%first_partref
   do while (associated(partref))
      write(stderr, '(A,Z4.4,A,I0,A,I0,A)') "  Part ", &
         address(partref%part_node), " (", partref%part_node%num_items1, &
         "/", partref%part_node%num_items2, " items)"
      partref => partref%next_partref
   end do

   ! Print all child branches
   child_branch => branch%first_child_branch
   do while (associated(child_branch))
      call print_branch_contents_recursive(child_branch)
      child_branch => child_branch%next_sibling_branch
   end do
end subroutine

end module
