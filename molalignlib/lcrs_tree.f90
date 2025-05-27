module lcrs_tree
use parameters
implicit none
private

! Item node
type, public :: item_node_t
   integer :: value
   type(item_node_t), pointer :: next_item
end type

! Chain root (simplified from chain_root_t)
type, public :: chain_root_t
   integer :: num_links
   integer :: tot_items1
   integer :: tot_items2
   type(link_node_t), pointer :: first_link
   type(link_node_t), pointer :: last_link
end type

! Part reference node - new indirection layer
type, public :: partref_node_t
   type(part_node_t), pointer :: part_ptr
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
   integer :: index
   integer :: num_items1
   integer :: num_items2
   integer :: num_children
   integer :: num_leaves
   type(item_node_t), pointer :: first_item1
   type(item_node_t), pointer :: first_item2
   type(item_node_t), pointer :: last_item1
   type(item_node_t), pointer :: last_item2
   type(link_node_t), pointer :: partition_root
   type(part_node_t), pointer :: parent
   type(part_node_t), pointer :: first_child
   type(part_node_t), pointer :: last_child
   type(part_node_t), pointer :: next_sibling
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

! Make types and procedures public
public address
public make_chain_root
public make_new_part
public add_new_link
public add_part
public add_new_part
public add_new_item1
public add_new_item2
public move_first_item1
public move_first_item2
public move_node_items
public delete_chain
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

subroutine update_leaf_counts_on_new_leaf(new_leaf)
   type(part_node_t), pointer, intent(inout) :: new_leaf
   type(part_node_t), pointer :: ancestor

   ! When a new leaf is added, increment all ancestors' leaf counts
   ancestor => new_leaf%parent
   do while (associated(ancestor))
      ancestor%num_leaves = ancestor%num_leaves + 1
      ancestor => ancestor%parent
   end do
end subroutine

function make_new_part() result(new_part)
   type(part_node_t), pointer :: new_part

   allocate (new_part)

   new_part%index = 0
   new_part%num_children = 0
   new_part%num_items1 = 0
   new_part%num_items2 = 0
   new_part%num_leaves = 1
   new_part%parent => null()
   new_part%first_item1 => null()
   new_part%last_item1 => null()
   new_part%first_item2 => null()
   new_part%last_item2 => null()
   new_part%first_child => null()
   new_part%last_child => null()
   new_part%next_sibling => null()
   new_part%partition_root => null()
   new_part%signature = [part_nodeptr_t::]
end function

subroutine add_part(link, part)
   type(link_node_t), target, intent(inout) :: link
   type(part_node_t), target, intent(inout) :: part
   type(partref_node_t), pointer :: new_partref
   type(item_node_t), pointer :: item

   part%partition_root => link

   allocate(new_partref)
   new_partref%part_ptr => part
   new_partref%next_partref => null()

   if (.not. associated(link%first_partref)) then
      link%first_partref => new_partref
   else
      link%last_partref%next_partref => new_partref
   end if
   link%last_partref => new_partref
   link%num_parts = link%num_parts + 1

   if (part%index == 0) then
      part%index = link%num_parts
   end if

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

function add_new_part(partition_root, parent_part) result(new_part)
   type(link_node_t), pointer, intent(inout) :: partition_root
   type(part_node_t), pointer, intent(inout) :: parent_part
   type(part_node_t), pointer :: new_part

   ! Create new part but don't add it yet
   new_part => make_new_part()

   ! Set parent BEFORE calling add_part
   new_part%parent => parent_part
   new_part%next_sibling => null()

   ! Now add to partition
   call add_part(partition_root, new_part)

   ! Set up parent-child relationships
   if (.not. associated(parent_part%first_child)) then
      parent_part%first_child => new_part
   else
      parent_part%last_child%next_sibling => new_part
   end if
   parent_part%last_child => new_part
   parent_part%num_children = parent_part%num_children + 1

   ! Update leaf counts
   if (parent_part%num_children >= 2) then
      call update_leaf_counts_on_new_leaf(new_part)
   end if
end function

subroutine add_new_item1(part, value)
   type(part_node_t), target, intent(inout) :: part
   integer, intent(in) :: value
   type(item_node_t), pointer :: new_item

   allocate(new_item)
   new_item%value = value
   new_item%next_item => null()
   part%partition_root%itemdir1(value)%ptr => part

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
   part%partition_root%itemdir2(value)%ptr => part

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
   dest%partition_root%itemdir1(item_to_move%value)%ptr => dest

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
   dest%partition_root%itemdir2(item_to_move%value)%ptr => dest

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

subroutine delete_partition(partition_root)
   type(link_node_t), pointer, intent(inout) :: partition_root
   type(partref_node_t), pointer :: partref, next_partref

   if ((.not. associated(partition_root))) error stop

   ! Delete all part references and only delete parts created by this partition
   partref => partition_root%first_partref
   do while (associated(partref))
      next_partref => partref%next_partref
      
      ! Only delete the part if this partition is its original creator
      if (associated(partref%part_ptr)) then
         if (associated(partref%part_ptr%partition_root, partition_root)) then
            ! This partition created this part, so we can safely delete it
            call delete_part(partref%part_ptr)
         end if
      end if
      
      ! Always delete the part reference
      deallocate(partref)
      partref => next_partref
   end do

   ! Deallocate directories (always allocated)
   deallocate(partition_root%itemdir1)
   deallocate(partition_root%itemdir2)

   ! Deallocate partition root
   deallocate(partition_root)
   partition_root => null()
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

subroutine partition_to_partitionarray(partition_root, partition)
   type(link_node_t), target, intent(in) :: partition_root
   type(partitionarray_t), intent(out) :: partition
   type(partref_node_t), pointer :: partref
   type(item_node_t), pointer :: item
   integer :: i, j

   partition%num_parts = partition_root%num_parts
   allocate(partition%parts(partition%num_parts))
   allocate(partition%itemdir1(size(partition_root%itemdir1)))
   allocate(partition%itemdir2(size(partition_root%itemdir2)))

   partref => partition_root%first_partref
   i = 1
   do while (associated(partref))
      if (.not. associated(partref%part_ptr)) error stop 'Unexpected null part'

      partition%parts(i)%num_items1 = partref%part_ptr%num_items1
      partition%parts(i)%num_items2 = partref%part_ptr%num_items2
      allocate(partition%parts(i)%items1(partref%part_ptr%num_items1))
      allocate(partition%parts(i)%items2(partref%part_ptr%num_items2))

      item => partref%part_ptr%first_item1
      do j = 1, partref%part_ptr%num_items1
         partition%parts(i)%items1(j) = item%value
         partition%itemdir1(item%value) = i
         item => item%next_item
      end do

      item => partref%part_ptr%first_item2
      do j = 1, partref%part_ptr%num_items2
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

subroutine print_itemdirs(partition_root)
   type(link_node_t), intent(in) :: partition_root
   integer :: i

   write(stderr, *)
   write(stderr,'(A)') "Item Directory 1:"
   do i = 1, size(partition_root%itemdir1)
      if (associated(partition_root%itemdir1(i)%ptr)) then
         write(stderr,'(2X,I0,A,Z4.4)') i, " -> ", address(partition_root%itemdir1(i)%ptr)
      else
         write(stderr,'(2X,I0,A)') i, " -> Not associated"
      end if
   end do

   write(stderr, *)
   write(stderr,'(A)') "Item Directory 2:"
   do i = 1, size(partition_root%itemdir2)
      if (associated(partition_root%itemdir2(i)%ptr)) then
         write(stderr,'(2X,I0,A,Z4.4)') i, " -> ", address(partition_root%itemdir2(i)%ptr)
      else
         write(stderr,'(2X,I0,A)') i, " -> Not associated"
      end if
   end do
end subroutine

function make_chain_root(tot_items1, tot_items2) result(chain_root)
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

function find_child_part_node(parent_part, signature) result(child_part)
   type(part_node_t), intent(in) :: parent_part
   type(part_nodeptr_t), dimension(:), intent(in) :: signature
   type(part_node_t), pointer :: child_part

   child_part => parent_part%first_child
   do while (associated(child_part))
      if (child_part%signature .equiv. signature) then
         return
      end if
      child_part => child_part%next_sibling
   end do

   child_part => null()
end function

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

subroutine print_partition(partition_root, link_idx)
   type(link_node_t), pointer, intent(in) :: partition_root
   integer, intent(in) :: link_idx

   if (.not. associated(partition_root)) error stop

   write(stderr,*)
   write(stderr,'(A,I0,A)') "( Level ", link_idx, " )"

   ! Print part items
   call print_partition_items(partition_root)

   ! Print part signatures
   call print_partition_signatures(partition_root)

   ! Print part children
   call print_partition_children(partition_root)
end subroutine

subroutine print_partition_items(partition_root)
   type(link_node_t), pointer, intent(in) :: partition_root
   type(partref_node_t), pointer :: partref

   if (.not. associated(partition_root)) error stop

   write(stderr, *)
   write(stderr,'(A)') "Items"
   write(stderr,'(A)') repeat("-", 6)

   partref => partition_root%first_partref
   do while (associated(partref))
      call print_items(partref%part_ptr)
      partref => partref%next_partref
   end do
end subroutine

subroutine print_partition_children(partition_root)
   type(link_node_t), pointer, intent(in) :: partition_root
   type(partref_node_t), pointer :: partref
   type(part_node_t), pointer :: child_part

   if (.not. associated(partition_root)) error stop

   write(stderr,*)
   write(stderr,'(A)') "Children"
   write(stderr,'(A)') repeat("-", 9)

   partref => partition_root%first_partref
   do while (associated(partref))
      ! Print part address
      write(stderr,'(Z4.4,A)',advance='no') address(partref%part_ptr), ":"

      ! Print part children
      child_part => partref%part_ptr%first_child
      do while (associated(child_part))
         ! Print which part in next partition this child_part points to
         write(stderr,'(1X,Z4.4)',advance='no') address(child_part)
         child_part => child_part%next_sibling
      end do
      write(stderr,*)

      partref => partref%next_partref
   end do
end subroutine

subroutine print_partition_signatures(partition_root)
   type(link_node_t), pointer, intent(in) :: partition_root
   type(partref_node_t), pointer :: partref
   integer :: i

   if (.not. associated(partition_root)) error stop

   write(stderr,*)
   write(stderr,'(A)') "Signatures"
   write(stderr,'(A)') repeat("-", 11)

   partref => partition_root%first_partref
   do while (associated(partref))
      ! Print part address
      write(stderr,'(Z4.4,A)',advance='no') address(partref%part_ptr), ":"

      ! Print part signature
      do i = 1, size(partref%part_ptr%signature)
         write(stderr,'(*(1X,Z4.4))',advance='no') address(partref%part_ptr%signature(i)%ptr)
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

function chain_from_partitionarray(partition) result(chain_root)
   type(partitionarray_t), intent(in) :: partition
   type(chain_root_t), pointer :: chain_root
   type(link_node_t), pointer :: new_link
   type(part_node_t), pointer :: part
   integer :: i, j

   ! Create root chain
   chain_root => make_chain_root(size(partition%itemdir1), size(partition%itemdir2))
   new_link => add_new_link(chain_root)

   ! Create parts and add items
   do i = 1, partition%num_parts
      part => make_new_part()
      call add_part(new_link, part)

      do j = 1, partition%parts(i)%num_items1
         call add_new_item1(part, partition%parts(i)%items1(j))
      end do

      do j = 1, partition%parts(i)%num_items2
         call add_new_item2(part, partition%parts(i)%items2(j))
      end do
   end do
end function

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

subroutine sort_parts_by_size(partition_root)
   type(link_node_t), pointer, intent(inout) :: partition_root
   type(partref_node_t), pointer :: partref, next_partref, prev_partref
   logical :: swapped

   ! Bubble sort implementation for part reference linked list
   do
      swapped = .false.
      partref => partition_root%first_partref
      prev_partref => null()
      do while (associated(partref) .and. associated(partref%next_partref))
         next_partref => partref%next_partref
         ! Check if we need to swap (current has more items than next)
         if (partref%part_ptr%num_items1 > next_partref%part_ptr%num_items1) then
            swapped = .true.
            ! Perform the swap
            partref%next_partref => next_partref%next_partref
            next_partref%next_partref => partref
            if (associated(prev_partref)) then
               prev_partref%next_partref => next_partref
            else
               ! Update first_partref if we're swapping the first element
               partition_root%first_partref => next_partref
            end if
            ! Update last_partref if necessary
            if (.not. associated(partref%next_partref)) then
               partition_root%last_partref => partref
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

end module
