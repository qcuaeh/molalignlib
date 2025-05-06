module lcrs_tree
use parameters
implicit none
private

! Item node
type, public :: item_node_t
   integer :: value
   type(item_node_t), pointer :: next_item
end type

! Tree node
type, public :: tree_node_t
   integer :: num_links
   integer :: tot_items1
   integer :: tot_items2
   type(tree_node_t), pointer :: first_sub_branch
   type(tree_node_t), pointer :: last_sub_branch
   type(tree_node_t), pointer :: next_branch
   type(link_node_t), pointer :: first_link
   type(link_node_t), pointer :: last_link
end type

! Link node
type, public :: link_node_t
   integer :: num_parts
   type(link_node_t), pointer :: next_link
   type(tree_node_t), pointer :: chain_root
   type(part_node_t), pointer :: first_part
   type(part_node_t), pointer :: last_part
   type(part_nodeptr_t), dimension(:), allocatable :: itemdir1
   type(part_nodeptr_t), dimension(:), allocatable :: itemdir2
end type

! Part node
type, public :: part_node_t
   integer :: index
   integer :: num_items1
   integer :: num_items2
   integer :: num_children
   integer :: num_neighbors
   integer :: num_leaves
   type(item_node_t), pointer :: first_item1
   type(item_node_t), pointer :: first_item2
   type(item_node_t), pointer :: last_item1
   type(item_node_t), pointer :: last_item2
   type(link_node_t), pointer :: partition_root
   type(part_node_t), pointer :: next_part
   type(part_node_t), pointer :: parent
   type(part_node_t), pointer :: first_child
   type(part_node_t), pointer :: last_child
   type(part_node_t), pointer :: next_sibling
   type(part_nodeptr_t), dimension(:), allocatable :: neighbors
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
   integer :: num_neighbors
   integer :: num_items1
   integer :: num_items2
   integer, dimension(:), allocatable :: children
   integer, dimension(:), allocatable :: items1
   integer, dimension(:), allocatable :: items2
   integer, dimension(:), allocatable :: neighbors
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
   module procedure neighbor_nodes_equivalence
end interface

interface add_new_part
   module procedure add_new_part_root
   module procedure add_new_part_root_parent
   module procedure add_new_part_root_parent_signature
end interface

! Make types and procedures public
public address
public make_new_tree
public add_new_branch
public add_new_link
public add_new_part
public add_new_item1
public add_new_item2
public move_first_item1
public move_first_item2
public move_node_items
public delete_tree
public print_tree
public print_chain
public print_chainarray
public print_items
public print_itemdirs
public sort_parts_by_size
public partition_to_partitionarray
public tree_from_partitionarray
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

function neighbor_nodes_equivalence(array1, array2) result(equiv)
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

function add_new_part_root(partition_root) result(new_part)
   type(link_node_t), target, intent(inout) :: partition_root
   type(part_node_t), pointer :: new_part

   allocate (new_part)

   new_part%num_neighbors = 0
   new_part%num_children = 0
   new_part%num_items1 = 0
   new_part%num_items2 = 0
   new_part%num_leaves = 1
   new_part%first_item1 => null()
   new_part%last_item1 => null()
   new_part%first_item2 => null()
   new_part%last_item2 => null()
   new_part%first_child => null()
   new_part%last_child => null()
   new_part%next_sibling => null()
   new_part%partition_root => partition_root
   new_part%parent => null()
   new_part%next_part => null()

   if (.not. associated(partition_root%first_part)) then
      partition_root%first_part => new_part
   else
      partition_root%last_part%next_part => new_part
   end if
   partition_root%last_part => new_part
   partition_root%num_parts = partition_root%num_parts + 1
   new_part%index = partition_root%num_parts
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

function add_new_part_root_parent(partition_root, parent_part) result(new_part)
   type(link_node_t), pointer, intent(inout) :: partition_root
   type(part_node_t), pointer, intent(inout) :: parent_part
   type(part_node_t), pointer :: new_part

   new_part => add_new_part_root(partition_root)
   new_part%parent => parent_part
   new_part%next_sibling => null()

   if (.not. associated(parent_part%first_child)) then
      parent_part%first_child => new_part
   else
      parent_part%last_child%next_sibling => new_part
   end if
   parent_part%last_child => new_part
   parent_part%num_children = parent_part%num_children + 1

   ! Only update leaf counts starting from the second child
   ! First child: parent was leaf (1), now has 1 child (1) - no net change needed
   ! Second+ child: each additional child adds 1 to the leaf count
   if (parent_part%num_children >= 2) then
      call update_leaf_counts_on_new_leaf(new_part)
   end if
end function

function add_new_part_root_parent_signature(partition_root, parent_part, neighbors) result(new_part)
   type(link_node_t), pointer, intent(inout) :: partition_root
   type(part_node_t), pointer, intent(inout) :: parent_part
   type(part_nodeptr_t), dimension(:), intent(in):: neighbors
   type(part_node_t), pointer :: new_part

   new_part => add_new_part_root_parent(partition_root, parent_part)
   new_part%neighbors = neighbors
   new_part%num_neighbors = size(neighbors)
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

subroutine delete_tree(tree_root)
   type(tree_node_t), pointer, intent(inout) :: tree_root

   if ((.not. associated(tree_root))) error stop

   call delete_tree_hierarchy(tree_root)
   tree_root => null()
end subroutine

recursive subroutine delete_tree_hierarchy(tree_node)
   type(tree_node_t), pointer, intent(inout) :: tree_node
   type(tree_node_t), pointer :: branch, next_branch

   if ((.not. associated(tree_node))) error stop

   ! First delete all branches
   branch => tree_node%first_sub_branch
   do while (associated(branch))
      next_branch => branch%next_branch
      call delete_tree_hierarchy(branch)
      branch => next_branch
   end do

   ! Then delete this node
   call delete_chain(tree_node)
end subroutine

subroutine delete_chain(chain_root)
   type(tree_node_t), pointer, intent(inout) :: chain_root
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
   type(part_node_t), pointer :: part, next_part

   if ((.not. associated(partition_root))) error stop

   ! Delete all parts
   part => partition_root%first_part
   do while (associated(part))
      next_part => part%next_part
      call delete_part(part)
      part => next_part
   end do

   ! Deallocate directories
   deallocate(partition_root%itemdir1)
   deallocate(partition_root%itemdir2)

   ! Deallocate partition root
   deallocate(partition_root)
   partition_root => null()
end subroutine

subroutine delete_part(part_root)
   type(part_node_t), pointer, intent(inout) :: part_root

   if ((.not. associated(part_root))) error stop

   call deallocate_items(part_root%first_item1)
   call deallocate_items(part_root%first_item2)

   deallocate(part_root)
   part_root => null()
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
   type(part_node_t), pointer :: part
   type(item_node_t), pointer :: item
   integer :: i, j

   partition%num_parts = partition_root%num_parts
   allocate(partition%parts(partition%num_parts))
   allocate(partition%itemdir1(size(partition_root%itemdir1)))
   allocate(partition%itemdir2(size(partition_root%itemdir2)))

   part => partition_root%first_part
   do i = 1, partition%num_parts
      if (.not. associated(part)) error stop 'Unexpected null part'
      if (.not. associated(part%partition_root, partition_root)) error stop 'Leaf points to wrong partition root'

      partition%parts(i)%num_neighbors = part%num_neighbors
      partition%parts(i)%num_items1 = part%num_items1
      partition%parts(i)%num_items2 = part%num_items2
      allocate(partition%parts(i)%items1(part%num_items1))
      allocate(partition%parts(i)%items2(part%num_items2))

      item => part%first_item1
      do j = 1, part%num_items1
         partition%parts(i)%items1(j) = item%value
         partition%itemdir1(item%value) = i
         item => item%next_item
      end do

      item => part%first_item2
      do j = 1, part%num_items2
         partition%parts(i)%items2(j) = item%value
         partition%itemdir2(item%value) = i
         item => item%next_item
      end do

      part => part%next_part
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

function make_new_tree(tot_items1, tot_items2) result(tree_root)
   integer, intent(in) :: tot_items1, tot_items2
   type(tree_node_t), pointer :: tree_root

   allocate(tree_root)
   tree_root%num_links = 0
   tree_root%tot_items1 = tot_items1
   tree_root%tot_items2 = tot_items2
   tree_root%first_link => null()
   tree_root%last_sub_branch => null()
   tree_root%next_branch => null()
   tree_root%first_sub_branch => null()
   tree_root%last_link => null()
end function

function add_new_link(chain_root) result(new_link)
   type(tree_node_t), target, intent(inout) :: chain_root
   type(link_node_t), pointer :: new_link
   integer :: i

   allocate(new_link)
   new_link%num_parts = 0
   new_link%first_part => null()
   new_link%last_part => null()
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

function add_new_branch(tree_node) result(branch)
   type(tree_node_t), target, intent(inout) :: tree_node
   type(tree_node_t), pointer :: branch

   allocate(branch)
   branch%num_links = 0
   branch%tot_items1 = tree_node%tot_items1
   branch%tot_items2 = tree_node%tot_items2
   branch%first_link => null()
   branch%last_link => null()
   branch%first_sub_branch => null()
   branch%last_sub_branch => null()
   branch%next_branch => null()

   ! Add branch to tree node
   if (.not. associated(tree_node%first_sub_branch)) then
      ! First branch
      tree_node%first_sub_branch => branch
      tree_node%last_sub_branch => branch
   else
      ! Add to end of the linked list
      tree_node%last_sub_branch%next_branch => branch
      tree_node%last_sub_branch => branch
   end if
end function

function find_child_part_node(parent_part, neighbors) result(child_part)
   type(part_node_t), intent(in) :: parent_part
   type(part_nodeptr_t), dimension(:), intent(in) :: neighbors
   type(part_node_t), pointer :: child_part

   child_part => parent_part%first_child
   do while (associated(child_part))
      if (child_part%neighbors .equiv. neighbors) then
         return
      end if
      child_part => child_part%next_sibling
   end do

   child_part => null()
end function

subroutine print_tree(tree_root)
   type(tree_node_t), pointer, intent(in) :: tree_root

   if (.not. associated(tree_root)) error stop

   write(stderr,*)
   write(stderr,'(A)') repeat("=", 30)
   write(stderr,'(A)') " Tree hierarchy (depth first)"
   write(stderr,'(A)') repeat("=", 30)

   ! Start recursive printing from the root branch
   call print_tree_hierarchy(tree_root, [1])

   write(stderr,*)  ! Final newline
end subroutine

recursive subroutine print_tree_hierarchy(tree_node, path)
   type(tree_node_t), pointer, intent(in) :: tree_node
   integer, dimension(:), intent(in) :: path
   type(tree_node_t), pointer :: branch
   integer :: i
   integer, dimension(:), allocatable :: current_path

   if (.not. associated(tree_node)) error stop

   ! Print this tree node (chain and all its partitions)
   call print_chain(tree_node, path)

   ! Recursively process branches in depth-first order
   branch => tree_node%first_sub_branch
   i = 1
   do while (associated(branch))
      ! Create new path array for branch
      allocate(current_path(size(path) + 1))
      current_path(1:size(path)) = path
      current_path(size(path)+1) = i

      ! Recursive call for branch with updated path
      call print_tree_hierarchy(branch, current_path)

      deallocate(current_path)

      ! Move to next sibling
      branch => branch%next_branch
      i = i + 1
   end do
end subroutine

subroutine print_chain(chain_root, path)
   type(tree_node_t), pointer, intent(in) :: chain_root
   integer, dimension(:), intent(in) :: path
   type(link_node_t), pointer :: link
   integer :: i, link_idx
   character(len=200) :: path_str

   if (.not. associated(chain_root)) error stop

   ! Create current path string representation
   path_str = ''
   do i = 1, size(path)
      if (i > 1) path_str = trim(path_str) // '.'
      write(path_str(len_trim(path_str)+1:), '(I0)') path(i)
   end do

   write(stderr,*)
   write(stderr,'(A)') repeat("-", len_trim(path_str)+7)
   write(stderr,'(A)') " Node " // trim(path_str)
   write(stderr,'(A)') repeat("-", len_trim(path_str)+7)

   ! Print all partitions in this chain
   link => chain_root%first_link
   link_idx = 1
   do while (associated(link))
      call print_partition(link, link_idx)
      link => link%next_link
      link_idx = link_idx + 1
   end do
end subroutine

subroutine print_partition(partition_root, link_idx)
   type(link_node_t), pointer, intent(in) :: partition_root
   integer, intent(in) :: link_idx

   if (.not. associated(partition_root)) error stop

   write(stderr,*)
   write(stderr,'(A,I0,A)') "( Level ", link_idx, " )"

   ! Print part neighbors
   call print_partition_neighbors(partition_root)

   ! Print part items
   call print_partition_items(partition_root)

   ! Print part children
   call print_partition_children(partition_root)
end subroutine

subroutine print_partition_items(partition_root)
   type(link_node_t), pointer, intent(in) :: partition_root
   type(part_node_t), pointer :: part

   if (.not. associated(partition_root)) error stop

   write(stderr, *)
   write(stderr,'(A)') "Items"
   write(stderr,'(A)') repeat("-", 6)

   part => partition_root%first_part
   do while (associated(part))
      call print_items(part)
      part => part%next_part
   end do
end subroutine

subroutine print_partition_children(partition_root)
   type(link_node_t), pointer, intent(in) :: partition_root
   type(part_node_t), pointer :: part, child_part

   if (.not. associated(partition_root)) error stop

   write(stderr,*)
   write(stderr,'(A)') "Children"
   write(stderr,'(A)') repeat("-", 9)

   part => partition_root%first_part
   do while (associated(part))
      ! Print part address
      write(stderr,'(Z4.4,A)',advance='no') address(part), ":"

      ! Print part children
      child_part => part%first_child
      do while (associated(child_part))
         ! Print which part in next partition this child_part points to
         write(stderr,'(1X,Z4.4)',advance='no') address(child_part)
         child_part => child_part%next_sibling
      end do
      write(stderr,*)

      part => part%next_part
   end do
end subroutine

subroutine print_partition_neighbors(partition_root)
   type(link_node_t), pointer, intent(in) :: partition_root
   type(part_node_t), pointer :: part
   integer :: i

   if (.not. associated(partition_root)) error stop

   write(stderr,*)
   write(stderr,'(A)') "Neighbors"
   write(stderr,'(A)') repeat("-", 10)

   part => partition_root%first_part
   do while (associated(part))
      ! Print part address
      write(stderr,'(Z4.4,A)',advance='no') address(part), ":"

      ! Print part neighbors
      do i = 1, part%num_neighbors
         write(stderr,'(*(1X,Z4.4))',advance='no') address(part%neighbors(i)%ptr)
      end do
      write(stderr,*)

      part => part%next_part
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

   ! Print part neighbors
   call print_partitionarray_neighbors(chainarray, link_idx)

   ! Print part items
   call print_partitionarray_items(chainarray, link_idx)

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

subroutine print_partitionarray_neighbors(chainarray, link_idx)
   type(chainarray_t), intent(in) :: chainarray
   integer, intent(in) :: link_idx
   integer :: i, j

   write(stderr,*)
   write(stderr,'(A)') "Neighbors"
   write(stderr,'(A)') repeat("-", 10)

   do i = 1, chainarray%partitions(link_idx)%num_parts
      write(stderr,'(I3,A)',advance='no') i, ":"
      do j = 1, chainarray%partitions(link_idx)%parts(i)%num_neighbors
         write(stderr,'(1X,I3)',advance='no') &
            chainarray%partitions(link_idx)%parts(i)%neighbors(j)
      end do
      write(stderr,*)
   end do
end subroutine

function tree_from_partitionarray(partition) result(tree_root)
   type(partitionarray_t), intent(in) :: partition
   type(tree_node_t), pointer :: tree_root
   type(link_node_t), pointer :: new_link
   type(part_node_t), pointer :: part
   integer :: i, j

   ! Create root branch
   tree_root => make_new_tree(size(partition%itemdir1), size(partition%itemdir2))
   new_link => add_new_link(tree_root)

   ! Create parts and add items
   do i = 1, partition%num_parts
      part => add_new_part(new_link)
      part%num_neighbors = partition%parts(i)%num_neighbors

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
   type(part_node_t), pointer :: current, next_part, prev
   logical :: swapped

   ! Bubble sort implementation for linked list
   do
      swapped = .false.
      current => partition_root%first_part
      prev => null()
      do while (associated(current%next_part))
         next_part => current%next_part
         ! Check if we need to swap (current has more items than next)
         if (current%num_items1 > next_part%num_items1) then
            swapped = .true.
            ! Perform the swap
            current%next_part => next_part%next_part
            next_part%next_part => current
            if (associated(prev)) then
               prev%next_part => next_part
            else
               ! Update first_part if we're swapping the first element
               partition_root%first_part => next_part
            end if
            ! Update last_part if necessary
            if (.not. associated(current%next_part)) then
               partition_root%last_part => current
            end if
            ! Update prev for next iteration
            prev => next_part
         else
            ! No swap needed, just advance
            prev => current
            current => next_part
         end if
      end do
      ! If no swaps occurred, the list is sorted
      if (.not. swapped) exit
   end do
end subroutine

end module
