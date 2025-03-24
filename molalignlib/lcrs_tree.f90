module lcrs_tree
use iso_fortran_env, only: stdout => error_unit
use iso_c_binding, only: c_loc, c_intptr_t
implicit none
private

type, public :: partition_part
   integer :: num_items
   integer, dimension(:), allocatable :: items
end type

type, public :: partition_container
   integer :: num_parts
   type(partition_part), dimension(:), allocatable :: parts
   integer, dimension(:), allocatable :: itemdir
end type

type, public :: bipartition_part
   integer :: num_items1
   integer :: num_items2
   integer, dimension(:), allocatable :: indices1
   integer, dimension(:), allocatable :: indices2
end type

type, public :: bipartition_container
   integer :: num_parts
   type(bipartition_part), dimension(:), allocatable :: parts
   integer, dimension(:), allocatable :: itemdir1
   integer, dimension(:), allocatable :: itemdir2
end type

type, public :: item_node
   integer :: index
   type(item_node), pointer :: next
end type

type, public :: tree_node
   integer :: num_childs
   integer :: num_items1
   integer :: num_items2
   integer, pointer :: num_leaves
   type(tree_node), pointer :: first_child
   type(tree_node), pointer :: next_sibling
   type(item_node), pointer :: first_item1
   type(item_node), pointer :: first_item2
   type(tree_node_ptr), dimension(:), pointer :: itemdir1
   type(tree_node_ptr), dimension(:), pointer :: itemdir2
end type

type, public :: polytree_node
   integer :: num_roots
   integer :: size_itemdir1
   integer :: size_itemdir2
   type(tree_node), pointer :: first_root
end type

type, public :: tree_node_ptr
   type(tree_node), pointer :: ptr
end type

interface operator(==)
   module procedure treenodeptr_equality
end interface

interface assignment(=)
   module procedure flat_tree_assign
end interface

! Make types and procedures public
public make_new_root
public add_new_child
public flatten_tree
public delete_tree
public delete_subtree
public add_new_item1
public add_new_item2
public add_item1
public add_item2
public move_next_item1
public move_next_item2
public move_node_items
public partition_from_tree
public tree_from_partition
public first_partition
public second_partition
public print_items
public print_itemdir
public print_subtree
public print_tree
public print_partition
public operator(==)
public assignment(=)

! New public types and procedures
public :: new_polytree
public :: add_new_root
public :: delete_polytree
public :: partitionlist_from_polytree
public :: polytree_from_partition
contains

elemental function treenodeptr_equality(left, right) result(equality)
   type(tree_node_ptr), intent(in) :: left, right
   logical :: equality
   equality = associated(left%ptr, right%ptr)
end function

function make_new_root(tot_items1, tot_items2) result(new_root)
   integer, intent(in) :: tot_items1, tot_items2
   type(tree_node), pointer :: new_root
   integer :: i

   allocate(new_root)
   new_root%num_childs = 0
   new_root%num_items1 = 0
   new_root%num_items2 = 0
   new_root%first_child => null()
   new_root%next_sibling => null()
   new_root%first_item1 => null()
   new_root%first_item2 => null()

   ! Allocate num_leaves
   allocate(new_root%num_leaves)
   new_root%num_leaves = 1  ! Root starts as a leaf

   ! Allocate item directories
   allocate(new_root%itemdir1(tot_items1))
   allocate(new_root%itemdir2(tot_items2))

   ! Initialize all pointers to null
   do i = 1, tot_items1
      new_root%itemdir1(i)%ptr => null()
   end do
   do i = 1, tot_items2
      new_root%itemdir2(i)%ptr => null()
   end do
end function

function add_new_child(parent) result(new_child)
   type(tree_node), target, intent(inout) :: parent
   type(tree_node), pointer :: new_child

   allocate(new_child)
   new_child%num_childs = 0
   new_child%num_items1 = 0
   new_child%num_items2 = 0
   new_child%first_child => null()
   new_child%next_sibling => null()
   new_child%first_item1 => null()
   new_child%first_item2 => null()
   new_child%num_leaves => parent%num_leaves
   new_child%itemdir1 => parent%itemdir1
   new_child%itemdir2 => parent%itemdir2

   ! Add as first child
   new_child%next_sibling => parent%first_child
   parent%first_child => new_child
   parent%num_childs = parent%num_childs + 1

   ! Update leaf count
   if (parent%num_childs > 1) then
      new_child%num_leaves = new_child%num_leaves + 1
   end if
end function

subroutine add_new_item1(node, index)
   type(tree_node), target, intent(inout) :: node
   integer, intent(in) :: index
   type(item_node), pointer :: new_item

   allocate(new_item)
   new_item%index = index
   node%itemdir1(index)%ptr => node
   new_item%next => node%first_item1
   node%first_item1 => new_item
   node%num_items1 = node%num_items1 + 1
end subroutine

subroutine add_new_item2(node, index)
   type(tree_node), target, intent(inout) :: node
   integer, intent(in) :: index
   type(item_node), pointer :: new_item

   allocate(new_item)
   new_item%index = index
   node%itemdir2(index)%ptr => node
   new_item%next => node%first_item2
   node%first_item2 => new_item
   node%num_items2 = node%num_items2 + 1
end subroutine

subroutine add_item1(node, item)
   type(tree_node), target, intent(inout) :: node
   type(item_node), pointer, intent(inout) :: item

   node%itemdir1(item%index)%ptr => node
   item%next => node%first_item1
   node%first_item1 => item
   node%num_items1 = node%num_items1 + 1
end subroutine

subroutine add_item2(node, item)
   type(tree_node), target, intent(inout) :: node
   type(item_node), pointer, intent(inout) :: item

   node%itemdir2(item%index)%ptr => node
   item%next => node%first_item2
   node%first_item2 => item
   node%num_items2 = node%num_items2 + 1
end subroutine

subroutine move_next_item1(srce, dest)
   type(tree_node), intent(inout) :: srce
   type(tree_node), target, intent(inout) :: dest
   type(item_node), pointer :: srce_second_item1, dest_first_item1

   srce_second_item1 => srce%first_item1%next
   dest_first_item1 => dest%first_item1
   dest%itemdir1(srce%first_item1%index)%ptr => dest
   dest%first_item1 => srce%first_item1
   dest%first_item1%next => dest_first_item1
   srce%first_item1 => srce_second_item1
   srce%num_items1 = srce%num_items1 - 1
   dest%num_items1 = dest%num_items1 + 1
end subroutine

subroutine move_next_item2(srce, dest)
   type(tree_node), intent(inout) :: srce
   type(tree_node), target, intent(inout) :: dest
   type(item_node), pointer :: srce_second_item2, dest_first_item2

   srce_second_item2 => srce%first_item2%next
   dest_first_item2 => dest%first_item2
   dest%itemdir2(srce%first_item2%index)%ptr => dest
   dest%first_item2 => srce%first_item2
   dest%first_item2%next => dest_first_item2
   srce%first_item2 => srce_second_item2
   srce%num_items2 = srce%num_items2 - 1
   dest%num_items2 = dest%num_items2 + 1
end subroutine

subroutine move_node_items(srce, dest)
   type(tree_node), intent(inout) :: srce, dest

   do while (associated(srce%first_item1))
      call move_next_item1(srce, dest)
   end do

   do while (associated(srce%first_item2))
      call move_next_item2(srce, dest)
   end do
end subroutine

subroutine deallocate_items(first_item)
   type(item_node), pointer, intent(inout) :: first_item
   type(item_node), pointer :: item, next_item

   ! Traverse the list and deallocate each item
   item => first_item
   do while (associated(item))
      next_item => item%next
      deallocate(item)
      item => next_item
   end do

   ! Set the head pointer to null
   first_item => null()
end subroutine deallocate_items

recursive subroutine delete_subtree(node)
   type(tree_node), pointer, intent(inout) :: node
   type(tree_node), pointer :: child, next_child

   if (.not. associated(node)) error stop 'Node not associated'

   if (associated(node%first_child)) then
      ! Internal node - Delete all children recursively
      child => node%first_child
      do
         next_child => child%next_sibling
         call delete_subtree(child)
         if (.not. associated(next_child)) exit
         child => next_child
!         node%num_leaves = node%num_leaves - 1
      end do
   end if

   ! Delete all items of the node
   call deallocate_items(node%first_item1)
   call deallocate_items(node%first_item2)

   ! Finally, deallocate the node itself
   deallocate(node)
   node => null()
end subroutine delete_subtree

subroutine delete_tree(root)
   type(tree_node), pointer, intent(inout) :: root

   if (.not. associated(root)) error stop 'Root not associated'

   ! Deallocate shared resources
   deallocate(root%num_leaves)
   deallocate(root%itemdir1)
   deallocate(root%itemdir2)

   ! Delete the entire tree structure
   call delete_subtree(root)
end subroutine

subroutine flatten_tree(root)
   type(tree_node), pointer, intent(inout) :: root
   type(tree_node), pointer :: temp_root

   if (.not. associated(root)) error stop 'Root not associated'

   allocate(temp_root)
   temp_root%first_child => root%first_child
   root%first_child => null()

   ! Process the temporary root subtree and collect leaf nodes
   call traverse_subtree(temp_root, root)

   root%num_childs = root%num_leaves
   deallocate(temp_root)

contains
   recursive subroutine traverse_subtree(node, root)
      type(tree_node), pointer, intent(inout) :: node, root
      type(tree_node), pointer :: child, next_child

      ! Process all children of this node
      child => node%first_child
      do while (associated(child))
         next_child => child%next_sibling

         if (associated(child%first_child)) then
            ! Internal node - process its children recursively
            call traverse_subtree(child, root)

            ! After processing, deallocate this internal node's items
            call deallocate_items(child%first_item1)
            call deallocate_items(child%first_item2)

            ! Deallocate the internal node
            deallocate(child)
         else
            ! Leaf node - add to the beginning of the leaf list
            child%next_sibling => root%first_child
            root%first_child => child
         end if

         child => next_child
      end do
   end subroutine traverse_subtree
end subroutine flatten_tree

subroutine flat_tree_assign(flat_root, root)
   type(tree_node), pointer, intent(out) :: flat_root
   type(tree_node), target, intent(in) :: root
   type(tree_node), pointer :: child, new_child
   type(item_node), pointer :: item
   
   ! Create root node with same directory sizes as original
   flat_root => make_new_root(size(root%itemdir1), size(root%itemdir2))
   
   ! Process each leaf node and create a corresponding child in flat tree
   child => root%first_child
   do while (associated(child))
      ! Create a new child in the flat tree
      new_child => add_new_child(flat_root)
      
      ! Copy items from first list
      item => child%first_item1
      do while (associated(item))
         call add_new_item1(new_child, item%index)
         item => item%next
      end do
      
      ! Copy items from second list
      item => child%first_item2
      do while (associated(item))
         call add_new_item2(new_child, item%index)
         item => item%next
      end do
      
      ! Move to next child
      child => child%next_sibling
   end do
end subroutine

subroutine partition_from_tree(root, partition)
   type(tree_node), intent(in) :: root
   type(bipartition_container), intent(out) :: partition
   integer :: i, j
   type(tree_node), pointer :: child
   type(item_node), pointer :: item

   ! Set number of parts equal to number of leaf nodes (children of root)
   partition%num_parts = root%num_childs
   allocate(partition%parts(partition%num_parts))

   ! Allocate item directories same size as tree's directories
   allocate(partition%itemdir1(size(root%itemdir1)))
   allocate(partition%itemdir2(size(root%itemdir2)))

   ! Process each child (leaf) directly
   child => root%first_child
   do i = 1, partition%num_parts
      if (.not. associated(child)) error stop 'Unexpected null child'

      ! Set item counts for this part
      partition%parts(i)%num_items1 = child%num_items1
      partition%parts(i)%num_items2 = child%num_items2

      ! Allocate arrays for indices
      allocate(partition%parts(i)%indices1(child%num_items1))
      allocate(partition%parts(i)%indices2(child%num_items2))

      ! Copy indices1 and update directory
      item => child%first_item1
      do j = 1, child%num_items1
         partition%parts(i)%indices1(j) = item%index
         partition%itemdir1(item%index) = i
         item => item%next
      end do

      ! Copy indices2 and update directory
      item => child%first_item2
      do j = 1, child%num_items2
         partition%parts(i)%indices2(j) = item%index
         partition%itemdir2(item%index) = i
         item => item%next
      end do

      ! Move to next child
      child => child%next_sibling
   end do
end subroutine

subroutine tree_from_partition(partition, root)
   type(bipartition_container), intent(in) :: partition
   type(tree_node), pointer, intent(out) :: root
   integer :: i, j
   type(tree_node), pointer :: curr_node

   ! Create root node using partition's directory sizes
   root => make_new_root(size(partition%itemdir1), size(partition%itemdir2))

   ! For each partition part, create a leaf and add items
   do i = 1, partition%num_parts
      curr_node => add_new_child(root)

      ! Add indices1
      do j = 1, partition%parts(i)%num_items1
         call add_new_item1(curr_node, partition%parts(i)%indices1(j))
      end do

      ! Add indices2
      do j = 1, partition%parts(i)%num_items2
         call add_new_item2(curr_node, partition%parts(i)%indices2(j))
      end do
   end do
end subroutine

function first_partition(bipartition) result(partition)
   type(bipartition_container), intent(in) :: bipartition
   type(partition_container) :: partition
   integer :: i

   ! Set number of parts equal to bipartition's number of parts
   partition%num_parts = bipartition%num_parts
   allocate(partition%parts(partition%num_parts))

   ! Allocate item directory same size as bipartition's first directory
   allocate(partition%itemdir(size(bipartition%itemdir1)))

   ! Copy item directory
   partition%itemdir = bipartition%itemdir1

   ! Copy parts - only first items
   do i = 1, partition%num_parts
      partition%parts(i)%num_items = bipartition%parts(i)%num_items1

      ! Allocate and copy items array
      allocate(partition%parts(i)%items(partition%parts(i)%num_items))
      partition%parts(i)%items = bipartition%parts(i)%indices1
   end do
end function

function second_partition(bipartition) result(partition)
   type(bipartition_container), intent(in) :: bipartition
   type(partition_container) :: partition
   integer :: i

   ! Set number of parts equal to bipartition's number of parts
   partition%num_parts = bipartition%num_parts
   allocate(partition%parts(partition%num_parts))

   ! Allocate item directory same size as bipartition's second directory
   allocate(partition%itemdir(size(bipartition%itemdir2)))

   ! Copy item directory
   partition%itemdir = bipartition%itemdir2

   ! Copy parts - only second items
   do i = 1, partition%num_parts
      partition%parts(i)%num_items = bipartition%parts(i)%num_items2

      ! Allocate and copy items array
      allocate(partition%parts(i)%items(partition%parts(i)%num_items))
      partition%parts(i)%items = bipartition%parts(i)%indices2
   end do
end function

subroutine print_items(node)
   type(tree_node), intent(in) :: node
   type(item_node), pointer :: item

   write(stdout, '(A)', advance='no') '('

   ! Print first list
   item => node%first_item1
   do while (associated(item))
      write(stdout, '(1X,I2)', advance='no') item%index
      item => item%next
   end do

   write(stdout, '(A)', advance='no') '|'

   ! Print second list
   item => node%first_item2
   do while (associated(item))
      write(stdout, '(1X,I2)', advance='no') item%index
      item => item%next
   end do

   write(stdout, '(A)', advance='no') ')'
end subroutine

recursive subroutine print_subtree(node, indent)
   type(tree_node), intent(in) :: node
   integer, intent(in) :: indent
   type(tree_node), pointer :: child

   if (associated(node%first_child)) then
      child => node%first_child
      call print_items(node)
      write(stdout, '(A)', advance='no') '---'
      call print_subtree(child, indent + 1)

      child => child%next_sibling
      do
         if (.not. associated(child)) return
         write(stdout, '(A,A)', advance='no') repeat('      ', indent)
         call print_subtree(child, indent + 1)
         child => child%next_sibling
      end do
   end if

   ! Leaf node
   call print_items(node)
   write (stdout, *)
end subroutine

subroutine print_itemdir(itemdir)
   type(tree_node_ptr), dimension(:), intent(in) :: itemdir
   integer :: i

   do i = 1, size(itemdir)
      if (associated(itemdir(i)%ptr)) then
         write(stdout,'(A,I2,A,Z8)') "  Item ", i, " -> Node ", transfer(c_loc(itemdir(i)%ptr), c_intptr_t)
      else
         write(stdout,'(A,I2,A)') "  Item ", i, " -> Not associated"
      end if
   end do
end subroutine

subroutine print_tree(root)
   type(tree_node), intent(in) :: root

   ! Print the tree structure
   write(stdout, *)
   write(stdout,'(A)') "Tree Structure:"
   call print_subtree(root, 1)

   ! Print item directories
!   write(stdout,'(A)') "Item Directory 1:"
!   call print_itemdir(root%itemdir1)
!   write(stdout,'(A)') "Item Directory 2:"
!   call print_itemdir(root%itemdir2)
end subroutine

subroutine print_partition(partition)
   type(bipartition_container), intent(in) :: partition
   integer :: i, j

   write(stdout,*)
   write(stdout,'(A)') "Partition Parts:"
   do i = 1, partition%num_parts
      write(stdout,'(A,I0,A)',advance='no') "Part ", i, ": ("

      ! Print indices1
      do j = 1, partition%parts(i)%num_items1
         write(stdout,'(1X,I2)',advance='no') partition%parts(i)%indices1(j)
      end do

      write(stdout,'(A)',advance='no') "|"

      ! Print indices2
      do j = 1, partition%parts(i)%num_items2
         write(stdout,'(1X,I2)',advance='no') partition%parts(i)%indices2(j)
      end do

      write(stdout,'(A)') ")"
   end do

!   ! Print item directories
!   write(stdout,*)
!   write(stdout,'(A)') "Item Directory 1:"
!   do i = 1, size(partition%itemdir1)
!      write(stdout,'(A,I0,A,I0)') "  Item ", i, " -> Part ", partition%itemdir1(i)
!   end do
!
!   write(stdout,*)
!   write(stdout,'(A)') "Item Directory 2:"
!   do i = 1, size(partition%itemdir2)
!      write(stdout,'(A,I0,A,I0)') "  Item ", i, " -> Part ", partition%itemdir2(i)
!   end do
end subroutine

! Makes a new polytree root with no trees
function new_polytree(size_itemdir1, size_itemdir2)
   integer, intent(in) :: size_itemdir1, size_itemdir2
   type(polytree_node), pointer :: new_polytree
   
   allocate(new_polytree)
   new_polytree%num_roots = 0
   new_polytree%size_itemdir1 = size_itemdir1
   new_polytree%size_itemdir2 = size_itemdir2
   new_polytree%first_root => null()
end function

! Adds a new root to a polytree
function add_new_root(polytree) result(new_root)
   type(polytree_node), target, intent(inout) :: polytree
   type(tree_node), pointer :: new_root
   
   ! Create new root with itemdir sizes inherited from polytree
   new_root => make_new_root(polytree%size_itemdir1, polytree%size_itemdir2)
   
   ! Add to beginning of root list (linked list of roots)
   new_root%next_sibling => polytree%first_root
   polytree%first_root => new_root
   polytree%num_roots = polytree%num_roots + 1
end function

! Deletes an entire polytree
subroutine delete_polytree(polytree)
   type(polytree_node), pointer, intent(inout) :: polytree
   type(tree_node), pointer :: root, next_root
   
   if (.not. associated(polytree)) error stop 'Polytree not associated'
   
   ! Delete all tree roots
   root => polytree%first_root
   do while (associated(root))
      next_root => root%next_sibling
      call delete_tree(root)
      root => next_root
   end do
   
   ! Finally, deallocate the polytree itself
   deallocate(polytree)
   polytree => null()
end subroutine

! Creates an array of partition containers from a polytree
subroutine partitionlist_from_polytree(polytree, partitions)
   type(polytree_node), intent(in) :: polytree
   type(bipartition_container), dimension(:), allocatable, intent(out) :: partitions
   type(tree_node), pointer :: root
   integer :: i
   
   ! Allocate array of bipartitions with size equal to number of roots
   allocate(partitions(polytree%num_roots))
   
   ! Process each root in the polytree
   root => polytree%first_root
   i = polytree%num_roots
   
   do while (associated(root))
      ! Convert the current tree to a partition (process in reverse order)
      call partition_from_tree(root, partitions(i))
      
      ! Move to next root
      root => root%next_sibling
      i = i - 1
   end do
end subroutine

! Creates a polytree with a single tree from a partition container
subroutine polytree_from_partition(partition, polytree)
   type(bipartition_container), intent(in) :: partition
   type(polytree_node), pointer, intent(out) :: polytree
   type(tree_node), pointer :: root
   
   ! Create a polytree with appropriate directory sizes
   polytree => new_polytree(size(partition%itemdir1), size(partition%itemdir2))
   
   ! Create a tree from the partition
   call tree_from_partition(partition, root)
   
   ! Add the tree to the polytree
   root%next_sibling => polytree%first_root
   polytree%first_root => root
   polytree%num_roots = 1
end subroutine

end module
