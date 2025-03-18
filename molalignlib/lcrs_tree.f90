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

   ! Create root node with same directory sizes as original
   flat_root => make_new_root(size(root%itemdir1), size(root%itemdir2))

   ! Find all leaves in the original tree and add them as direct children to the flat tree
   call collect_leaf_nodes(root, flat_root)

contains   
   recursive subroutine collect_leaf_nodes(src_node, dest_root)
      type(tree_node), intent(in) :: src_node
      type(tree_node), intent(inout) :: dest_root
      type(tree_node), pointer :: child
      type(item_node), pointer :: item
      type(tree_node), pointer :: new_child

      if (associated(src_node%first_child)) then
         ! Internal node - process children
         child => src_node%first_child
         do
            call collect_leaf_nodes(child, dest_root)
            if (.not. associated(child%next_sibling)) exit
            child => child%next_sibling
         end do
      else
         ! Leaf node - create a new child in the flat tree
         new_child => add_new_child(dest_root)
         ! Copy items from original leaf to new child
         item => src_node%first_item1
         do while (associated(item))
            call add_new_item1(new_child, item%index)
            item => item%next
         end do
         item => src_node%first_item2
         do while (associated(item))
            call add_new_item2(new_child, item%index)
            item => item%next
         end do
      end if
   end subroutine
end subroutine

subroutine partition_from_tree(root, partition)
   type(tree_node), intent(in) :: root
   type(bipartition_container), intent(out) :: partition
   integer :: leaf_index

   ! Set number of parts equal to number of leaves
   partition%num_parts = root%num_leaves
   allocate(partition%parts(partition%num_parts))

   ! Allocate item directories same size as tree's directories
   allocate(partition%itemdir1(size(root%itemdir1)))
   allocate(partition%itemdir2(size(root%itemdir2)))

   ! Initialize leaf index
   leaf_index = 1

   ! Traverse tree and collect items
   call collect_items(root, partition, leaf_index)

contains
   recursive subroutine collect_items(node, partition, leaf_idx)
      type(tree_node), intent(in) :: node
      type(bipartition_container), intent(inout) :: partition
      integer, intent(inout) :: leaf_idx
      type(tree_node), pointer :: child
      type(item_node), pointer :: item
      integer :: i

      if (.not. associated(node%first_child)) then
         ! This is a leaf - collect its items
         partition%parts(leaf_idx)%num_items1 = node%num_items1
         partition%parts(leaf_idx)%num_items2 = node%num_items2

         allocate(partition%parts(leaf_idx)%indices1(node%num_items1))
         allocate(partition%parts(leaf_idx)%indices2(node%num_items2))

         ! Collect indices1 and update directory
         i = 1
         item => node%first_item1
         do while (associated(item))
            partition%parts(leaf_idx)%indices1(i) = item%index
            partition%itemdir1(item%index) = leaf_idx
            i = i + 1
            item => item%next
         end do

         ! Collect indices2 and update directory
         i = 1
         item => node%first_item2
         do while (associated(item))
            partition%parts(leaf_idx)%indices2(i) = item%index
            partition%itemdir2(item%index) = leaf_idx
            i = i + 1
            item => item%next
         end do

         leaf_idx = leaf_idx + 1
      else
         ! Process children
         child => node%first_child
         do while (associated(child))
            call collect_items(child, partition, leaf_idx)
            child => child%next_sibling
         end do
      end if
   end subroutine
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

end module
