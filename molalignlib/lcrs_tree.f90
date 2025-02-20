module lcrs_tree
use iso_fortran_env, only: stdout => error_unit
implicit none
private

type, public :: partition_part
   integer :: num_items1
   integer :: num_items2
   integer, allocatable :: items1(:)
   integer, allocatable :: items2(:)
end type

type, public :: partition_container
   integer :: num_parts
   type(partition_part), allocatable :: parts(:)
   integer, allocatable :: itemdir1(:)
   integer, allocatable :: itemdir2(:)
end type

type, public :: item_node
   integer :: index
   type(item_node), pointer :: next_item
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
   type(tree_node_ptr), pointer :: itemdir1(:)
   type(tree_node_ptr), pointer :: itemdir2(:)
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
public delete_tree
public delete_descendants
public add_new_item1
public add_new_item2
public add_linked_item1
public add_linked_item2
public move_node_items
public print_tree
public print_subtree
public print_items
public print_itemdir
public partition_from_tree
public tree_from_partition
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
   new_item%next_item => node%first_item1
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
   new_item%next_item => node%first_item2
   node%first_item2 => new_item
   node%num_items2 = node%num_items2 + 1
end subroutine

subroutine add_linked_item1(node, item)
   type(tree_node), target, intent(inout) :: node
   type(item_node), pointer, intent(inout) :: item
   type(item_node), pointer :: next_item

   next_item => item%next_item
   node%itemdir1(item%index)%ptr => node
   item%next_item => node%first_item1
   node%first_item1 => item
   node%num_items1 = node%num_items1 + 1
   item => next_item
end subroutine

subroutine add_linked_item2(node, item)
   type(tree_node), target, intent(inout) :: node
   type(item_node), pointer, intent(inout) :: item
   type(item_node), pointer :: next_item

   next_item => item%next_item
   node%itemdir2(item%index)%ptr => node
   item%next_item => node%first_item2
   node%first_item2 => item
   node%num_items2 = node%num_items2 + 1
   item => next_item
end subroutine

subroutine move_node_items(from, dest)
   type(tree_node), intent(inout) :: from, dest
   type(item_node), pointer :: item

   item => from%first_item1
   do while (associated(item))
      call add_linked_item1(dest, item)
   end do

   item => from%first_item2
   do while (associated(item))
      call add_linked_item2(dest, item)
   end do

   from%num_items1 = 0
   from%num_items2 = 0
   from%first_item1 => null()
   from%first_item2 => null()
end subroutine

recursive subroutine delete_descendants(node)
   type(tree_node), intent(inout) :: node
   type(tree_node), pointer :: child, next_child
   type(item_node), pointer :: item, next_item

   child => node%first_child
   do while (associated(child))
      next_child => child%next_sibling

      ! Delete child's items
      item => child%first_item1
      do while (associated(item))
         next_item => item%next_item
         deallocate(item)
         item => next_item
      end do

      item => child%first_item2
      do while (associated(item))
         next_item => item%next_item
         deallocate(item)
         item => next_item
      end do

      call delete_descendants(child)

      if (node%num_childs > 1) then
         node%num_leaves = node%num_leaves - 1
      end if

      deallocate(child)
      child => next_child
      node%num_childs = node%num_childs - 1
   end do

   node%first_child => null()
end subroutine

subroutine delete_tree(root)
   type(tree_node), pointer, intent(inout) :: root

   ! Delete all descendants and their items
   call delete_descendants(root)

   ! Delete root's own items and shared resources
   deallocate(root%num_leaves)
   deallocate(root%itemdir1)
   deallocate(root%itemdir2)
   deallocate(root)
end subroutine

subroutine flat_tree_assign(flat_root, root)
   type(tree_node), pointer, intent(out) :: flat_root
   type(tree_node), target, intent(in) :: root

   ! Create root node with same directory sizes as original
   flat_root => make_new_root(size(root%itemdir1), size(root%itemdir2))

   ! Collect all items from the original tree into the root
   call collect_all_items(root, flat_root)

contains   
   recursive subroutine collect_all_items(src_node, dest_node)
      type(tree_node), intent(in) :: src_node
      type(tree_node), intent(inout) :: dest_node
      type(tree_node), pointer :: child
      type(item_node), pointer :: item

      ! Add items from current node
      item => src_node%first_item1
      do while (associated(item))
         call add_new_item1(dest_node, item%index)
         item => item%next_item
      end do

      item => src_node%first_item2
      do while (associated(item))
         call add_new_item2(dest_node, item%index)
         item => item%next_item
      end do

      ! Recursively collect items from children
      child => src_node%first_child
      do while (associated(child))
         call collect_all_items(child, dest_node)
         child => child%next_sibling
      end do
   end subroutine
end subroutine

subroutine partition_from_tree(root, partition)
   type(tree_node), intent(in) :: root
   type(partition_container), intent(out) :: partition
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
      type(partition_container), intent(inout) :: partition
      integer, intent(inout) :: leaf_idx
      type(tree_node), pointer :: child
      type(item_node), pointer :: item
      integer :: i

      if (.not. associated(node%first_child)) then
         ! This is a leaf - collect its items
         partition%parts(leaf_idx)%num_items1 = node%num_items1
         partition%parts(leaf_idx)%num_items2 = node%num_items2

         allocate(partition%parts(leaf_idx)%items1(node%num_items1))
         allocate(partition%parts(leaf_idx)%items2(node%num_items2))

         ! Collect items1 and update directory
         i = 1
         item => node%first_item1
         do while (associated(item))
            partition%parts(leaf_idx)%items1(i) = item%index
            partition%itemdir1(item%index) = leaf_idx
            i = i + 1
            item => item%next_item
         end do

         ! Collect items2 and update directory
         i = 1
         item => node%first_item2
         do while (associated(item))
            partition%parts(leaf_idx)%items2(i) = item%index
            partition%itemdir2(item%index) = leaf_idx
            i = i + 1
            item => item%next_item
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
   type(partition_container), intent(in) :: partition
   type(tree_node), pointer, intent(out) :: root
   integer :: i, j
   type(tree_node), pointer :: curr_node

   ! Create root node using partition's directory sizes
   root => make_new_root(size(partition%itemdir1), size(partition%itemdir2))

   ! For each partition part, create a leaf and add items
   do i = 1, partition%num_parts
      curr_node => add_new_child(root)

      ! Add items1
      do j = 1, partition%parts(i)%num_items1
         call add_new_item1(curr_node, partition%parts(i)%items1(j))
      end do

      ! Add items2
      do j = 1, partition%parts(i)%num_items2
         call add_new_item2(curr_node, partition%parts(i)%items2(j))
      end do
   end do
end subroutine

subroutine print_items(node)
   type(tree_node), intent(in) :: node
   type(item_node), pointer :: item

   if (.not. associated(node%first_item1) .and. &
       .not. associated(node%first_item2)) then
      write(stdout, '(A)', advance='no') '()'
      return
   end if

   write(stdout, '(A)', advance='no') '('

   ! Print first list
   item => node%first_item1
   do while (associated(item))
      write(stdout, '(1X,I0)', advance='no') item%index
      item => item%next_item
   end do

   write(stdout, '(A)', advance='no') ' |'

   ! Print second list
   item => node%first_item2
   do while (associated(item))
      write(stdout, '(1X,I0)', advance='no') item%index
      item => item%next_item
   end do

   write(stdout, '(A)', advance='no') ' )'
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
         write(stdout, '(A,A)', advance='no') repeat('     ', indent)
         call print_subtree(child, indent + 1)
         child => child%next_sibling
      end do
   end if

   ! Leaf node
   call print_items(node)
   write (stdout, *)
end subroutine

subroutine print_itemdir(itemdir)
   type(tree_node_ptr), intent(in) :: itemdir(:)
   integer :: i

   write(stdout,*)
   do i = 1, size(itemdir)
      if (associated(itemdir(i)%ptr)) then
         write(stdout,'(A,I0,A)') "  Item ", i, " -> Node ("
         call print_items(itemdir(i)%ptr)
         write(stdout,'(A)') ")"
      else
         write(stdout,'(A,I0,A)') "  Item ", i, " -> Not associated"
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
   type(partition_container), intent(in) :: partition
   integer :: i, j

   write(stdout,*)
   write(stdout,'(A)') "Partition Parts:"
   do i = 1, partition%num_parts
      write(stdout,'(A,I0,A)',advance='no') "Part ", i, ": ("

      ! Print items1
      do j = 1, partition%parts(i)%num_items1
         write(stdout,'(1X,I0)',advance='no') partition%parts(i)%items1(j)
      end do

      write(stdout,'(A)',advance='no') " |"

      ! Print items2
      do j = 1, partition%parts(i)%num_items2
         write(stdout,'(1X,I0)',advance='no') partition%parts(i)%items2(j)
      end do

      write(stdout,'(A)') " )"
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
