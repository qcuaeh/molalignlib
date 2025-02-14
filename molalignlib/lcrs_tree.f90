module lcrs_tree
use iso_fortran_env, only: stdout => error_unit
implicit none
private

type, public :: item_node
   integer :: index
   type(item_node), pointer :: next_item
   type(tree_node_ptr), pointer :: dirloc
end type

type, public :: tree_node
   type(tree_node), pointer :: first_child
   type(tree_node), pointer :: last_child
   type(tree_node), pointer :: next_sibling
   type(tree_node), pointer :: prev_sibling
   type(item_node), pointer :: first_item1
   type(item_node), pointer :: first_item2
   integer :: child_count
   integer :: item_count
end type

type, public :: tree_node_ptr
   type(tree_node), pointer :: ptr
end type

type, public :: type_tree
   type(tree_node), pointer :: tree_root
   type(tree_node_ptr), pointer :: itemdir1(:)
   type(tree_node_ptr), pointer :: itemdir2(:)
end type

interface operator (==)
   module procedure treenodeptr_equality
end interface

! Make types and procedures public
public make_new_root
public add_new_child
public prune_branch
public delete_tree
public flatten_tree
public add_new_item1
public add_new_item2
public move_item1
public move_item2
public move_node_items
public print_tree
public print_items
public operator (==)

contains

elemental function treenodeptr_equality(left, right) result(equality)
   type(tree_node_ptr), intent(in) :: left, right
   logical :: equality
   equality = associated(left%ptr, right%ptr)
end function

function make_new_root() result(new_root)
   type(tree_node), pointer :: new_root
   allocate(new_root)
   ! Initialize all pointers and counters
   new_root%first_child => null()
   new_root%last_child => null()
   new_root%next_sibling => null()
   new_root%prev_sibling => null()
   new_root%first_item1 => null()
   new_root%first_item2 => null()
   new_root%child_count = 0
   new_root%item_count = 0
end function

function add_new_child(parent) result(new_child)
   type(tree_node), target, intent(inout) :: parent
   type(tree_node), pointer :: new_child

   allocate(new_child)
   ! Initialize common elements
   new_child%first_child => null()
   new_child%last_child => null()
   new_child%first_item1 => null()
   new_child%first_item2 => null()
   new_child%child_count = 0
   new_child%item_count = 0

   ! If parent has no children, make this the first child
   if (.not. associated(parent%first_child)) then
      ! First child initialization
      new_child%next_sibling => null()
      new_child%prev_sibling => null()
      parent%first_child => new_child
      parent%last_child => new_child
   else
      ! Sibling initialization
      new_child%next_sibling => null()
      new_child%prev_sibling => parent%last_child
      parent%last_child%next_sibling => new_child
      parent%last_child => new_child
   end if

   parent%child_count = parent%child_count + 1
end function

subroutine add_new_item1(node, index, itemdir)
   type(tree_node), target, intent(inout) :: node
   type(tree_node_ptr), target, intent(inout) :: itemdir(:)
   integer, intent(in) :: index
   type(item_node), pointer :: new_item

   allocate(new_item)
   new_item%index = index
   itemdir(index)%ptr => node
   new_item%dirloc => itemdir(index)
   new_item%next_item => node%first_item1
   node%first_item1 => new_item
   node%item_count = node%item_count + 1
end subroutine

subroutine add_new_item2(node, index, itemdir)
   type(tree_node), target, intent(inout) :: node
   type(tree_node_ptr), target, intent(inout) :: itemdir(:)
   integer, intent(in) :: index
   type(item_node), pointer :: new_item

   allocate(new_item)
   new_item%index = index
   itemdir(index)%ptr => node
   new_item%dirloc => itemdir(index)
   new_item%next_item => node%first_item2
   node%first_item2 => new_item
   node%item_count = node%item_count + 1
end subroutine

subroutine move_item1(node, item)
   type(tree_node), target, intent(inout) :: node
   type(item_node), pointer, intent(inout) :: item
   type(item_node), pointer :: next_item

   next_item => item%next_item
   item%dirloc%ptr => node
   item%next_item => node%first_item1
   node%first_item1 => item
   node%item_count = node%item_count + 1
   item => next_item
end subroutine

subroutine move_item2(node, item)
   type(tree_node), target, intent(inout) :: node
   type(item_node), pointer, intent(inout) :: item
   type(item_node), pointer :: next_item

   next_item => item%next_item
   item%dirloc%ptr => node
   item%next_item => node%first_item2
   node%first_item2 => item
   node%item_count = node%item_count + 1
   item => next_item
end subroutine

subroutine move_node_items(from, dest)
   type(tree_node), intent(inout) :: from, dest
   type(item_node), pointer :: item

   item => from%first_item1
   from%first_item1 => null()
   do while (associated(item))
      call move_item1(dest, item)
   end do

   item => from%first_item2
   from%first_item2 => null()
   do while (associated(item))
      call move_item2(dest, item)
   end do
end subroutine

recursive subroutine delete_tree(node)
   type(tree_node), pointer :: node
   type(item_node), pointer :: item, next_item
   type(tree_node), pointer :: child, next_child

   if (.not. associated(node)) return

   ! Delete first list items
   item => node%first_item1
   do while (associated(item))
      next_item => item%next_item
      deallocate(item)
      item => next_item
   end do

   ! Delete second list items
   item => node%first_item2
   do while (associated(item))
      next_item => item%next_item
      deallocate(item)
      item => next_item
   end do

   ! Delete children
   child => node%first_child
   do while (associated(child))
      next_child => child%next_sibling
      call delete_tree(child)
      child => next_child
   end do

   deallocate(node)
   node => null()
end subroutine

subroutine prune_branch(node)
   type(tree_node) :: node
   type(tree_node), pointer :: child, next_child

   child => node%first_child
   do while (associated(child))
      next_child => child%next_sibling
      call delete_tree(child)
      child => next_child
   end do

   node%first_child => null()
   node%child_count = 0
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

   write(stdout, '(1X,A)', advance='no') '|'

   ! Print second list
   item => node%first_item2
   do while (associated(item))
      write(stdout, '(1X,I0)', advance='no') item%index
      item => item%next_item
   end do

   write(stdout, '(1X,A)', advance='no') ')'
end subroutine

subroutine print_tree(node)
   type(tree_node), intent(in) :: node
   write(stdout, *)
   call traverse_tree(node, 1)

contains
   recursive subroutine traverse_tree(node, indent)
      type(tree_node), intent(in) :: node
      integer, intent(in) :: indent
      type(tree_node), pointer :: child

      if (associated(node%first_child)) then
         child => node%first_child
         call print_items(node)
         write(stdout, '(A)', advance='no') '---'
         call traverse_tree(child, indent + 1)

         child => child%next_sibling
         do
            if (.not. associated(child)) return
            write(stdout, '(A,A)', advance='no') repeat('     ', indent)
            call traverse_tree(child, indent + 1)
            child => child%next_sibling
         end do
      end if

      ! Leaf node
      call print_items(node)
      write (stdout, *)
   end subroutine
end subroutine

subroutine flatten_tree(tree)
   type(type_tree), target, intent(inout) :: tree
   type(tree_node), pointer :: new_root
   integer :: total_comparisons

   ! Create temporary root for flattened tree
   new_root => make_new_root()
   total_comparisons = 0

   ! Collect leaves and sort them by total item count
   call collect_and_sort_leaves(tree%tree_root, new_root, total_comparisons)
!   write(stdout,*) "Total comparisons performed:", total_comparisons

   ! Clean up original tree and update root
   call delete_tree(tree%tree_root)
   tree%tree_root => new_root

contains
   recursive subroutine collect_and_sort_leaves(node, flat_parent, comparisons)
      type(tree_node), target, intent(inout) :: node
      type(tree_node), target, intent(inout) :: flat_parent
      integer, intent(inout) :: comparisons
      type(tree_node), pointer :: child, next_child, new_leaf
      type(tree_node), pointer :: curr

      if (.not. associated(node%first_child)) then
         ! Leaf node - create new leaf
         new_leaf => make_new_root()
         call move_node_items(node, new_leaf)

         ! First element case
         if (.not. associated(flat_parent%first_child)) then
            flat_parent%first_child => new_leaf
            flat_parent%last_child => new_leaf
            new_leaf%prev_sibling => null()
         else
            ! Compare with last element first
            comparisons = comparisons + 1
            curr => flat_parent%last_child
            if (curr%item_count <= new_leaf%item_count) then
               ! Insert at end
               curr%next_sibling => new_leaf
               new_leaf%prev_sibling => curr
               flat_parent%last_child => new_leaf
            else
               ! Scan backwards for insertion point
               do while (associated(curr))
                  if (.not. associated(curr%prev_sibling) .or. &
                      curr%prev_sibling%item_count <= new_leaf%item_count) exit
                  comparisons = comparisons + 1
                  curr => curr%prev_sibling
               end do

               ! Insert before curr
               new_leaf%next_sibling => curr
               new_leaf%prev_sibling => curr%prev_sibling
               if (associated(curr%prev_sibling)) then
                  curr%prev_sibling%next_sibling => new_leaf
               else
                  flat_parent%first_child => new_leaf
               end if
               curr%prev_sibling => new_leaf
            end if
         end if
         flat_parent%child_count = flat_parent%child_count + 1
      else
         ! Process children
         child => node%first_child
         do while (associated(child))
            next_child => child%next_sibling
            call collect_and_sort_leaves(child, flat_parent, comparisons)
            child => next_child
         end do
         node%first_child => null()
         node%last_child => null()
      end if
   end subroutine
end subroutine

end module
