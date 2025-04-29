module lcrs_tree
use parameters
implicit none
private

! Item node
type, public :: item_node
   integer :: index
   type(item_node), pointer :: next_item
end type

! Root node
type, public :: root_node
   integer :: num_leaves
   type(leaf_node), pointer :: first_leaf
   type(leaf_node_ptr), dimension(:), allocatable :: itemdir1
   type(leaf_node_ptr), dimension(:), allocatable :: itemdir2
   type(root_node), pointer :: next_root
end type

! Leaf node
type, public :: leaf_node
   integer :: num_heirs
   integer :: num_items1
   integer :: num_items2
   type(root_node), pointer :: tree_root
   type(leaf_node), pointer :: next_leaf
   type(leaf_node), pointer :: first_heir
   type(leaf_node), pointer :: next_heir
   type(item_node), pointer :: first_item1
   type(item_node), pointer :: first_item2
   type(leaf_node_ptr), dimension(:), allocatable :: typehood
end type

! Leaf node pointer
type, public :: leaf_node_ptr
   type(leaf_node), pointer :: ptr
end type

! Polytree node
type, public :: poly_node
   integer :: num_roots
   integer :: size_itemdir1
   integer :: size_itemdir2
   type(root_node), pointer :: first_root
end type

! Semipartition part
type, public :: semipartition_part
   integer :: num_items
   integer, dimension(:), allocatable :: items
end type

! Item semipartition
type, public :: item_semipartition
   integer :: num_parts
   type(semipartition_part), dimension(:), allocatable :: parts
   integer, dimension(:), allocatable :: itemdir
end type

! Partition part
type, public :: partition_part
   integer :: num_heirs
   integer :: num_items1
   integer :: num_items2
   integer, dimension(:), allocatable :: heirs
   integer, dimension(:), allocatable :: items1
   integer, dimension(:), allocatable :: items2
   integer, dimension(:), allocatable :: typehood
end type

! Item partition
type, public :: item_partition
   integer :: num_parts
   type(partition_part), dimension(:), allocatable :: parts
   integer, dimension(:), allocatable :: itemdir1
   integer, dimension(:), allocatable :: itemdir2
end type

interface assignment(=)
   module procedure copy_tree
end interface

interface operator(==)
   module procedure leafnodeptr_equality
end interface

interface operator (.equiv.)
   module procedure typehood_equivalence
end interface

! Make types and procedures public
public make_new_root
public add_new_leaf
public delete_root
public add_heir
public find_heir
public add_new_item1
public add_new_item2
public move_next_item1
public move_next_item2
public move_node_items
public partition_from_tree
public tree_from_partition
public first_partition
public second_partition
public print_items
public print_itemdirs
public print_tree
public print_polytree
public print_partition
public make_new_poly
public add_new_root
public delete_polytree
public reverse_polytree
public assignment(=)
public operator(==)
public operator(.equiv.)

contains

! Copy a leaf node and all its items
function copy_leaf(src_leaf, dest_root) result(new_leaf)
   type(leaf_node), intent(in) :: src_leaf
   type(root_node), intent(inout) :: dest_root
   type(leaf_node), pointer :: new_leaf
   type(item_node), pointer :: src_item
   type(leaf_node_ptr) :: typehood(0)

   ! Create new leaf
   new_leaf => add_new_leaf(dest_root, typehood)
   new_leaf%num_heirs = src_leaf%num_heirs

   ! Copy first list of items
   src_item => src_leaf%first_item1
   do while (associated(src_item))
      call add_new_item1(new_leaf, src_item%index)
      src_item => src_item%next_item
   end do

   ! Copy second list of items
   src_item => src_leaf%first_item2
   do while (associated(src_item))
      call add_new_item2(new_leaf, src_item%index)
      src_item => src_item%next_item
   end do
end function

subroutine copy_tree(dest, src)
   type(root_node), pointer, intent(out) :: dest
   type(root_node), intent(in) :: src
   ! Local variables
   type(leaf_node), pointer :: src_leaf, new_leaf

   ! Create root node with same directory sizes as original
   dest => make_new_root(size(src%itemdir1), size(src%itemdir2))

   ! Copy all leaves
   src_leaf => src%first_leaf
   do while (associated(src_leaf))
      new_leaf => copy_leaf(src_leaf, dest)
      src_leaf => src_leaf%next_leaf
   end do
end subroutine

integer function address(node)
   use iso_c_binding, only: c_loc, c_intptr_t
   type(leaf_node), target, intent(in) :: node
   address = modulo(transfer(c_loc(node), c_intptr_t), 4096)
end function

elemental function leafnodeptr_equality(left, right) result(equality)
   type(leaf_node_ptr), intent(in) :: left, right
   logical :: equality
   equality = associated(left%ptr, right%ptr)
end function

function typehood_equivalence(array1, array2) result(equiv)
   type(leaf_node_ptr), dimension(:), intent(in) :: array1, array2
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

function make_new_root(tot_items1, tot_items2) result(new_root)
   integer, intent(in) :: tot_items1, tot_items2
   type(root_node), pointer :: new_root
   integer :: i

   allocate(new_root)
   new_root%num_leaves = 0
   new_root%first_leaf => null()
   new_root%next_root => null()

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

function add_new_leaf(root, typehood) result(new_leaf)
   type(root_node), target, intent(inout) :: root
   type(leaf_node_ptr), dimension(:), intent(in) :: typehood
   type(leaf_node), pointer :: new_leaf

   allocate (new_leaf)
   allocate (new_leaf%typehood, source=typehood)

   new_leaf%num_items1 = 0
   new_leaf%num_items2 = 0
   new_leaf%num_heirs = 0
   new_leaf%first_item1 => null()
   new_leaf%first_item2 => null()
   new_leaf%first_heir => null()
   new_leaf%next_heir => null()
   new_leaf%tree_root => root

   ! Add as first leaf
   new_leaf%next_leaf => root%first_leaf
   root%first_leaf => new_leaf
   root%num_leaves = root%num_leaves + 1
end function

subroutine add_new_item1(leaf, index)
   type(leaf_node), target, intent(inout) :: leaf
   integer, intent(in) :: index
   type(item_node), pointer :: new_item

   allocate(new_item)
   new_item%index = index
   leaf%tree_root%itemdir1(index)%ptr => leaf
   new_item%next_item => leaf%first_item1
   leaf%first_item1 => new_item
   leaf%num_items1 = leaf%num_items1 + 1
end subroutine

subroutine add_new_item2(leaf, index)
   type(leaf_node), target, intent(inout) :: leaf
   integer, intent(in) :: index
   type(item_node), pointer :: new_item

   allocate(new_item)
   new_item%index = index
   leaf%tree_root%itemdir2(index)%ptr => leaf
   new_item%next_item => leaf%first_item2
   leaf%first_item2 => new_item
   leaf%num_items2 = leaf%num_items2 + 1
end subroutine

subroutine move_next_item1(src, dest)
   type(leaf_node), intent(inout) :: src
   type(leaf_node), target, intent(inout) :: dest
   type(item_node), pointer :: src_second_item1, dest_first_item1

   src_second_item1 => src%first_item1%next_item
   dest_first_item1 => dest%first_item1
   dest%tree_root%itemdir1(src%first_item1%index)%ptr => dest
   dest%first_item1 => src%first_item1
   dest%first_item1%next_item => dest_first_item1
   src%first_item1 => src_second_item1
   src%num_items1 = src%num_items1 - 1
   dest%num_items1 = dest%num_items1 + 1
end subroutine

subroutine move_next_item2(src, dest)
   type(leaf_node), intent(inout) :: src
   type(leaf_node), target, intent(inout) :: dest
   type(item_node), pointer :: src_second_item2, dest_first_item2

   src_second_item2 => src%first_item2%next_item
   dest_first_item2 => dest%first_item2
   dest%tree_root%itemdir2(src%first_item2%index)%ptr => dest
   dest%first_item2 => src%first_item2
   dest%first_item2%next_item => dest_first_item2
   src%first_item2 => src_second_item2
   src%num_items2 = src%num_items2 - 1
   dest%num_items2 = dest%num_items2 + 1
end subroutine

subroutine move_node_items(src, dest)
   type(leaf_node), intent(inout) :: src, dest

   do while (associated(src%first_item1))
      call move_next_item1(src, dest)
   end do

   do while (associated(src%first_item2))
      call move_next_item2(src, dest)
   end do
end subroutine

subroutine deallocate_items(first_item)
   type(item_node), pointer, intent(inout) :: first_item
   type(item_node), pointer :: item, next_item

   item => first_item
   do while (associated(item))
      next_item => item%next_item
      deallocate(item)
      item => next_item
   end do

   first_item => null()
end subroutine

subroutine delete_leaf(leaf)
   type(leaf_node), pointer, intent(inout) :: leaf

   if (.not. associated(leaf)) return

   call deallocate_items(leaf%first_item1)
   call deallocate_items(leaf%first_item2)

   deallocate(leaf)
   leaf => null()
end subroutine

subroutine delete_root(root)
   type(root_node), pointer, intent(inout) :: root
   type(leaf_node), pointer :: leaf, next_leaf

   if (.not. associated(root)) return

   ! Delete all leaves
   leaf => root%first_leaf
   do while (associated(leaf))
      next_leaf => leaf%next_leaf
      call delete_leaf(leaf)
      leaf => next_leaf
   end do

   ! Deallocate directories
   deallocate(root%itemdir1)
   deallocate(root%itemdir2)

   ! Deallocate root
   deallocate(root)
   root => null()
end subroutine

function partition_from_tree(root) result(partition)
   type(root_node), target, intent(in) :: root
   type(item_partition) :: partition
   type(leaf_node), pointer :: leaf
   type(item_node), pointer :: item
   integer :: i, j

   partition%num_parts = root%num_leaves
   allocate(partition%parts(partition%num_parts))
   allocate(partition%itemdir1(size(root%itemdir1)))
   allocate(partition%itemdir2(size(root%itemdir2)))

   leaf => root%first_leaf
   do i = 1, partition%num_parts
      if (.not. associated(leaf)) error stop 'Unexpected null leaf'
      if (.not. associated(leaf%tree_root, root)) error stop 'Leaf points to wrong root'

      partition%parts(i)%num_items1 = leaf%num_items1
      partition%parts(i)%num_items2 = leaf%num_items2

      allocate(partition%parts(i)%items1(leaf%num_items1))
      allocate(partition%parts(i)%items2(leaf%num_items2))

      item => leaf%first_item1
      do j = 1, leaf%num_items1
         partition%parts(i)%items1(j) = item%index
         partition%itemdir1(item%index) = i
         item => item%next_item
      end do

      item => leaf%first_item2
      do j = 1, leaf%num_items2
         partition%parts(i)%items2(j) = item%index
         partition%itemdir2(item%index) = i
         item => item%next_item
      end do

      leaf => leaf%next_leaf
   end do
end function

function tree_from_partition(partition) result(root)
   type(item_partition), intent(in) :: partition
   type(root_node), pointer :: root
   type(leaf_node), pointer :: leaf
   type(leaf_node_ptr) :: typehood(0)
   integer :: i, j

   root => make_new_root(size(partition%itemdir1), size(partition%itemdir2))

   do i = 1, partition%num_parts
      leaf => add_new_leaf(root, typehood)

      do j = 1, partition%parts(i)%num_items1
         call add_new_item1(leaf, partition%parts(i)%items1(j))
      end do

      do j = 1, partition%parts(i)%num_items2
         call add_new_item2(leaf, partition%parts(i)%items2(j))
      end do
   end do
end function

function first_partition(bipartition) result(partition)
   type(item_partition), intent(in) :: bipartition
   type(item_semipartition) :: partition
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
   type(item_partition), intent(in) :: bipartition
   type(item_semipartition) :: partition
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

subroutine print_items(leaf)
   type(leaf_node), intent(in) :: leaf
   type(item_node), pointer :: item

   write(stderr,'(Z3.3)', advance='no') address(leaf)
   write(stderr, '(A)', advance='no') '('

   item => leaf%first_item1
   do while (associated(item))
      write(stderr, '(1X,I0)', advance='no') item%index
      item => item%next_item
   end do

   write(stderr, '(A)', advance='no') '|'

   item => leaf%first_item2
   do while (associated(item))
      write(stderr, '(1X,I0)', advance='no') item%index
      item => item%next_item
   end do

   write(stderr, '(A)', advance='no') ')'
end subroutine

subroutine print_itemdirs(root)
   type(root_node), intent(in) :: root
   integer :: i

   write(stderr, *)
   write(stderr,'(A)') "Item Directory 1:"
   do i = 1, size(root%itemdir1)
      if (associated(root%itemdir1(i)%ptr)) then
         write(stderr,'(2X,I0,A,Z3.3)') i, " -> ", address(root%itemdir1(i)%ptr)
      else
         write(stderr,'(2X,I0,A)') i, " -> Not associated"
      end if
   end do

   write(stderr, *)
   write(stderr,'(A)') "Item Directory 2:"
   do i = 1, size(root%itemdir2)
      if (associated(root%itemdir2(i)%ptr)) then
         write(stderr,'(2X,I0,A,Z3.3)') i, " -> ", address(root%itemdir2(i)%ptr)
      else
         write(stderr,'(2X,I0,A)') i, " -> Not associated"
      end if
   end do
end subroutine

subroutine print_tree(root)
   type(root_node), intent(in) :: root
   type(leaf_node), pointer :: leaf

   write(stderr, *)
   write(stderr,'(A)') "/types/"

   leaf => root%first_leaf
   do while (associated(leaf))
      call print_items(leaf)
      write(stderr, *)
      leaf => leaf%next_leaf
   end do
end subroutine

subroutine print_partition(partition)
   type(item_partition), intent(in) :: partition
   integer :: i, j

   write(stderr,*)
   write(stderr,'(A)') "Partition Parts:"
   do i = 1, partition%num_parts
      write(stderr,'(A,I0,A)',advance='no') "Part ", i, ": ("

      do j = 1, partition%parts(i)%num_items1
         write(stderr,'(1X,I0)',advance='no') partition%parts(i)%items1(j)
      end do

      write(stderr,'(A)',advance='no') "|"

      do j = 1, partition%parts(i)%num_items2
         write(stderr,'(1X,I0)',advance='no') partition%parts(i)%items2(j)
      end do

      write(stderr,'(A)') ")"
   end do
end subroutine

function make_new_poly(size_itemdir1, size_itemdir2)
   integer, intent(in) :: size_itemdir1, size_itemdir2
   type(poly_node), pointer :: make_new_poly

   allocate(make_new_poly)
   make_new_poly%num_roots = 0
   make_new_poly%size_itemdir1 = size_itemdir1
   make_new_poly%size_itemdir2 = size_itemdir2
   make_new_poly%first_root => null()
end function

function add_new_root(poly) result(new_root)
   type(poly_node), target, intent(inout) :: poly
   type(root_node), pointer :: new_root

   new_root => make_new_root(poly%size_itemdir1, poly%size_itemdir2)
   new_root%next_root => poly%first_root
   poly%first_root => new_root
   poly%num_roots = poly%num_roots + 1
end function

subroutine delete_polytree(poly)
   type(poly_node), pointer, intent(inout) :: poly
   type(root_node), pointer :: root, next_root

   if (.not. associated(poly)) return

   root => poly%first_root
   do while (associated(root))
      next_root => root%next_root
      call delete_root(root)
      root => next_root
   end do

   deallocate(poly)
   poly => null()
end subroutine

subroutine add_heir(leaf, heir_leaf)
   type(leaf_node), intent(inout) :: leaf
   type(leaf_node), pointer, intent(in) :: heir_leaf

   if (associated(leaf%first_heir)) then
      heir_leaf%next_heir => leaf%first_heir
   end if
   leaf%first_heir => heir_leaf
   leaf%num_heirs = leaf%num_heirs + 1
end subroutine

function find_heir(leaf, typehood) result(heir_leaf)
   type(leaf_node), intent(in) :: leaf
   type(leaf_node_ptr), dimension(:), intent(in) :: typehood
   type(leaf_node), pointer :: heir_leaf
   type(leaf_node), pointer :: heir

   heir => leaf%first_heir
   do while (associated(heir))
      if (heir%typehood .equiv. typehood) then
         heir_leaf => heir
         return
      end if
      heir => heir%next_heir
   end do

   heir_leaf => null()
end function

subroutine print_polytree(poly)
   type(poly_node), intent(in) :: poly
   type(root_node), pointer :: root
   type(leaf_node), pointer :: leaf
   type(leaf_node), pointer :: heir
   integer :: level_idx

   write(stderr,*)
   write(stderr,'(A)') "MNA Polytree Structure:"

   level_idx = 1
   root => poly%first_root

   ! Print each level
   do while (associated(root))
      write(stderr,*)
      write(stderr,'(A,I0)') "Level ", level_idx
      write(stderr,'(A)') "----------------"

      ! Print basic tree structure
      call print_tree(root)

      ! Print heir information for each leaf
      write(stderr,*)
      write(stderr,'(A)') "/subtypes/"

      leaf => root%first_leaf
      do while (associated(leaf))
         ! Print leaf address
         write(stderr,'(Z3.3,A)',advance='no') address(leaf), ":"

         ! Print each heir's connections
         heir => leaf%first_heir
         do while (associated(heir))
            ! Print which leaf in next level this heir points to
            write(stderr,'(1X,Z3.3)',advance='no') address(heir)
            heir => heir%next_heir
         end do
         write(stderr,*)

         leaf => leaf%next_leaf
      end do

      root => root%next_root
      level_idx = level_idx + 1
   end do

   write(stderr,*)  ! Final newline
end subroutine

subroutine reverse_polytree(poly)
    type(poly_node), intent(inout) :: poly
    type(root_node), pointer :: prev, curr, next

    ! Initialize pointers for reversal
    prev => null()
    curr => poly%first_root

    ! Reverse the linked list of roots
    do while (associated(curr))
        ! Store next root before changing links
        next => curr%next_root
        ! Reverse the current root's pointer
        curr%next_root => prev
        ! Move pointers one step forward
        prev => curr
        curr => next
    end do

    ! Update polytree's first root pointer to last root (which is now first)
    poly%first_root => prev
end subroutine

end module
