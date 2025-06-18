module array_trees
use parameters
use lcrs_tree
implicit none
private

! Arrays of derived types with scalar components
! REMOVED: item_array_t - replaced with pure arrays

type :: part_array_t
   integer :: depth
   integer :: num_children
   ! Relationships (0 = null)
   integer :: parent_part_idx
   integer :: first_child_idx
   integer :: last_child_idx
   integer :: next_sibling_idx
   ! Item segments in flattened arrays
   integer :: items1_start_idx, items1_count
   integer :: items2_start_idx, items2_count
   integer :: items1_fill_count, items2_fill_count  ! current fill level during redistribution
   ! Fixed-size signature (much simpler!)
   integer :: signature(MAX_COORD)
   integer :: signature_length
end type

! REMOVED: partref_array_t - replaced with pure array

type :: link_array_t
   integer :: num_parts
   integer :: parent_chain_idx    ! which chain owns this link
   ! Part reference segment in flattened array
   integer :: partref_start_idx   ! start of this link's partrefs in flattened array
   ! Item directories (flattened storage)
   integer :: itemdir1_start_idx
   integer :: itemdir2_start_idx
end type

type :: chain_array_t
   integer :: tot_items1, tot_items2
   integer :: num_links, num_children
   ! Cross-tree reference (0 = null)
   integer :: split_part_idx      ! points to part array
   ! Chain tree relationships (0 = null)
   integer :: parent_chain_idx
   integer :: first_child_idx
   integer :: last_child_idx
   integer :: next_sibling_idx
   ! Link segment in flattened array
   integer :: link_start_idx      ! start of this chain's links in flattened array
end type

! Complete array-based representation
type, public :: array_trees_t
   ! Pure arrays for item values (no linked lists!)
   integer, allocatable :: item1_values(:)
   integer, allocatable :: item2_values(:)
   type(part_array_t), allocatable :: parts(:)
   type(link_array_t), allocatable :: links(:)
   type(chain_array_t), allocatable :: chains(:)

   ! Flattened variable-length data - all pure integer arrays!
   integer, allocatable :: itemdir_entries(:)   ! Part indices for itemdir (no wrapper type!)
   integer, allocatable :: partref_entries(:)   ! Part indices for partrefs

   ! Metadata
   integer :: total_items1, total_items2, total_parts
   integer :: total_links, total_chains
   integer :: total_itemdir_entries, total_partref_entries
   integer :: itemdir_size1, itemdir_size2  ! size of each itemdir
end type

public convert_trees_to_arrays
public validate_conversion
public print_tree_items_array
public print_part_tree_array
public print_leaf_items_array
public print_first_level_items_array
public print_chain_tree_array
public print_part_signatures_array
public print_chain_details_array

contains

subroutine convert_trees_to_arrays(root_part, root_chain, array_trees)
   type(part_node_t), pointer, intent(in) :: root_part
   type(chain_node_t), pointer, intent(in) :: root_chain
   type(array_trees_t), intent(out) :: array_trees
   integer :: i

   ! Get totals from the tree counters
   array_trees%total_parts = root_part%total_parts
   array_trees%total_items1 = root_part%total_items1
   array_trees%total_items2 = root_part%total_items2
   array_trees%total_chains = root_chain%total_chains
   array_trees%total_links = root_chain%total_links
   array_trees%total_partref_entries = root_chain%total_partrefs  ! Same count, just stored differently
   array_trees%itemdir_size1 = root_chain%tot_items1
   array_trees%itemdir_size2 = root_chain%tot_items2

   ! Calculate itemdir storage needs
   array_trees%total_itemdir_entries = array_trees%total_links * &
                                      (array_trees%itemdir_size1 + array_trees%itemdir_size2)

   ! Allocate all arrays with exact sizes
   allocate(array_trees%item1_values(array_trees%total_items1))
   allocate(array_trees%item2_values(array_trees%total_items2))
   allocate(array_trees%parts(array_trees%total_parts))
   allocate(array_trees%links(array_trees%total_links))
   allocate(array_trees%chains(array_trees%total_chains))
   allocate(array_trees%itemdir_entries(array_trees%total_itemdir_entries))
   allocate(array_trees%partref_entries(array_trees%total_partref_entries))

   ! Initialize arrays to 0
   do i = 1, array_trees%total_itemdir_entries
      array_trees%itemdir_entries(i) = 0
   end do

   do i = 1, array_trees%total_partref_entries
      array_trees%partref_entries(i) = 0
   end do

   do i = 1, array_trees%total_items1
      array_trees%item1_values(i) = 0
   end do

   do i = 1, array_trees%total_items2
      array_trees%item2_values(i) = 0
   end do

   ! Convert using global indices
   call populate_arrays_direct(root_part, root_chain, array_trees)
end subroutine

subroutine populate_arrays_direct(root_part, root_chain, array_trees)
   type(part_node_t), pointer, intent(in) :: root_part
   type(chain_node_t), pointer, intent(in) :: root_chain
   type(array_trees_t), intent(inout) :: array_trees
   integer :: itemdir_idx, partref_idx, link_idx, item1_idx, item2_idx

   ! Convert part tree starting from root with global item tracking
   item1_idx = 0
   item2_idx = 0
   call convert_parts_recursive(root_part, array_trees, item1_idx, item2_idx)

   ! Convert chain tree starting from root with global tracking
   itemdir_idx = 0
   partref_idx = 0
   link_idx = 0
   call convert_chains_recursive(root_chain, array_trees, itemdir_idx, partref_idx, link_idx)
end subroutine

recursive subroutine convert_parts_recursive(part, array_trees, item1_idx, item2_idx)
   type(part_node_t), pointer, intent(in) :: part
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(inout) :: item1_idx, item2_idx
   type(part_node_t), pointer :: child_part
   type(item_node_t), pointer :: item
   integer :: part_idx, i

   if (.not. associated(part)) return

   ! Convert this part (global indices always start at 1)
   part_idx = part%global_index

   array_trees%parts(part_idx)%depth = part%depth
   array_trees%parts(part_idx)%num_children = part%num_children

   ! Relationships using global indices
   if (associated(part%parent_part)) then
      array_trees%parts(part_idx)%parent_part_idx = part%parent_part%global_index
   else
      array_trees%parts(part_idx)%parent_part_idx = 0
   end if

   if (associated(part%first_child_part)) then
      array_trees%parts(part_idx)%first_child_idx = part%first_child_part%global_index
   else
      array_trees%parts(part_idx)%first_child_idx = 0
   end if

   if (associated(part%last_child_part)) then
      array_trees%parts(part_idx)%last_child_idx = part%last_child_part%global_index
   else
      array_trees%parts(part_idx)%last_child_idx = 0
   end if

   if (associated(part%next_sibling_part)) then
      array_trees%parts(part_idx)%next_sibling_idx = part%next_sibling_part%global_index
   else
      array_trees%parts(part_idx)%next_sibling_idx = 0
   end if

   ! Set up item segments using pure array approach
   array_trees%parts(part_idx)%items1_start_idx = item1_idx + 1
   array_trees%parts(part_idx)%items1_count = part%num_items1
   array_trees%parts(part_idx)%items1_fill_count = part%num_items1  ! Initially fully filled

   array_trees%parts(part_idx)%items2_start_idx = item2_idx + 1
   array_trees%parts(part_idx)%items2_count = part%num_items2
   array_trees%parts(part_idx)%items2_fill_count = part%num_items2  ! Initially fully filled

   ! Convert signature to fixed-size array
   array_trees%parts(part_idx)%signature_length = size(part%signature)
   do i = 1, size(part%signature)
      if (associated(part%signature(i)%ptr)) then
         array_trees%parts(part_idx)%signature(i) = part%signature(i)%ptr%global_index
      else
         array_trees%parts(part_idx)%signature(i) = 0
      end if
   end do
   ! Zero out unused signature entries
   do i = size(part%signature) + 1, MAX_COORD
      array_trees%parts(part_idx)%signature(i) = 0
   end do

   ! Convert items1 to pure array format
   item => part%first_item1
   i = 0
   do while (associated(item))
      i = i + 1
      item1_idx = item1_idx + 1
      array_trees%item1_values(item1_idx) = item%value
      item => item%next_item
   end do

   ! Convert items2 to pure array format
   item => part%first_item2
   i = 0
   do while (associated(item))
      i = i + 1
      item2_idx = item2_idx + 1
      array_trees%item2_values(item2_idx) = item%value
      item => item%next_item
   end do

   ! Recursively convert all children
   child_part => part%first_child_part
   do while (associated(child_part))
      call convert_parts_recursive(child_part, array_trees, item1_idx, item2_idx)
      child_part => child_part%next_sibling_part
   end do
end subroutine

recursive subroutine convert_chains_recursive(chain, array_trees, itemdir_idx, partref_idx, link_idx)
   type(chain_node_t), pointer, intent(in) :: chain
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(inout) :: itemdir_idx, partref_idx, link_idx
   type(chain_node_t), pointer :: child_chain
   type(link_node_t), pointer :: link
   type(partref_node_t), pointer :: partref
   integer :: chain_idx, current_link_idx

   if (.not. associated(chain)) return

   ! Convert this chain
   chain_idx = chain%global_index

   array_trees%chains(chain_idx)%tot_items1 = chain%tot_items1
   array_trees%chains(chain_idx)%tot_items2 = chain%tot_items2
   array_trees%chains(chain_idx)%num_links = chain%num_links
   array_trees%chains(chain_idx)%num_children = chain%num_children

   ! Cross-tree reference
   if (associated(chain%split_part)) then
      array_trees%chains(chain_idx)%split_part_idx = chain%split_part%global_index
   else
      array_trees%chains(chain_idx)%split_part_idx = 0
   end if

   ! Chain relationships
   if (associated(chain%parent_chain)) then
      array_trees%chains(chain_idx)%parent_chain_idx = chain%parent_chain%global_index
   else
      array_trees%chains(chain_idx)%parent_chain_idx = 0
   end if

   if (associated(chain%first_child_chain)) then
      array_trees%chains(chain_idx)%first_child_idx = chain%first_child_chain%global_index
   else
      array_trees%chains(chain_idx)%first_child_idx = 0
   end if

   if (associated(chain%last_child_chain)) then
      array_trees%chains(chain_idx)%last_child_idx = chain%last_child_chain%global_index
   else
      array_trees%chains(chain_idx)%last_child_idx = 0
   end if

   if (associated(chain%next_sibling_chain)) then
      array_trees%chains(chain_idx)%next_sibling_idx = chain%next_sibling_chain%global_index
   else
      array_trees%chains(chain_idx)%next_sibling_idx = 0
   end if

   ! Set link start index for pure array approach
   array_trees%chains(chain_idx)%link_start_idx = link_idx + 1

   ! Convert links in this chain using pure array approach
   link => chain%first_link
   do while (associated(link))
      link_idx = link_idx + 1
      current_link_idx = link_idx

      array_trees%links(current_link_idx)%num_parts = link%num_parts
      array_trees%links(current_link_idx)%parent_chain_idx = chain%global_index

      ! Set partref start index (pure array approach)
      array_trees%links(current_link_idx)%partref_start_idx = partref_idx + 1

      ! Convert partrefs to pure array format
      partref => link%first_partref
      do while (associated(partref))
         partref_idx = partref_idx + 1
         array_trees%partref_entries(partref_idx) = partref%part%global_index
         partref => partref%nextref
      end do

      ! Set itemdir start indices
      array_trees%links(current_link_idx)%itemdir1_start_idx = itemdir_idx + 1
      itemdir_idx = itemdir_idx + array_trees%itemdir_size1

      array_trees%links(current_link_idx)%itemdir2_start_idx = itemdir_idx + 1
      itemdir_idx = itemdir_idx + array_trees%itemdir_size2

      link => link%next_link
   end do

   ! Recursively convert child chains
   child_chain => chain%first_child_chain
   do while (associated(child_chain))
      call convert_chains_recursive(child_chain, array_trees, itemdir_idx, partref_idx, link_idx)
      child_chain => child_chain%next_sibling_chain
   end do
end subroutine

! Simplified validation procedure
subroutine validate_conversion(root_part, root_chain, array_trees)
   type(part_node_t), pointer, intent(in) :: root_part
   type(chain_node_t), pointer, intent(in) :: root_chain
   type(array_trees_t), intent(in) :: array_trees
   logical :: validation_passed

   validation_passed = .true.

   write(stderr, '(A)') "=== CONVERSION VALIDATION ==="

   ! Validate counts
   if (array_trees%total_parts /= root_part%total_parts) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Part count mismatch: ", &
         array_trees%total_parts, " vs ", root_part%total_parts
      validation_passed = .false.
   end if

   if (array_trees%total_items1 /= root_part%total_items1) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Items1 count mismatch: ", &
         array_trees%total_items1, " vs ", root_part%total_items1
      validation_passed = .false.
   end if

   if (array_trees%total_items2 /= root_part%total_items2) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Items2 count mismatch: ", &
         array_trees%total_items2, " vs ", root_part%total_items2
      validation_passed = .false.
   end if

   if (array_trees%total_chains /= root_chain%total_chains) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Chain count mismatch: ", &
         array_trees%total_chains, " vs ", root_chain%total_chains
      validation_passed = .false.
   end if

   if (array_trees%total_links /= root_chain%total_links) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Link count mismatch: ", &
         array_trees%total_links, " vs ", root_chain%total_links
      validation_passed = .false.
   end if

   if (array_trees%total_partref_entries /= root_chain%total_partrefs) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Partref count mismatch: ", &
         array_trees%total_partref_entries, " vs ", root_chain%total_partrefs
      validation_passed = .false.
   end if

   if (validation_passed) then
      write(stderr, '(A)') "✓ Conversion validation PASSED"
   else
      write(stderr, '(A)') "✗ Conversion validation FAILED"
   end if
end subroutine

subroutine print_tree_items_array(array_trees)
   type(array_trees_t), intent(in) :: array_trees

   if (array_trees%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "      Part Items (Array)"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Root part is always at index 1, print its children recursively
   call print_items_recursive_array(array_trees, 1)

   write(stderr, *)
end subroutine

recursive subroutine print_items_recursive_array(array_trees, part_idx)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   integer :: child_idx

   ! Traverse children using first_child_idx and next_sibling_idx
   child_idx = array_trees%parts(part_idx)%first_child_idx
   do while (child_idx > 0 .and. child_idx <= array_trees%total_parts)
      ! Print the child items with part index prefix
      write(stderr, '(A,I0,A)', advance='no') "Part ", child_idx, ':'
      call print_part_items_array(array_trees, child_idx)

      ! Recursively print this child's children
      call print_items_recursive_array(array_trees, child_idx)

      ! Move to next sibling
      child_idx = array_trees%parts(child_idx)%next_sibling_idx
   end do
end subroutine

subroutine print_part_items_array(array_trees, part_idx)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   integer :: i, start_idx

   ! Print items1 using pure array access
   start_idx = array_trees%parts(part_idx)%items1_start_idx
   do i = 1, array_trees%parts(part_idx)%items1_count
      write(stderr, '(1X,I0)', advance='no') array_trees%item1_values(start_idx + i - 1)
   end do

   write(stderr, '(A)', advance='no') ' /'

   ! Print items2 using pure array access
   start_idx = array_trees%parts(part_idx)%items2_start_idx
   do i = 1, array_trees%parts(part_idx)%items2_count
      write(stderr, '(1X,I0)', advance='no') array_trees%item2_values(start_idx + i - 1)
   end do

   write(stderr, *)
end subroutine

! Array-based tree printing procedures

subroutine print_part_tree_array(array_trees)
   type(array_trees_t), intent(in) :: array_trees
   logical, dimension(:), allocatable :: is_last_child

   if (array_trees%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "    Part Tree (Array)"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Allocate tracking array for tree lines (max depth 100)
   allocate(is_last_child(100))
   is_last_child = .false.

   ! Print root line
   write(stderr, '(A)') 'ROOT'

   ! Print children recursively (root is always at index 1)
   call print_part_recursive_array(array_trees, 1, 0, is_last_child)

   deallocate(is_last_child)
   write(stderr, *)
end subroutine

recursive subroutine print_part_recursive_array(array_trees, part_idx, depth, is_last_child)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx, depth
   logical, dimension(:), intent(inout) :: is_last_child
   integer :: child_idx, next_child_idx, i, pos
   character(len=200) :: prefix

   if (part_idx == 0) return

   ! Process all children
   child_idx = array_trees%parts(part_idx)%first_child_idx
   do while (child_idx > 0)
      ! Check if this is the last child
      next_child_idx = array_trees%parts(child_idx)%next_sibling_idx
      is_last_child(depth + 1) = (next_child_idx == 0)

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

      ! Add branch characters
      if (is_last_child(depth + 1)) then
         prefix(pos:pos+2) = "`--"
      else
         prefix(pos:pos+2) = "|--"
      end if
      pos = pos + 3

      ! Print part index with item counts
      write(stderr, '(A,I0,1X,A,I0,A,I0,A)') prefix(1:pos-1), child_idx, &
         '(', array_trees%parts(child_idx)%items1_count, '/', &
         array_trees%parts(child_idx)%items2_count, ')'

      ! Recursively print this child's children
      call print_part_recursive_array(array_trees, child_idx, depth + 1, is_last_child)

      child_idx = next_child_idx
   end do
end subroutine

subroutine print_chain_tree_array(array_trees)
   type(array_trees_t), intent(in) :: array_trees
   logical, dimension(:), allocatable :: is_last_child

   if (array_trees%total_chains == 0) then
      write(stderr, '(A)') "Chain tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "   Chain Tree (Array)"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Allocate tracking array for tree lines (max depth 100)
   allocate(is_last_child(100))
   is_last_child = .false.

   ! Print root line
   write(stderr, '(A)') 'ROOT'

   ! Print children recursively (root is always at index 1)
   call print_chain_recursive_array(array_trees, 1, 0, is_last_child)

   deallocate(is_last_child)
   write(stderr, *)
end subroutine

recursive subroutine print_chain_recursive_array(array_trees, chain_idx, depth, is_last_child)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: chain_idx, depth
   logical, dimension(:), intent(inout) :: is_last_child
   integer :: child_idx, next_child_idx, split_part_idx, i, pos
   character(len=200) :: prefix

   if (chain_idx == 0) return

   ! Process all children
   child_idx = array_trees%chains(chain_idx)%first_child_idx
   do while (child_idx > 0)
      ! Check if this is the last child
      next_child_idx = array_trees%chains(child_idx)%next_sibling_idx
      is_last_child(depth + 1) = (next_child_idx == 0)

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

      ! Add branch characters
      if (is_last_child(depth + 1)) then
         prefix(pos:pos+2) = "`--"
      else
         prefix(pos:pos+2) = "|--"
      end if
      pos = pos + 3

      ! Print the split part index with item counts
      split_part_idx = array_trees%chains(child_idx)%split_part_idx
      if (split_part_idx > 0) then
         write(stderr, '(A,I0,1X,A,I0,A,I0,A)') prefix(1:pos-1), split_part_idx, &
            '(', array_trees%parts(split_part_idx)%items1_count, '/', &
            array_trees%parts(split_part_idx)%items2_count, ')'
      else
         write(stderr, '(A,A)') prefix(1:pos-1), '(no split part)'
      end if

      ! Recursively print this child's children
      call print_chain_recursive_array(array_trees, child_idx, depth + 1, is_last_child)

      child_idx = next_child_idx
   end do
end subroutine

subroutine print_part_signatures_array(array_trees)
   type(array_trees_t), intent(in) :: array_trees

   if (array_trees%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "  Part Signatures (Array)"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Print signatures for all parts (root is always at index 1)
   call print_signatures_recursive_array(array_trees, 1)

   write(stderr, *)
end subroutine

recursive subroutine print_signatures_recursive_array(array_trees, part_idx)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   integer :: child_idx, i

   if (part_idx == 0) return

   ! Process all children in the same order as print_part_tree
   child_idx = array_trees%parts(part_idx)%first_child_idx
   do while (child_idx > 0)
      ! Print the child signature
      write(stderr,'(A,I0,A)',advance='no') 'Part ', child_idx, ':'

      ! Print signature using direct array access
      if (array_trees%parts(child_idx)%signature_length > 0) then
         write(stderr, '(A)', advance='no') ' ['
         do i = 1, array_trees%parts(child_idx)%signature_length
            if (i > 1) write(stderr, '(A)', advance='no') ', '
            write(stderr, '(I0)', advance='no') array_trees%parts(child_idx)%signature(i)
         end do
         write(stderr, '(A)') ']'
      else
         write(stderr, '(A)') ' []'
      end if

      ! Recursively print this child's children
      call print_signatures_recursive_array(array_trees, child_idx)

      child_idx = array_trees%parts(child_idx)%next_sibling_idx
   end do
end subroutine

subroutine print_leaf_items_array(array_trees)
   type(array_trees_t), intent(in) :: array_trees

   if (array_trees%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "   Leaf Items (Array)"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Print leaf items recursively (root is always at index 1)
   call print_leaf_items_recursive_array(array_trees, 1)

   write(stderr, *)
end subroutine

recursive subroutine print_leaf_items_recursive_array(array_trees, part_idx)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   integer :: child_idx

   if (part_idx == 0) return

   ! Process all children
   child_idx = array_trees%parts(part_idx)%first_child_idx
   do while (child_idx > 0)
      ! Only print items if this is a leaf part (no children)
      if (array_trees%parts(child_idx)%num_children == 0) then
         write(stderr, '(A,I0,A)', advance='no') 'Part ', child_idx, ':'
         call print_part_items_array(array_trees, child_idx)
      end if

      ! Recursively traverse this child's children to find more leaves
      call print_leaf_items_recursive_array(array_trees, child_idx)

      child_idx = array_trees%parts(child_idx)%next_sibling_idx
   end do
end subroutine

subroutine print_chain_details_array(array_trees)
   type(array_trees_t), intent(in) :: array_trees
   integer :: i, link_idx, j, k, part_idx

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 40)
   write(stderr, '(A)') "        Chain Details (Array)"
   write(stderr, '(A)') repeat("=", 40)
   write(stderr, *)

   do i = 1, array_trees%total_chains
      write(stderr, '(A,I0,A,I0,A,I0,A)') 'Chain ', i, ': ', &
         array_trees%chains(i)%num_links, ' links, ', &
         array_trees%chains(i)%num_children, ' children'

      if (array_trees%chains(i)%split_part_idx > 0) then
         write(stderr, '(A,I0)') '  Split part: ', array_trees%chains(i)%split_part_idx
      end if

      ! Show links in this chain using pure array access
      do j = 1, array_trees%chains(i)%num_links
         link_idx = array_trees%chains(i)%link_start_idx + j - 1
         write(stderr, '(A,I0,A,I0,A)', advance='no') '  Link ', link_idx, &
            ' (', array_trees%links(link_idx)%num_parts, ' parts): '

         ! Show parts in this link using pure array access
         do k = 1, array_trees%links(link_idx)%num_parts
            part_idx = array_trees%partref_entries(array_trees%links(link_idx)%partref_start_idx + k - 1)
            write(stderr, '(I0)', advance='no') part_idx
            if (k < array_trees%links(link_idx)%num_parts) write(stderr, '(A)', advance='no') ', '
         end do
         write(stderr, *)
      end do

      if (i < array_trees%total_chains) write(stderr, *)
   end do

   write(stderr, *)
end subroutine

subroutine print_first_level_items_array(array_trees)
   type(array_trees_t), intent(in) :: array_trees
   integer :: child_idx

   if (array_trees%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "  First Level Items (Array)"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Get first child of root (root is always at index 1)
   child_idx = array_trees%parts(1)%first_child_idx

   ! Print items for all parts at first level (direct children of root)
   do while (child_idx > 0)
      write(stderr, '(A,I0,A)', advance='no') 'Part ', child_idx, ':'
      call print_part_items_array(array_trees, child_idx)
      child_idx = array_trees%parts(child_idx)%next_sibling_idx
   end do

   write(stderr, *)
end subroutine

end module
