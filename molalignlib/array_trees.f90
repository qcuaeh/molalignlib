module array_trees
use parameters
use lcrs_tree
implicit none
private

! Arrays of derived types with scalar components
! REMOVED: item_array_t - replaced with pure arrays

type, public :: part_array_t
   integer :: depth
   integer :: num_children
   ! Relationships (0 = null)
   integer :: parent_part_idx
   integer :: first_child_idx
   integer :: last_child_idx
   integer :: next_sibling_idx
   ! OPTIMIZATION: Direct child access array - eliminates linked traversal
   integer, allocatable :: child_indices(:)  ! Direct array of child part indices
   ! Item segments in flattened arrays (using offsets)
   integer :: items1_offset, items1_count
   integer :: items2_offset, items2_count
   ! OPTIMIZED: Store unique signature values, frequencies, and total length
   integer :: signature_values(MAX_COORD)       ! unique values in signature
   integer :: signature_frequencies(MAX_COORD)  ! frequency of each unique value
   integer :: signature_unique_count            ! number of unique values
   integer :: signature_length            ! total signature length (sum of frequencies)
end type

! REMOVED: partref_array_t - replaced with pure array

type, public :: link_array_t
   integer :: num_parts
   integer :: parent_chain_idx    ! which chain owns this link
   ! Part reference segment in flattened array (using offset)
   integer :: partref_offset      ! offset into partref_entries array
   ! Item directories (flattened storage, using offsets)
   integer :: itemdir1_offset
   integer :: itemdir2_offset
end type

type, public :: chain_array_t
   integer :: tot_items1, tot_items2
   integer :: num_links, num_children
   ! Cross-tree reference (0 = null)
   integer :: split_part_idx      ! points to part array
   ! Chain tree relationships (0 = null)
   integer :: parent_chain_idx
   integer :: first_child_idx
   integer :: last_child_idx
   integer :: next_sibling_idx
   ! OPTIMIZATION: Direct child access array - eliminates linked traversal
   integer, allocatable :: child_indices(:)  ! Direct array of child chain indices
   ! Link segment in flattened array (using offset)
   integer :: link_offset         ! offset into links array
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

   ! OPTIMIZATION 4: Use intrinsic array operations instead of explicit loops
   ! These are highly optimized by the compiler and much faster than manual loops
   array_trees%itemdir_entries = 0      ! O(1) intrinsic vs O(n) explicit loop
   array_trees%partref_entries = 0      ! O(1) intrinsic vs O(n) explicit loop
   array_trees%item1_values = 0         ! O(1) intrinsic vs O(n) explicit loop
   array_trees%item2_values = 0         ! O(1) intrinsic vs O(n) explicit loop

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

! Modified signature conversion procedure
subroutine convert_signature(part, array_trees, part_idx)
   type(part_node_t), pointer, intent(in) :: part
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx
   integer :: temp_values(MAX_COORD)
   integer :: temp_count, i, j, value
   logical :: found

   ! First pass: collect all non-null signature values
   temp_count = 0
   do i = 1, size(part%signature)
      if (associated(part%signature(i)%ptr)) then
         temp_count = temp_count + 1
         temp_values(temp_count) = part%signature(i)%ptr%global_index
      end if
   end do

   ! Store total signature length
   array_trees%parts(part_idx)%signature_length = temp_count

   ! Second pass: compute unique values and their frequencies
   array_trees%parts(part_idx)%signature_unique_count = 0
   do i = 1, temp_count
      value = temp_values(i)
      found = .false.

      ! Check if this value is already in unique list
      do j = 1, array_trees%parts(part_idx)%signature_unique_count
         if (array_trees%parts(part_idx)%signature_values(j) == value) then
            array_trees%parts(part_idx)%signature_frequencies(j) = &
               array_trees%parts(part_idx)%signature_frequencies(j) + 1
            found = .true.
            exit
         end if
      end do

      ! If not found, add as new unique value
      if (.not. found) then
         array_trees%parts(part_idx)%signature_unique_count = &
            array_trees%parts(part_idx)%signature_unique_count + 1
         array_trees%parts(part_idx)%signature_values(array_trees%parts(part_idx)%signature_unique_count) = value
         array_trees%parts(part_idx)%signature_frequencies(array_trees%parts(part_idx)%signature_unique_count) = 1
      end if
   end do

   ! Zero out unused entries using intrinsic operation
   array_trees%parts(part_idx)%signature_values(array_trees%parts(part_idx)%signature_unique_count + 1:MAX_COORD) = 0
   array_trees%parts(part_idx)%signature_frequencies(array_trees%parts(part_idx)%signature_unique_count + 1:MAX_COORD) = 0
end subroutine

! Modified convert_parts_recursive - replace signature conversion section
recursive subroutine convert_parts_recursive(part, array_trees, item1_idx, item2_idx)
   type(part_node_t), pointer, intent(in) :: part
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(inout) :: item1_idx, item2_idx
   type(part_node_t), pointer :: child_part
   type(item_node_t), pointer :: item
   integer :: part_idx, i, child_count

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

   ! Set up item segments using offset approach (offset = start_idx - 1)
   array_trees%parts(part_idx)%items1_offset = item1_idx  ! item1_idx tracks the last used index
   array_trees%parts(part_idx)%items1_count = part%num_items1

   array_trees%parts(part_idx)%items2_offset = item2_idx  ! item2_idx tracks the last used index
   array_trees%parts(part_idx)%items2_count = part%num_items2

   ! OPTIMIZED: Convert signature to unique values, frequencies, and total length
   call convert_signature(part, array_trees, part_idx)

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

   ! OPTIMIZATION: Populate direct child access array for faster traversal
   if (part%num_children > 0) then
      allocate(array_trees%parts(part_idx)%child_indices(part%num_children))
      child_part => part%first_child_part
      child_count = 0
      do while (associated(child_part))
         child_count = child_count + 1
         array_trees%parts(part_idx)%child_indices(child_count) = child_part%global_index
         child_part => child_part%next_sibling_part
      end do
   end if

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
   integer :: chain_idx, current_link_idx, child_count

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

   ! OPTIMIZATION: Populate direct child access array for faster traversal
   if (chain%num_children > 0) then
      allocate(array_trees%chains(chain_idx)%child_indices(chain%num_children))
      child_chain => chain%first_child_chain
      child_count = 0
      do while (associated(child_chain))
         child_count = child_count + 1
         array_trees%chains(chain_idx)%child_indices(child_count) = child_chain%global_index
         child_chain => child_chain%next_sibling_chain
      end do
   end if

   ! Set link offset (offset = start_idx - 1)
   array_trees%chains(chain_idx)%link_offset = link_idx  ! link_idx tracks the last used index

   ! Convert links in this chain using offset approach
   link => chain%first_link
   do while (associated(link))
      link_idx = link_idx + 1
      current_link_idx = link_idx

      array_trees%links(current_link_idx)%num_parts = link%num_parts
      array_trees%links(current_link_idx)%parent_chain_idx = chain%global_index

      ! Set partref offset (offset = start_idx - 1)
      array_trees%links(current_link_idx)%partref_offset = partref_idx  ! partref_idx tracks the last used index

      ! Convert partrefs to pure array format
      partref => link%first_partref
      do while (associated(partref))
         partref_idx = partref_idx + 1
         array_trees%partref_entries(partref_idx) = partref%part%global_index
         partref => partref%nextref
      end do

      ! Set itemdir offsets (offset = start_idx - 1)
      array_trees%links(current_link_idx)%itemdir1_offset = itemdir_idx  ! itemdir_idx tracks the last used index
      itemdir_idx = itemdir_idx + array_trees%itemdir_size1

      array_trees%links(current_link_idx)%itemdir2_offset = itemdir_idx  ! itemdir_idx tracks the last used index
      itemdir_idx = itemdir_idx + array_trees%itemdir_size2

      link => link%next_link
   end do

   ! Recursively convert child chains using direct array access
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
   integer :: child_idx, i

   ! Use direct array access instead of linked traversal for better performance
   do i = 1, array_trees%parts(part_idx)%num_children
      child_idx = array_trees%parts(part_idx)%child_indices(i)

      ! Print the child items with part index prefix
      write(stderr, '(A,I0,A)', advance='no') "Part ", child_idx, ':'
      call print_part_items_array(array_trees, child_idx)

      ! Recursively print this child's children
      call print_items_recursive_array(array_trees, child_idx)
   end do
end subroutine

subroutine print_part_items_array(array_trees, part_idx)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   integer :: i

   ! Print items1 using offset-based access
   do i = 1, array_trees%parts(part_idx)%items1_count
      write(stderr, '(1X,I0)', advance='no') array_trees%item1_values(array_trees%parts(part_idx)%items1_offset + i)
   end do

   write(stderr, '(A)', advance='no') ' /'

   ! Print items2 using offset-based access
   do i = 1, array_trees%parts(part_idx)%items2_count
      write(stderr, '(1X,I0)', advance='no') array_trees%item2_values(array_trees%parts(part_idx)%items2_offset + i)
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
   integer :: child_idx, i, j, pos
   character(len=200) :: prefix

   if (part_idx == 0) return

   ! Process all children using direct array access
   do i = 1, array_trees%parts(part_idx)%num_children
      child_idx = array_trees%parts(part_idx)%child_indices(i)

      ! Check if this is the last child
      is_last_child(depth + 1) = (i == array_trees%parts(part_idx)%num_children)

      ! Build prefix for this level
      prefix = " "
      pos = 2
      do j = 1, depth
         if (is_last_child(j)) then
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
   integer :: child_idx, split_part_idx, i, j, pos
   character(len=200) :: prefix

   if (chain_idx == 0) return

   ! Process all children using direct array access instead of linked traversal
   do i = 1, array_trees%chains(chain_idx)%num_children
      child_idx = array_trees%chains(chain_idx)%child_indices(i)

      ! Check if this is the last child
      is_last_child(depth + 1) = (i == array_trees%chains(chain_idx)%num_children)

      ! Build prefix for this level
      prefix = " "
      pos = 2
      do j = 1, depth
         if (is_last_child(j)) then
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
   integer :: child_idx, i, j

   if (part_idx == 0) return

   ! Process all children using direct array access
   do i = 1, array_trees%parts(part_idx)%num_children
      child_idx = array_trees%parts(part_idx)%child_indices(i)

      ! Print the child signature with frequencies and total length
      write(stderr,'(A,I0,A,I0,A)',advance='no') 'Part ', child_idx, ' (len=', &
         array_trees%parts(child_idx)%signature_length, '):'

      ! Print unique signature values with frequencies
      if (array_trees%parts(child_idx)%signature_unique_count > 0) then
         write(stderr, '(A)', advance='no') ' ['
         do j = 1, array_trees%parts(child_idx)%signature_unique_count
            if (j > 1) write(stderr, '(A)', advance='no') ', '
            write(stderr, '(I0,A,I0)', advance='no') &
               array_trees%parts(child_idx)%signature_values(j), '×', &
               array_trees%parts(child_idx)%signature_frequencies(j)
         end do
         write(stderr, '(A)') ']'
      else
         write(stderr, '(A)') ' []'
      end if

      ! Recursively print this child's children
      call print_signatures_recursive_array(array_trees, child_idx)
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
   integer :: child_idx, i

   if (part_idx == 0) return

   ! Process all children using direct array access
   do i = 1, array_trees%parts(part_idx)%num_children
      child_idx = array_trees%parts(part_idx)%child_indices(i)

      ! Only print items if this is a leaf part (no children)
      if (array_trees%parts(child_idx)%num_children == 0) then
         write(stderr, '(A,I0,A)', advance='no') 'Part ', child_idx, ':'
         call print_part_items_array(array_trees, child_idx)
      end if

      ! Recursively traverse this child's children to find more leaves
      call print_leaf_items_recursive_array(array_trees, child_idx)
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

      ! Show links in this chain using offset-based access
      do j = 1, array_trees%chains(i)%num_links
         link_idx = array_trees%chains(i)%link_offset + j
         write(stderr, '(A,I0,A,I0,A)', advance='no') '  Link ', link_idx, &
            ' (', array_trees%links(link_idx)%num_parts, ' parts): '

         ! Show parts in this link using offset-based access
         do k = 1, array_trees%links(link_idx)%num_parts
            part_idx = array_trees%partref_entries(array_trees%links(link_idx)%partref_offset + k)
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
   integer :: child_idx, i

   if (array_trees%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "  First Level Items (Array)"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Print items for all parts at first level using direct array access
   do i = 1, array_trees%parts(1)%num_children
      child_idx = array_trees%parts(1)%child_indices(i)
      write(stderr, '(A,I0,A)', advance='no') 'Part ', child_idx, ':'
      call print_part_items_array(array_trees, child_idx)
   end do

   write(stderr, *)
end subroutine

end module
