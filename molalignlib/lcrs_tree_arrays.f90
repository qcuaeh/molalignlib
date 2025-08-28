module lcrs_tree_arrays
use parameters
use lcrs_tree
use molecule
implicit none
private
public convert_trees_to_arrays
public validate_conversion
public print_tree_items_array
public print_part_tree_array
public print_leaf_items_array
public print_first_level_items_array
public print_chain_tree_array
public print_part_signatures_array
public print_chain_details_array

type, public :: partree_item_t
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

type, public :: chain_item_t
   integer :: num_parts
   integer :: parent_chain_idx    ! which chain owns this link
   ! Part reference segment in flattened array (using offset)
   integer :: partref_offset      ! offset into partref_entries array
   ! NOTE: itemdir1_offset and itemdir2_offset REMOVED - no longer needed with 2D arrays
end type

type, public :: assigntree_item_t
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

! Array-based assignment tree with 2D itemdir arrays
type, public :: array_trees_t
   ! Pure arrays for item values (no linked lists!)
   integer, allocatable :: vertidcs1(:)
   integer, allocatable :: vertidcs2(:)
   type(chain_item_t), allocatable :: chain(:)
   type(partree_item_t), allocatable :: partree(:)
   type(assigntree_item_t), allocatable :: assigntree(:)
   ! Flattened variable-length data - all pure integer arrays!
   integer, allocatable :: itemdir1_entries(:,:)  ! [link_idx, atom_idx] - vertices1 itemdir 2D array
   integer, allocatable :: itemdir2_entries(:,:)  ! [link_idx, atom_idx] - vertices2 itemdir 2D array
   integer, allocatable :: partref_entries(:)     ! Part indices for partrefs
   ! Adjacency information stored directly for fastest access
   integer, allocatable :: adj_lists1(:,:)     ! Direct 2D adjacency lists for vertices1 [atom_idx, neighbor_idx]
   integer, allocatable :: adj_lists2(:,:)     ! Direct 2D adjacency lists for vertices2 [atom_idx, neighbor_idx]
   integer, allocatable :: adj_counts1(:)      ! Count for each vertices1 atom's adjacency list
   integer, allocatable :: adj_counts2(:)      ! Count for each vertices2 atom's adjacency list
   ! Metadata
   integer :: num_atoms1, num_atoms2  ! number of atoms in each molecule
   integer :: total_items1, total_items2, total_parts
   integer :: total_links, total_chains
   integer :: total_partref_entries
   ! Assignment statistics
   real(rk) :: local_combinations
   real(rk) :: global_combinations
end type

contains

subroutine convert_trees_to_arrays(vertices1, vertices2, part_tree, assign_tree, assign_arrays)
   type(vertex_t), dimension(:), intent(in) :: vertices1, vertices2
   type(partree_node_t), pointer, intent(in) :: part_tree
   type(assigntree_node_t), pointer, intent(in) :: assign_tree
   type(array_trees_t), intent(out) :: assign_arrays
   ! Local variables
   integer :: partref_idx, link_idx, item1_idx, item2_idx

   ! Get totals from the tree counters
   assign_arrays%total_parts = part_tree%total_parts
   assign_arrays%total_items1 = part_tree%total_items1
   assign_arrays%total_items2 = part_tree%total_items2
   assign_arrays%total_chains = assign_tree%total_chains
   assign_arrays%total_links = assign_tree%total_links
   assign_arrays%total_partref_entries = assign_tree%total_partrefs

   ! Store molecule sizes
   assign_arrays%num_atoms1 = size(vertices1)
   assign_arrays%num_atoms2 = size(vertices2)

   ! Allocate all arrays with exact sizes (existing allocation code)
   allocate(assign_arrays%vertidcs1(assign_arrays%total_items1))
   allocate(assign_arrays%vertidcs2(assign_arrays%total_items2))
   allocate(assign_arrays%partree(assign_arrays%total_parts))
   allocate(assign_arrays%chain(assign_arrays%total_links))
   allocate(assign_arrays%assigntree(assign_arrays%total_chains))
   allocate(assign_arrays%partref_entries(assign_arrays%total_partref_entries))

   ! NEW: Allocate 2D itemdir arrays - one for each molecule
   allocate(assign_arrays%itemdir1_entries(assign_arrays%total_links, assign_arrays%num_atoms1))
   allocate(assign_arrays%itemdir2_entries(assign_arrays%total_links, assign_arrays%num_atoms2))

   ! Allocate adjacency arrays - 2D format for direct access using MAX_COORD
   allocate(assign_arrays%adj_lists1(assign_arrays%num_atoms1, MAX_COORD))
   allocate(assign_arrays%adj_lists2(assign_arrays%num_atoms2, MAX_COORD))
   allocate(assign_arrays%adj_counts1(assign_arrays%num_atoms1))
   allocate(assign_arrays%adj_counts2(assign_arrays%num_atoms2))

   ! OPTIMIZATION: Use intrinsic array operations instead of explicit loops
   assign_arrays%partref_entries = 0
   assign_arrays%vertidcs1 = 0
   assign_arrays%vertidcs2 = 0

   ! NEW: Initialize 2D itemdir arrays
   assign_arrays%itemdir1_entries = 0
   assign_arrays%itemdir2_entries = 0

   ! Initialize adjacency arrays
   assign_arrays%adj_lists1 = 0
   assign_arrays%adj_lists2 = 0
   assign_arrays%adj_counts1 = 0
   assign_arrays%adj_counts2 = 0

   ! Populate adjacency information
   call populate_adjacency_arrays(vertices1, vertices2, assign_arrays)

   ! Convert part tree starting from root with global item tracking
   item1_idx = 0
   item2_idx = 0
   call convert_parts_recurse(part_tree, assign_arrays, item1_idx, item2_idx)

   ! Convert chain tree and count statistics in single traversal
   partref_idx = 0
   link_idx = 0
   call convert_chains_recurse(assign_tree, assign_arrays, partref_idx, link_idx, &
                               assign_arrays%local_combinations, assign_arrays%global_combinations)
end subroutine

subroutine populate_adjacency_arrays(vertices1, vertices2, assign_arrays)
   ! Populate the 2D adjacency arrays directly - each atom gets its own row
   type(vertex_t), dimension(:), intent(in) :: vertices1, vertices2
   type(array_trees_t), intent(inout) :: assign_arrays
   integer :: i, j

   ! Populate vertices1 adjacency information - direct 2D storage
   do i = 1, assign_arrays%num_atoms1
      assign_arrays%adj_counts1(i) = size(vertices1(i)%adjlist)

      ! Copy adjacency list directly to 2D array
      do j = 1, size(vertices1(i)%adjlist)
         assign_arrays%adj_lists1(i, j) = vertices1(i)%adjlist(j)
      end do

      ! Zero out unused entries (though not strictly necessary)
      do j = size(vertices1(i)%adjlist) + 1, MAX_COORD
         assign_arrays%adj_lists1(i, j) = 0
      end do
   end do

   ! Populate vertices2 adjacency information - direct 2D storage
   do i = 1, assign_arrays%num_atoms2
      assign_arrays%adj_counts2(i) = size(vertices2(i)%adjlist)

      ! Copy adjacency list directly to 2D array
      do j = 1, size(vertices2(i)%adjlist)
         assign_arrays%adj_lists2(i, j) = vertices2(i)%adjlist(j)
      end do

      ! Zero out unused entries (though not strictly necessary)
      do j = size(vertices2(i)%adjlist) + 1, MAX_COORD
         assign_arrays%adj_lists2(i, j) = 0
      end do
   end do
end subroutine

subroutine convert_signature(part, assign_arrays, part_idx)
   type(partree_node_t), pointer, intent(in) :: part
   type(array_trees_t), intent(inout) :: assign_arrays
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
   assign_arrays%partree(part_idx)%signature_length = temp_count

   ! Second pass: compute unique values and their frequencies
   assign_arrays%partree(part_idx)%signature_unique_count = 0
   do i = 1, temp_count
      value = temp_values(i)
      found = .false.

      ! Check if this value is already in unique list
      do j = 1, assign_arrays%partree(part_idx)%signature_unique_count
         if (assign_arrays%partree(part_idx)%signature_values(j) == value) then
            assign_arrays%partree(part_idx)%signature_frequencies(j) = &
               assign_arrays%partree(part_idx)%signature_frequencies(j) + 1
            found = .true.
            exit
         end if
      end do

      ! If not found, add as new unique value
      if (.not. found) then
         assign_arrays%partree(part_idx)%signature_unique_count = &
            assign_arrays%partree(part_idx)%signature_unique_count + 1
         assign_arrays%partree(part_idx)%signature_values(assign_arrays%partree(part_idx)%signature_unique_count) = value
         assign_arrays%partree(part_idx)%signature_frequencies(assign_arrays%partree(part_idx)%signature_unique_count) = 1
      end if
   end do

   ! Zero out unused entries using intrinsic operation
   assign_arrays%partree(part_idx)%signature_values(assign_arrays%partree(part_idx)%signature_unique_count + 1:MAX_COORD) = 0
   assign_arrays%partree(part_idx)%signature_frequencies(assign_arrays%partree(part_idx)%signature_unique_count + 1:MAX_COORD) = 0
end subroutine

recursive subroutine convert_parts_recurse(part, assign_arrays, item1_idx, item2_idx)
   type(partree_node_t), pointer, intent(in) :: part
   type(array_trees_t), intent(inout) :: assign_arrays
   integer, intent(inout) :: item1_idx, item2_idx
   type(partree_node_t), pointer :: child_part
   type(item_node_t), pointer :: item
   integer :: part_idx, i, child_count

   if (.not. associated(part)) return

   ! Convert this part (global indices always start at 1)
   part_idx = part%global_index

   assign_arrays%partree(part_idx)%depth = part%depth
   assign_arrays%partree(part_idx)%num_children = part%num_children

   ! Relationships using global indices
   if (associated(part%parent_part)) then
      assign_arrays%partree(part_idx)%parent_part_idx = part%parent_part%global_index
   else
      assign_arrays%partree(part_idx)%parent_part_idx = 0
   end if

   if (associated(part%first_child_part)) then
      assign_arrays%partree(part_idx)%first_child_idx = part%first_child_part%global_index
   else
      assign_arrays%partree(part_idx)%first_child_idx = 0
   end if

   if (associated(part%last_child_part)) then
      assign_arrays%partree(part_idx)%last_child_idx = part%last_child_part%global_index
   else
      assign_arrays%partree(part_idx)%last_child_idx = 0
   end if

   if (associated(part%next_sibling_part)) then
      assign_arrays%partree(part_idx)%next_sibling_idx = part%next_sibling_part%global_index
   else
      assign_arrays%partree(part_idx)%next_sibling_idx = 0
   end if

   ! Set up item segments using offset approach (offset = start_idx - 1)
   assign_arrays%partree(part_idx)%items1_offset = item1_idx  ! item1_idx tracks the last used index
   assign_arrays%partree(part_idx)%items1_count = part%num_items1

   assign_arrays%partree(part_idx)%items2_offset = item2_idx  ! item2_idx tracks the last used index
   assign_arrays%partree(part_idx)%items2_count = part%num_items2

   ! OPTIMIZED: Convert signature to unique values, frequencies, and total length
   call convert_signature(part, assign_arrays, part_idx)

   ! Convert items1 to pure array format
   item => part%first_item1
   i = 0
   do while (associated(item))
      i = i + 1
      item1_idx = item1_idx + 1
      assign_arrays%vertidcs1(item1_idx) = item%vertidx
      item => item%next_item
   end do

   ! Convert items2 to pure array format
   item => part%first_item2
   i = 0
   do while (associated(item))
      i = i + 1
      item2_idx = item2_idx + 1
      assign_arrays%vertidcs2(item2_idx) = item%vertidx
      item => item%next_item
   end do

   ! OPTIMIZATION: Populate direct child access array for faster traversal
   if (part%num_children > 0) then
      allocate(assign_arrays%partree(part_idx)%child_indices(part%num_children))
      child_part => part%first_child_part
      child_count = 0
      do while (associated(child_part))
         child_count = child_count + 1
         assign_arrays%partree(part_idx)%child_indices(child_count) = child_part%global_index
         child_part => child_part%next_sibling_part
      end do
   end if

   ! Recursively convert all children
   child_part => part%first_child_part
   do while (associated(child_part))
      call convert_parts_recurse(child_part, assign_arrays, item1_idx, item2_idx)
      child_part => child_part%next_sibling_part
   end do
end subroutine

recursive subroutine convert_chains_recurse(chain, assign_arrays, partref_idx, link_idx, &
                                            local_combinations, global_combinations)
   type(assigntree_node_t), pointer, intent(in) :: chain
   type(array_trees_t), intent(inout) :: assign_arrays
   integer, intent(inout) :: partref_idx, link_idx
   real(rk), intent(out) :: local_combinations, global_combinations
   type(assigntree_node_t), pointer :: child_chain
   type(chain_node_t), pointer :: link
   type(partref_node_t), pointer :: partref
   integer :: chain_idx, current_link_idx, child_count
   real(rk) :: child_combinations, child_product
   integer :: items2_count

   if (.not. associated(chain)) return

   ! If this is a leaf level (no children), both values are 1
   if (chain%num_children == 0) then
      local_combinations = 1
      global_combinations = 1
   else
      ! Initialize accumulators for non-leaf nodes
      local_combinations = 0
      global_combinations = 1
   end if

   ! Convert this chain (existing conversion logic)
   chain_idx = chain%global_index

   assign_arrays%assigntree(chain_idx)%tot_items1 = chain%tot_items1
   assign_arrays%assigntree(chain_idx)%tot_items2 = chain%tot_items2
   assign_arrays%assigntree(chain_idx)%num_links = chain%num_links
   assign_arrays%assigntree(chain_idx)%num_children = chain%num_children

   ! Cross-tree reference
   if (associated(chain%split_part)) then
      assign_arrays%assigntree(chain_idx)%split_part_idx = chain%split_part%global_index
   else
      assign_arrays%assigntree(chain_idx)%split_part_idx = 0
   end if

   ! Chain relationships
   if (associated(chain%parent_chain)) then
      assign_arrays%assigntree(chain_idx)%parent_chain_idx = chain%parent_chain%global_index
   else
      assign_arrays%assigntree(chain_idx)%parent_chain_idx = 0
   end if

   if (associated(chain%first_child_chain)) then
      assign_arrays%assigntree(chain_idx)%first_child_idx = chain%first_child_chain%global_index
   else
      assign_arrays%assigntree(chain_idx)%first_child_idx = 0
   end if

   if (associated(chain%last_child_chain)) then
      assign_arrays%assigntree(chain_idx)%last_child_idx = chain%last_child_chain%global_index
   else
      assign_arrays%assigntree(chain_idx)%last_child_idx = 0
   end if

   if (associated(chain%next_sibling_chain)) then
      assign_arrays%assigntree(chain_idx)%next_sibling_idx = chain%next_sibling_chain%global_index
   else
      assign_arrays%assigntree(chain_idx)%next_sibling_idx = 0
   end if

   ! OPTIMIZATION: Populate direct child access array for faster traversal
   if (chain%num_children > 0) then
      allocate(assign_arrays%assigntree(chain_idx)%child_indices(chain%num_children))
      child_chain => chain%first_child_chain
      child_count = 0
      do while (associated(child_chain))
         child_count = child_count + 1
         assign_arrays%assigntree(chain_idx)%child_indices(child_count) = child_chain%global_index
         child_chain => child_chain%next_sibling_chain
      end do
   end if

   ! Set link offset (offset = start_idx - 1)
   assign_arrays%assigntree(chain_idx)%link_offset = link_idx

   ! Convert links in this chain using offset approach
   link => chain%first_link
   do while (associated(link))
      link_idx = link_idx + 1
      current_link_idx = link_idx

      assign_arrays%chain(current_link_idx)%num_parts = link%num_parts
      assign_arrays%chain(current_link_idx)%parent_chain_idx = chain%global_index

      ! Set partref offset (offset = start_idx - 1)
      assign_arrays%chain(current_link_idx)%partref_offset = partref_idx

      ! Convert partrefs to pure array format
      partref => link%first_partref
      do while (associated(partref))
         partref_idx = partref_idx + 1
         assign_arrays%partref_entries(partref_idx) = partref%part%global_index
         partref => partref%nextref
      end do

      link => link%next_link
   end do

   ! Recursively convert child chains and accumulate statistics
   child_chain => chain%first_child_chain
   do while (associated(child_chain))
      call convert_chains_recurse(child_chain, assign_arrays, partref_idx, link_idx, &
                                  child_combinations, child_product)

      ! Count statistics if this is not a leaf
      if (chain%num_children > 0) then
         items2_count = child_chain%split_part%num_items2

         ! Update combinations (sum): first item1 with each item2
         local_combinations = local_combinations + (items2_count * child_combinations)

         ! Update product (multiply): split sizes only
         global_combinations = global_combinations * (items2_count * child_product)
      end if

      child_chain => child_chain%next_sibling_chain
   end do
end subroutine

subroutine validate_conversion(part_tree, assign_tree, assign_arrays)
   type(partree_node_t), pointer, intent(in) :: part_tree
   type(assigntree_node_t), pointer, intent(in) :: assign_tree
   type(array_trees_t), intent(in) :: assign_arrays
   logical :: validation_passed

   validation_passed = .true.

   write(stderr, '(A)') "=== CONVERSION VALIDATION ==="

   ! Validate counts
   if (assign_arrays%total_parts /= part_tree%total_parts) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Part count mismatch: ", &
         assign_arrays%total_parts, " vs ", part_tree%total_parts
      validation_passed = .false.
   end if

   if (assign_arrays%total_items1 /= part_tree%total_items1) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Items1 count mismatch: ", &
         assign_arrays%total_items1, " vs ", part_tree%total_items1
      validation_passed = .false.
   end if

   if (assign_arrays%total_items2 /= part_tree%total_items2) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Items2 count mismatch: ", &
         assign_arrays%total_items2, " vs ", part_tree%total_items2
      validation_passed = .false.
   end if

   if (assign_arrays%total_chains /= assign_tree%total_chains) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Chain count mismatch: ", &
         assign_arrays%total_chains, " vs ", assign_tree%total_chains
      validation_passed = .false.
   end if

   if (assign_arrays%total_links /= assign_tree%total_links) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Link count mismatch: ", &
         assign_arrays%total_links, " vs ", assign_tree%total_links
      validation_passed = .false.
   end if

   if (assign_arrays%total_partref_entries /= assign_tree%total_partrefs) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Partref count mismatch: ", &
         assign_arrays%total_partref_entries, " vs ", assign_tree%total_partrefs
      validation_passed = .false.
   end if

   ! NEW: Validate 2D itemdir array dimensions
   if (size(assign_arrays%itemdir1_entries, 1) /= assign_arrays%total_links) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Itemdir1 links dimension mismatch: ", &
         size(assign_arrays%itemdir1_entries, 1), " vs ", assign_arrays%total_links
      validation_passed = .false.
   end if

   if (size(assign_arrays%itemdir1_entries, 2) /= assign_arrays%num_atoms1) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Itemdir1 atoms dimension mismatch: ", &
         size(assign_arrays%itemdir1_entries, 2), " vs ", assign_arrays%num_atoms1
      validation_passed = .false.
   end if

   if (size(assign_arrays%itemdir2_entries, 1) /= assign_arrays%total_links) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Itemdir2 links dimension mismatch: ", &
         size(assign_arrays%itemdir2_entries, 1), " vs ", assign_arrays%total_links
      validation_passed = .false.
   end if

   if (size(assign_arrays%itemdir2_entries, 2) /= assign_arrays%num_atoms2) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Itemdir2 atoms dimension mismatch: ", &
         size(assign_arrays%itemdir2_entries, 2), " vs ", assign_arrays%num_atoms2
      validation_passed = .false.
   end if

   if (validation_passed) then
      write(stderr, '(A)') "✓ Conversion validation PASSED"
   else
      write(stderr, '(A)') "✗ Conversion validation FAILED"
   end if
end subroutine

subroutine print_tree_items_array(assign_arrays)
   type(array_trees_t), intent(in) :: assign_arrays

   if (assign_arrays%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "      Part Items"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Root part is always at index 1, print its children recursively
   call print_items_recursive_array(assign_arrays, 1)

   write(stderr, *)
end subroutine

recursive subroutine print_items_recursive_array(assign_arrays, part_idx)
   type(array_trees_t), intent(in) :: assign_arrays
   integer, intent(in) :: part_idx
   integer :: child_idx, i

   ! Use direct array access instead of linked traversal for better performance
   do i = 1, assign_arrays%partree(part_idx)%num_children
      child_idx = assign_arrays%partree(part_idx)%child_indices(i)

      ! Print the child items with part index prefix
      write(stderr, '(A,I0,A)', advance='no') "Part ", child_idx, ':'
      call print_part_items_array(assign_arrays, child_idx)

      ! Recursively print this child's children
      call print_items_recursive_array(assign_arrays, child_idx)
   end do
end subroutine

subroutine print_part_items_array(assign_arrays, part_idx)
   type(array_trees_t), intent(in) :: assign_arrays
   integer, intent(in) :: part_idx
   integer :: i

   ! Print items1 using offset-based access
   do i = 1, assign_arrays%partree(part_idx)%items1_count
      write(stderr, '(1X,I0)', advance='no') assign_arrays%vertidcs1(assign_arrays%partree(part_idx)%items1_offset + i)
   end do

   write(stderr, '(A)', advance='no') ' /'

   ! Print items2 using offset-based access
   do i = 1, assign_arrays%partree(part_idx)%items2_count
      write(stderr, '(1X,I0)', advance='no') assign_arrays%vertidcs2(assign_arrays%partree(part_idx)%items2_offset + i)
   end do

   write(stderr, *)
end subroutine

! Array-based tree printing procedures

subroutine print_part_tree_array(assign_arrays)
   type(array_trees_t), intent(in) :: assign_arrays
   logical, dimension(:), allocatable :: is_last_child

   if (assign_arrays%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "    Part Tree"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Allocate tracking array for tree lines (max depth 100)
   allocate(is_last_child(100))
   is_last_child = .false.

   ! Print root line
   write(stderr, '(A)') 'ROOT'

   ! Print children recursively (root is always at index 1)
   call print_part_recursive_array(assign_arrays, 1, 0, is_last_child)

   deallocate(is_last_child)
   write(stderr, *)
end subroutine

recursive subroutine print_part_recursive_array(assign_arrays, part_idx, depth, is_last_child)
   type(array_trees_t), intent(in) :: assign_arrays
   integer, intent(in) :: part_idx, depth
   logical, dimension(:), intent(inout) :: is_last_child
   integer :: child_idx, i, j, pos
   character(len=200) :: prefix

   if (part_idx == 0) return

   ! Process all children using direct array access
   do i = 1, assign_arrays%partree(part_idx)%num_children
      child_idx = assign_arrays%partree(part_idx)%child_indices(i)

      ! Check if this is the last child
      is_last_child(depth + 1) = (i == assign_arrays%partree(part_idx)%num_children)

      ! Build prefix for this level
      prefix = ""
      pos = 1
      do j = 1, depth
         if (is_last_child(j)) then
            prefix(pos:pos+3) = "   "
         else
            prefix(pos:pos+3) = "|  "
         end if
         pos = pos + 3
      end do

      ! Add branch characters
      if (is_last_child(depth + 1)) then
         prefix(pos:pos+2) = "`--"
      else
         prefix(pos:pos+2) = "|--"
      end if
      pos = pos + 3

      ! Print part index with item counts
      write(stderr, '(A,A,I0,A,I0,A)') prefix(1:pos-1), '* (', &
         assign_arrays%partree(child_idx)%items1_count, '/', &
         assign_arrays%partree(child_idx)%items2_count, ')'

      ! Recursively print this child's children
      call print_part_recursive_array(assign_arrays, child_idx, depth + 1, is_last_child)
   end do
end subroutine

subroutine print_chain_tree_array(assign_arrays)
   type(array_trees_t), intent(in) :: assign_arrays
   logical, dimension(:), allocatable :: is_last_child

   if (assign_arrays%total_chains == 0) then
      write(stderr, '(A)') "Assignment tree is empty"
      return
   end if

   ! Allocate tracking array for tree lines (max depth 100)
   allocate(is_last_child(100))
   is_last_child = .false.

   write(stderr, *)
   write(stderr, '(A)') '*'

   ! Print children recursively (root is always at index 1)
   call print_chain_recursive_array(assign_arrays, 1, 0, is_last_child)
   write(stderr, *)

   ! Print assignment statistics
   if (assign_arrays%global_combinations > 2**24) then
      write(stderr, '(A,ES8.2)') "Global combinations: ", assign_arrays%global_combinations
   else
      write(stderr, '(A,I0)') "Global combinations: ", int(assign_arrays%global_combinations)
   end if
   if (assign_arrays%local_combinations > 2**24) then
      write(stderr, '(A,ES8.2)') "Total local combinations: ", assign_arrays%local_combinations
   else
      write(stderr, '(A,I0)') "Total local combinations: ", int(assign_arrays%local_combinations)
   end if

   deallocate(is_last_child)
end subroutine

recursive subroutine print_chain_recursive_array(assign_arrays, chain_idx, depth, is_last_child)
   type(array_trees_t), intent(in) :: assign_arrays
   integer, intent(in) :: chain_idx, depth
   logical, dimension(:), intent(inout) :: is_last_child
   integer :: child_idx, split_part_idx, i, j, pos
   character(len=200) :: prefix

   if (chain_idx == 0) return

   ! Process all children using direct array access instead of linked traversal
   do i = 1, assign_arrays%assigntree(chain_idx)%num_children
      child_idx = assign_arrays%assigntree(chain_idx)%child_indices(i)

      ! Check if this is the last child
      is_last_child(depth + 1) = (i == assign_arrays%assigntree(chain_idx)%num_children)

      ! Build prefix for this level
      prefix = ""
      pos = 1
      do j = 1, depth
         if (is_last_child(j)) then
            prefix(pos:pos+3) = "   "
         else
            prefix(pos:pos+3) = "|  "
         end if
         pos = pos + 3
      end do

      ! Add branch characters
      if (is_last_child(depth + 1)) then
         prefix(pos:pos+2) = "`--"
      else
         prefix(pos:pos+2) = "|--"
      end if
      pos = pos + 3

      ! Print the split part index with item counts
      split_part_idx = assign_arrays%assigntree(child_idx)%split_part_idx
      if (split_part_idx > 0) then
         write(stderr, '(A,A,I0,A,I0,A)') prefix(1:pos-1), '* (', &
            assign_arrays%partree(split_part_idx)%items1_count, '/', &
            assign_arrays%partree(split_part_idx)%items2_count, ')'
      else
         write(stderr, '(A,A)') prefix(1:pos-1), '(no split part)'
      end if

      ! Recursively print this child's children
      call print_chain_recursive_array(assign_arrays, child_idx, depth + 1, is_last_child)
   end do
end subroutine

subroutine print_part_signatures_array(assign_arrays)
   type(array_trees_t), intent(in) :: assign_arrays

   if (assign_arrays%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "  Part Signatures"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Print signatures for all parts (root is always at index 1)
   call print_signatures_recursive_array(assign_arrays, 1)

   write(stderr, *)
end subroutine

recursive subroutine print_signatures_recursive_array(assign_arrays, part_idx)
   type(array_trees_t), intent(in) :: assign_arrays
   integer, intent(in) :: part_idx
   integer :: child_idx, i, j

   if (part_idx == 0) return

   ! Process all children using direct array access
   do i = 1, assign_arrays%partree(part_idx)%num_children
      child_idx = assign_arrays%partree(part_idx)%child_indices(i)

      ! Print the child signature with frequencies and total length
      write(stderr,'(A,I0,A,I0,A)',advance='no') 'Part ', child_idx, ' (len=', &
         assign_arrays%partree(child_idx)%signature_length, '):'

      ! Print unique signature values with frequencies
      if (assign_arrays%partree(child_idx)%signature_unique_count > 0) then
         write(stderr, '(A)', advance='no') ' ['
         do j = 1, assign_arrays%partree(child_idx)%signature_unique_count
            if (j > 1) write(stderr, '(A)', advance='no') ', '
            write(stderr, '(I0,A,I0)', advance='no') &
               assign_arrays%partree(child_idx)%signature_values(j), '×', &
               assign_arrays%partree(child_idx)%signature_frequencies(j)
         end do
         write(stderr, '(A)') ']'
      else
         write(stderr, '(A)') ' []'
      end if

      ! Recursively print this child's children
      call print_signatures_recursive_array(assign_arrays, child_idx)
   end do
end subroutine

subroutine print_leaf_items_array(assign_arrays)
   type(array_trees_t), intent(in) :: assign_arrays

   if (assign_arrays%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "   Leaf Items"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Print leaf items recursively (root is always at index 1)
   call print_leaf_items_recursive_array(assign_arrays, 1)

   write(stderr, *)
end subroutine

recursive subroutine print_leaf_items_recursive_array(assign_arrays, part_idx)
   type(array_trees_t), intent(in) :: assign_arrays
   integer, intent(in) :: part_idx
   integer :: child_idx, i

   if (part_idx == 0) return

   ! Process all children using direct array access
   do i = 1, assign_arrays%partree(part_idx)%num_children
      child_idx = assign_arrays%partree(part_idx)%child_indices(i)

      ! Only print items if this is a leaf part (no children)
      if (assign_arrays%partree(child_idx)%num_children == 0) then
         write(stderr, '(A,I0,A)', advance='no') 'Part ', child_idx, ':'
         call print_part_items_array(assign_arrays, child_idx)
      end if

      ! Recursively traverse this child's children to find more leaves
      call print_leaf_items_recursive_array(assign_arrays, child_idx)
   end do
end subroutine

subroutine print_chain_details_array(assign_arrays)
   type(array_trees_t), intent(in) :: assign_arrays
   integer :: i, link_idx, j, k, part_idx

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 40)
   write(stderr, '(A)') "        Chain Details"
   write(stderr, '(A)') repeat("=", 40)
   write(stderr, *)

   do i = 1, assign_arrays%total_chains
      write(stderr, '(A,I0,A,I0,A,I0,A)') 'Chain ', i, ': ', &
         assign_arrays%assigntree(i)%num_links, ' links, ', &
         assign_arrays%assigntree(i)%num_children, ' children'

      if (assign_arrays%assigntree(i)%split_part_idx > 0) then
         write(stderr, '(A,I0)') '  Split part: ', assign_arrays%assigntree(i)%split_part_idx
      end if

      ! Show links in this chain using offset-based access
      do j = 1, assign_arrays%assigntree(i)%num_links
         link_idx = assign_arrays%assigntree(i)%link_offset + j
         write(stderr, '(A,I0,A,I0,A)', advance='no') '  Link ', link_idx, &
            ' (', assign_arrays%chain(link_idx)%num_parts, ' parts): '

         ! Show parts in this link using offset-based access
         do k = 1, assign_arrays%chain(link_idx)%num_parts
            part_idx = assign_arrays%partref_entries(assign_arrays%chain(link_idx)%partref_offset + k)
            write(stderr, '(I0)', advance='no') part_idx
            if (k < assign_arrays%chain(link_idx)%num_parts) write(stderr, '(A)', advance='no') ', '
         end do
         write(stderr, *)
      end do

      if (i < assign_arrays%total_chains) write(stderr, *)
   end do

   write(stderr, *)
end subroutine

subroutine print_first_level_items_array(assign_arrays)
   type(array_trees_t), intent(in) :: assign_arrays
   integer :: child_idx, i

   if (assign_arrays%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "  First Level Items"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Print items for all parts at first level using direct array access
   do i = 1, assign_arrays%partree(1)%num_children
      child_idx = assign_arrays%partree(1)%child_indices(i)
      write(stderr, '(A,I0,A)', advance='no') 'Part ', child_idx, ':'
      call print_part_items_array(assign_arrays, child_idx)
   end do

   write(stderr, *)
end subroutine

end module
