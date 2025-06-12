module array_trees
use parameters
use lcrs_tree
implicit none

! Arrays of derived types with scalar components
type :: item_array_t
   integer :: value
   integer :: next_item_idx  ! 0 = null
   integer :: part_idx       ! which part owns this item
end type

type :: part_array_t
   integer :: depth
   integer :: num_items1, num_items2, num_children
   ! Relationships (0 = null)
   integer :: parent_part_idx
   integer :: first_child_idx
   integer :: last_child_idx
   integer :: next_sibling_idx
   ! Item lists (0 = null)
   integer :: first_item1_idx
   integer :: first_item2_idx
   integer :: last_item1_idx
   integer :: last_item2_idx
   ! Fixed-size signature (much simpler!)
   integer :: signature(MAX_COORD)
   integer :: signature_length
end type

type :: partref_array_t
   integer :: part_idx        ! which part is referenced
   integer :: next_ref_idx    ! 0 = null
   integer :: parent_link_idx ! which link owns this ref
end type

type :: link_array_t
   integer :: num_parts
   integer :: next_link_idx       ! 0 = null
   integer :: parent_chain_idx    ! which chain owns this link
   ! Part reference lists (0 = null)
   integer :: first_partref_idx
   integer :: last_partref_idx
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
   ! Link lists (0 = null)
   integer :: first_link_idx
   integer :: last_link_idx
end type

! Flattened arrays for variable-length data
type :: itemdir_entry_t
   integer :: part_idx  ! part index for this item in itemdir (0 = null)
end type

! Complete array-based representation
type, public :: array_trees_t
   ! Separate items1/items2 arrays
   type(item_array_t), allocatable :: items1(:)
   type(item_array_t), allocatable :: items2(:)
   type(part_array_t), allocatable :: parts(:)
   type(partref_array_t), allocatable :: partrefs(:)
   type(link_array_t), allocatable :: links(:)
   type(chain_array_t), allocatable :: chains(:)

   ! Flattened variable-length data (only itemdirs now)
   type(itemdir_entry_t), allocatable :: itemdir_entries(:)

   ! Metadata
   integer :: total_items1, total_items2, total_parts
   integer :: total_partrefs, total_links, total_chains
   integer :: total_itemdir_entries
   integer :: itemdir_size1, itemdir_size2  ! size of each itemdir
end type

public :: convert_trees_to_arrays
public :: validate_conversion
public :: print_array_summary
public :: print_tree_items_array

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
   array_trees%total_partrefs = root_chain%total_partrefs
   array_trees%itemdir_size1 = root_chain%tot_items1
   array_trees%itemdir_size2 = root_chain%tot_items2

   ! Calculate itemdir storage needs only
   array_trees%total_itemdir_entries = array_trees%total_links * &
                                      (array_trees%itemdir_size1 + array_trees%itemdir_size2)

   ! Allocate all arrays with exact sizes (no signature_entries!)
   allocate(array_trees%items1(array_trees%total_items1))
   allocate(array_trees%items2(array_trees%total_items2))
   allocate(array_trees%parts(array_trees%total_parts))
   allocate(array_trees%partrefs(array_trees%total_partrefs))
   allocate(array_trees%links(array_trees%total_links))
   allocate(array_trees%chains(array_trees%total_chains))
   allocate(array_trees%itemdir_entries(array_trees%total_itemdir_entries))

   ! Initialize itemdir entries to 0 (will be recomputed in mna_recompute_arrays)
   do i = 1, array_trees%total_itemdir_entries
      array_trees%itemdir_entries(i)%part_idx = 0
   end do

   ! Convert using global indices (much simpler now!)
   call populate_arrays_direct(root_part, root_chain, array_trees)
end subroutine

subroutine populate_arrays_direct(root_part, root_chain, array_trees)
   type(part_node_t), pointer, intent(in) :: root_part
   type(chain_node_t), pointer, intent(in) :: root_chain
   type(array_trees_t), intent(inout) :: array_trees
   integer :: itemdir_idx  ! Add global itemdir index tracker

   ! Convert part tree starting from root (no changes needed)
   call convert_parts_recursive(root_part, array_trees)

   ! Convert chain tree starting from root with global itemdir tracking
   itemdir_idx = 0
   call convert_chains_recursive(root_chain, array_trees, itemdir_idx)
end subroutine

recursive subroutine convert_parts_recursive(part, array_trees)
   type(part_node_t), pointer, intent(in) :: part
   type(array_trees_t), intent(inout) :: array_trees
   type(part_node_t), pointer :: child_part
   type(item_node_t), pointer :: item
   integer :: part_idx, i

   if (.not. associated(part)) return

   ! Convert this part (global indices always start at 1)
   part_idx = part%global_index

   array_trees%parts(part_idx)%depth = part%depth
   array_trees%parts(part_idx)%num_items1 = part%num_items1
   array_trees%parts(part_idx)%num_items2 = part%num_items2
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

   ! Item lists using global indices
   if (associated(part%first_item1)) then
      array_trees%parts(part_idx)%first_item1_idx = part%first_item1%global_index
   else
      array_trees%parts(part_idx)%first_item1_idx = 0
   end if

   if (associated(part%first_item2)) then
      array_trees%parts(part_idx)%first_item2_idx = part%first_item2%global_index
   else
      array_trees%parts(part_idx)%first_item2_idx = 0
   end if

   if (associated(part%last_item1)) then
      array_trees%parts(part_idx)%last_item1_idx = part%last_item1%global_index
   else
      array_trees%parts(part_idx)%last_item1_idx = 0
   end if

   if (associated(part%last_item2)) then
      array_trees%parts(part_idx)%last_item2_idx = part%last_item2%global_index
   else
      array_trees%parts(part_idx)%last_item2_idx = 0
   end if

   ! Convert signature to fixed-size array (much simpler!)
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

   ! Convert items1
   item => part%first_item1
   do while (associated(item))
      array_trees%items1(item%global_index)%value = item%value
      array_trees%items1(item%global_index)%part_idx = part_idx
      if (associated(item%next_item)) then
         array_trees%items1(item%global_index)%next_item_idx = item%next_item%global_index
      else
         array_trees%items1(item%global_index)%next_item_idx = 0
      end if
      item => item%next_item
   end do

   ! Convert items2
   item => part%first_item2
   do while (associated(item))
      array_trees%items2(item%global_index)%value = item%value
      array_trees%items2(item%global_index)%part_idx = part_idx
      if (associated(item%next_item)) then
         array_trees%items2(item%global_index)%next_item_idx = item%next_item%global_index
      else
         array_trees%items2(item%global_index)%next_item_idx = 0
      end if
      item => item%next_item
   end do

   ! Recursively convert all children
   child_part => part%first_child_part
   do while (associated(child_part))
      call convert_parts_recursive(child_part, array_trees)
      child_part => child_part%next_sibling_part
   end do
end subroutine

recursive subroutine convert_chains_recursive(chain, array_trees, itemdir_idx)
   type(chain_node_t), pointer, intent(in) :: chain
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(inout) :: itemdir_idx  ! Change to inout for global tracking
   type(chain_node_t), pointer :: child_chain
   type(link_node_t), pointer :: link
   type(partref_node_t), pointer :: partref
   integer :: chain_idx, link_idx

   if (.not. associated(chain)) return

   ! Convert this chain (global indices always start at 1)
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

   ! Link lists
   if (associated(chain%first_link)) then
      array_trees%chains(chain_idx)%first_link_idx = chain%first_link%global_index
   else
      array_trees%chains(chain_idx)%first_link_idx = 0
   end if

   if (associated(chain%last_link)) then
      array_trees%chains(chain_idx)%last_link_idx = chain%last_link%global_index
   else
      array_trees%chains(chain_idx)%last_link_idx = 0
   end if

   ! Convert links in this chain using GLOBAL itemdir tracking
   link => chain%first_link
   do while (associated(link))
      link_idx = link%global_index

      array_trees%links(link_idx)%num_parts = link%num_parts
      if (associated(link%next_link)) then
         array_trees%links(link_idx)%next_link_idx = link%next_link%global_index
      else
         array_trees%links(link_idx)%next_link_idx = 0
      end if

      array_trees%links(link_idx)%parent_chain_idx = chain%global_index

      ! Partref lists
      if (associated(link%first_partref)) then
         array_trees%links(link_idx)%first_partref_idx = link%first_partref%global_index
      else
         array_trees%links(link_idx)%first_partref_idx = 0
      end if

      if (associated(link%last_partref)) then
         array_trees%links(link_idx)%last_partref_idx = link%last_partref%global_index
      else
         array_trees%links(link_idx)%last_partref_idx = 0
      end if

      ! Convert partrefs
      partref => link%first_partref
      do while (associated(partref))
         array_trees%partrefs(partref%global_index)%part_idx = partref%part%global_index
         array_trees%partrefs(partref%global_index)%parent_link_idx = link_idx
         if (associated(partref%nextref)) then
            array_trees%partrefs(partref%global_index)%next_ref_idx = partref%nextref%global_index
         else
            array_trees%partrefs(partref%global_index)%next_ref_idx = 0
         end if
         partref => partref%nextref
      end do

      ! Set itemdir start indices using GLOBAL counter (fixed!)
      array_trees%links(link_idx)%itemdir1_start_idx = itemdir_idx + 1
      itemdir_idx = itemdir_idx + array_trees%itemdir_size1

      array_trees%links(link_idx)%itemdir2_start_idx = itemdir_idx + 1
      itemdir_idx = itemdir_idx + array_trees%itemdir_size2

      link => link%next_link
   end do

   ! Recursively convert child chains using GLOBAL counter
   child_chain => chain%first_child_chain
   do while (associated(child_chain))
      call convert_chains_recursive(child_chain, array_trees, itemdir_idx)
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

   if (array_trees%total_partrefs /= root_chain%total_partrefs) then
      write(stderr, '(A,I0,A,I0)') "ERROR: Partref count mismatch: ", &
         array_trees%total_partrefs, " vs ", root_chain%total_partrefs
      validation_passed = .false.
   end if

   if (validation_passed) then
      write(stderr, '(A)') "✓ Conversion validation PASSED"
   else
      write(stderr, '(A)') "✗ Conversion validation FAILED"
   end if
end subroutine

subroutine print_array_summary(array_trees)
   type(array_trees_t), intent(in) :: array_trees

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 50)
   write(stderr, '(A)') "           ARRAY CONVERSION SUMMARY"
   write(stderr, '(A)') repeat("=", 50)
   write(stderr, '(A,I0)') "Parts:              ", array_trees%total_parts
   write(stderr, '(A,I0)') "Items1:             ", array_trees%total_items1
   write(stderr, '(A,I0)') "Items2:             ", array_trees%total_items2
   write(stderr, '(A,I0)') "Chains:             ", array_trees%total_chains
   write(stderr, '(A,I0)') "Links:              ", array_trees%total_links
   write(stderr, '(A,I0)') "Partrefs:           ", array_trees%total_partrefs
   write(stderr, '(A,I0)') "Itemdir entries:    ", array_trees%total_itemdir_entries
   write(stderr, '(A,I0,A,I0,A)') "Itemdir size:       ", array_trees%itemdir_size1, &
                                  " / ", array_trees%itemdir_size2, ""
   write(stderr, '(A)') "Signatures:         Fixed-size arrays (simplified!)"
   write(stderr, *)
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
   integer :: item_idx

   ! Print items1
   item_idx = array_trees%parts(part_idx)%first_item1_idx
   do while (item_idx > 0 .and. item_idx <= array_trees%total_items1)
      write(stderr, '(1X,I0)', advance='no') array_trees%items1(item_idx)%value
      item_idx = array_trees%items1(item_idx)%next_item_idx
   end do

   write(stderr, '(A)', advance='no') ' /'

   ! Print items2
   item_idx = array_trees%parts(part_idx)%first_item2_idx
   do while (item_idx > 0 .and. item_idx <= array_trees%total_items2)
      write(stderr, '(1X,I0)', advance='no') array_trees%items2(item_idx)%value
      item_idx = array_trees%items2(item_idx)%next_item_idx
   end do

   write(stderr, *)
end subroutine

! Array-based tree printing procedures
! Add these to the array_trees module

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

      ! Print part index with item counts (analogous to address + counts)
      write(stderr, '(A,I0,1X,A,I0,A,I0,A)') prefix(1:pos-1), child_idx, &
         '(', array_trees%parts(child_idx)%num_items1, '/', &
         array_trees%parts(child_idx)%num_items2, ')'

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

      ! Print the split part index with item counts (no chain index)
      split_part_idx = array_trees%chains(child_idx)%split_part_idx
      if (split_part_idx > 0) then
         write(stderr, '(A,I0,1X,A,I0,A,I0,A)') prefix(1:pos-1), split_part_idx, &
            '(', array_trees%parts(split_part_idx)%num_items1, '/', &
            array_trees%parts(split_part_idx)%num_items2, ')'
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
   integer :: i, link_idx, partref_idx

   write(stderr, *)
   write(stderr, '(A)') repeat("=", 40)
   write(stderr, '(A)') "        Chain Details (Array)"
   write(stderr, '(A)') repeat("=", 40)
   write(stderr, *)

   do i = 1, array_trees%total_chains
      write(stderr, '(A,I0,A,I0,A,I0,A,I0,A)') 'Chain ', i, ': ', &
         array_trees%chains(i)%num_links, ' links, ', &
         array_trees%chains(i)%num_children, ' children'

      if (array_trees%chains(i)%split_part_idx > 0) then
         write(stderr, '(A,I0)') '  Split part: ', array_trees%chains(i)%split_part_idx
      end if

      ! Show links in this chain
      link_idx = array_trees%chains(i)%first_link_idx
      do while (link_idx > 0)
         write(stderr, '(A,I0,A,I0,A)', advance='no') '  Link ', link_idx, &
            ' (', array_trees%links(link_idx)%num_parts, ' parts): '

         ! Show parts in this link
         partref_idx = array_trees%links(link_idx)%first_partref_idx
         do while (partref_idx > 0)
            write(stderr, '(I0)', advance='no') array_trees%partrefs(partref_idx)%part_idx
            partref_idx = array_trees%partrefs(partref_idx)%next_ref_idx
            if (partref_idx > 0) write(stderr, '(A)', advance='no') ', '
         end do
         write(stderr, *)

         link_idx = array_trees%links(link_idx)%next_link_idx
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
