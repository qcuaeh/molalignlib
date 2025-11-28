! MolAlignLib
! Copyright (C) 2022 José M. Vásquez

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module types_indexed
use parameters
use types_basic
use types_linked
use chemistry
use adjacency
implicit none
private
public cache_adjacency_lists
public cache_partition_tree
public cache_assignment_tree
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
   integer :: signature_size            ! total signature length (sum of frequencies)
end type

type, public :: chain_item_t
   integer :: num_parts
   integer :: parent_chain_idx    ! which chain owns this link
   ! Part reference segment in flattened array (using offset)
   integer :: partref_offset      ! offset into partref_entries array
   ! NOTE: itemdir1_offset and itemdir2_offset REMOVED - no longer needed with 2D arrays
end type

type, public :: assigntree_item_t
   integer :: atoms1_size, atoms2_size
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
   integer, allocatable :: atomidcs1(:)
   integer, allocatable :: atomidcs2(:)
   type(chain_item_t), allocatable :: chain(:)
   type(partree_item_t), allocatable :: partree(:)
   type(assigntree_item_t), allocatable :: assigntree(:)
   ! Flattened variable-length data - all pure integer arrays!
   integer, allocatable :: itemdir1_entries(:,:)  ! [link_idx, atom_idx]
   integer, allocatable :: itemdir2_entries(:,:)  ! [link_idx, atom_idx]
   integer, allocatable :: partref_entries(:)     ! Part indices for partrefs
   ! Adjacency information stored directly for fastest access
   integer, allocatable :: adjcs1_cn(:)         ! Count for each adjcs1 atom's adjacency list
   integer, allocatable :: adjcs2_cn(:)         ! Count for each adjcs2 atom's adjacency list
   integer, allocatable :: adjcs1_list(:,:)     ! Direct 2D adjacency lists for adjcs1 [atom_idx, neighbor_idx]
   integer, allocatable :: adjcs2_list(:,:)     ! Direct 2D adjacency lists for adjcs2 [atom_idx, neighbor_idx]
   ! Metadata
   integer :: atoms1_size, atoms2_size  ! number of atoms in each molecule
   integer :: total_items1, total_items2, total_parts
   integer :: total_links, total_chains
   integer :: total_partref_entries
   ! Assignment statistics
   real(rk) :: partial_combinations
   real(rk) :: total_combinations
end type

contains

subroutine cache_adjacency_lists(adjcs1, adjcs2, cache_arrays)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(array_trees_t), intent(inout) :: cache_arrays
   integer :: i, atoms1_size, atoms2_size

   atoms1_size = size(adjcs1)
   atoms2_size = size(adjcs2)

   ! Allocate adjacency arrays
   allocate(cache_arrays%adjcs1_cn(atoms1_size))
   allocate(cache_arrays%adjcs2_cn(atoms2_size))
   allocate(cache_arrays%adjcs1_list(atoms1_size, MAX_COORD))
   allocate(cache_arrays%adjcs2_list(atoms2_size, MAX_COORD))

   ! Initialize adjacency lists to zero
   cache_arrays%adjcs1_list = 0
   cache_arrays%adjcs2_list = 0

   ! Copy adjacency data for molecule 1
   do i = 1, atoms1_size
      cache_arrays%adjcs1_cn(i) = adjcs1(i)%cn
      cache_arrays%adjcs1_list(i, 1:adjcs1(i)%cn) = adjcs1(i)%list(1:adjcs1(i)%cn)
   end do

   ! Copy adjacency data for molecule 2
   do i = 1, atoms2_size
      cache_arrays%adjcs2_cn(i) = adjcs2(i)%cn
      cache_arrays%adjcs2_list(i, 1:adjcs2(i)%cn) = adjcs2(i)%list(1:adjcs2(i)%cn)
   end do
end subroutine

subroutine cache_partition_tree(partition_tree, cache_arrays)
   type(partition_node_t), pointer, intent(in) :: partition_tree
   type(array_trees_t), intent(inout) :: cache_arrays
   integer :: item1_idx, item2_idx

   ! Set part tree metadata
   cache_arrays%total_parts = partition_tree%total_parts
   cache_arrays%total_items1 = partition_tree%total_items1
   cache_arrays%total_items2 = partition_tree%total_items2

   ! Allocate part tree arrays
   allocate(cache_arrays%atomidcs1(cache_arrays%total_items1))
   allocate(cache_arrays%atomidcs2(cache_arrays%total_items2))
   allocate(cache_arrays%partree(cache_arrays%total_parts))

   ! Initialize arrays
   cache_arrays%atomidcs1 = 0
   cache_arrays%atomidcs2 = 0

   ! Convert part tree
   item1_idx = 0
   item2_idx = 0
   call convert_parts_recurse(partition_tree, cache_arrays, item1_idx, item2_idx)
end subroutine

subroutine cache_assignment_tree(assignment_tree, cache_arrays)
   type(chaintree_node_t), pointer, intent(in) :: assignment_tree
   type(array_trees_t), intent(inout) :: cache_arrays
   integer :: partref_idx, link_idx

   ! Set assignment tree metadata
   cache_arrays%atoms1_size = assignment_tree%atoms1_size
   cache_arrays%atoms2_size = assignment_tree%atoms2_size
   cache_arrays%total_chains = assignment_tree%total_chains
   cache_arrays%total_links = assignment_tree%total_links
   cache_arrays%total_partref_entries = assignment_tree%total_partrefs

   ! Allocate assignment tree arrays
   allocate(cache_arrays%chain(cache_arrays%total_links))
   allocate(cache_arrays%assigntree(cache_arrays%total_chains))
   allocate(cache_arrays%partref_entries(cache_arrays%total_partref_entries))
   allocate(cache_arrays%itemdir1_entries(cache_arrays%total_links, cache_arrays%atoms1_size))
   allocate(cache_arrays%itemdir2_entries(cache_arrays%total_links, cache_arrays%atoms2_size))

   ! Initialize arrays
   cache_arrays%partref_entries = 0
   cache_arrays%itemdir1_entries = 0
   cache_arrays%itemdir2_entries = 0

   ! Convert assignment tree
   partref_idx = 0
   link_idx = 0
   call convert_chains_recurse(assignment_tree, cache_arrays, partref_idx, link_idx, &
         cache_arrays%partial_combinations, cache_arrays%total_combinations)
end subroutine

subroutine convert_signature(part, cache_arrays, part_idx)
   type(partition_node_t), pointer, intent(in) :: part
   type(array_trees_t), intent(inout) :: cache_arrays
   integer, intent(in) :: part_idx
   integer :: temp_values(MAX_COORD)
   integer :: temp_count, i, j, value
   logical :: found

   ! First pass: collect all non-null signature values
   temp_count = 0
   do i = 1, size(part%signature)
      if (associated(part%signature(i)%ptr)) then
         temp_count = temp_count + 1
         temp_values(temp_count) = part%signature(i)%ptr%global_idx
      end if
   end do

   ! Store total signature length
   cache_arrays%partree(part_idx)%signature_size = temp_count

   ! Second pass: compute unique values and their frequencies
   cache_arrays%partree(part_idx)%signature_unique_count = 0
   do i = 1, temp_count
      value = temp_values(i)
      found = .false.

      ! Check if this value is already in unique list
      do j = 1, cache_arrays%partree(part_idx)%signature_unique_count
         if (cache_arrays%partree(part_idx)%signature_values(j) == value) then
            cache_arrays%partree(part_idx)%signature_frequencies(j) = &
               cache_arrays%partree(part_idx)%signature_frequencies(j) + 1
            found = .true.
            exit
         end if
      end do

      ! If not found, add as new unique value
      if (.not. found) then
         cache_arrays%partree(part_idx)%signature_unique_count = &
            cache_arrays%partree(part_idx)%signature_unique_count + 1
         cache_arrays%partree(part_idx)%signature_values(cache_arrays%partree(part_idx)%signature_unique_count) = value
         cache_arrays%partree(part_idx)%signature_frequencies(cache_arrays%partree(part_idx)%signature_unique_count) = 1
      end if
   end do

   ! Zero out unused entries using intrinsic operation
   cache_arrays%partree(part_idx)%signature_values(cache_arrays%partree(part_idx)%signature_unique_count + 1:MAX_COORD) = 0
   cache_arrays%partree(part_idx)%signature_frequencies(cache_arrays%partree(part_idx)%signature_unique_count + 1:MAX_COORD) = 0
end subroutine

recursive subroutine convert_parts_recurse(part, cache_arrays, item1_idx, item2_idx)
   type(partition_node_t), pointer, intent(in) :: part
   type(array_trees_t), intent(inout) :: cache_arrays
   integer, intent(inout) :: item1_idx, item2_idx
   type(partition_node_t), pointer :: child_part
   type(item_node_t), pointer :: item
   integer :: part_idx, i, child_count

   if (.not. associated(part)) return

   ! Convert this part (global indices always start at 1)
   part_idx = part%global_idx

   cache_arrays%partree(part_idx)%depth = part%depth
   cache_arrays%partree(part_idx)%num_children = part%num_children

   ! Relationships using global indices
   if (associated(part%parent_part)) then
      cache_arrays%partree(part_idx)%parent_part_idx = part%parent_part%global_idx
   else
      cache_arrays%partree(part_idx)%parent_part_idx = 0
   end if

   if (associated(part%first_child_part)) then
      cache_arrays%partree(part_idx)%first_child_idx = part%first_child_part%global_idx
   else
      cache_arrays%partree(part_idx)%first_child_idx = 0
   end if

   if (associated(part%last_child_part)) then
      cache_arrays%partree(part_idx)%last_child_idx = part%last_child_part%global_idx
   else
      cache_arrays%partree(part_idx)%last_child_idx = 0
   end if

   if (associated(part%next_sibling_part)) then
      cache_arrays%partree(part_idx)%next_sibling_idx = part%next_sibling_part%global_idx
   else
      cache_arrays%partree(part_idx)%next_sibling_idx = 0
   end if

   ! Set up item segments using offset approach (offset = start_idx - 1)
   cache_arrays%partree(part_idx)%items1_offset = item1_idx  ! item1_idx tracks the last used index
   cache_arrays%partree(part_idx)%items1_count = part%num_items1

   cache_arrays%partree(part_idx)%items2_offset = item2_idx  ! item2_idx tracks the last used index
   cache_arrays%partree(part_idx)%items2_count = part%num_items2

   ! OPTIMIZED: Convert signature to unique values, frequencies, and total length
   call convert_signature(part, cache_arrays, part_idx)

   ! Convert items1 to pure array format
   item => part%first_item1
   i = 0
   do while (associated(item))
      i = i + 1
      item1_idx = item1_idx + 1
      cache_arrays%atomidcs1(item1_idx) = item%idx
      item => item%next_item
   end do

   ! Convert items2 to pure array format
   item => part%first_item2
   i = 0
   do while (associated(item))
      i = i + 1
      item2_idx = item2_idx + 1
      cache_arrays%atomidcs2(item2_idx) = item%idx
      item => item%next_item
   end do

   ! OPTIMIZATION: Populate direct child access array for faster traversal
   if (part%num_children > 0) then
      allocate(cache_arrays%partree(part_idx)%child_indices(part%num_children))
      child_part => part%first_child_part
      child_count = 0
      do while (associated(child_part))
         child_count = child_count + 1
         cache_arrays%partree(part_idx)%child_indices(child_count) = child_part%global_idx
         child_part => child_part%next_sibling_part
      end do
   end if

   ! Recursively convert all children
   child_part => part%first_child_part
   do while (associated(child_part))
      call convert_parts_recurse(child_part, cache_arrays, item1_idx, item2_idx)
      child_part => child_part%next_sibling_part
   end do
end subroutine

recursive subroutine convert_chains_recurse(chain, cache_arrays, partref_idx, link_idx, &
                                            partial_combinations, total_combinations)
   type(chaintree_node_t), pointer, intent(in) :: chain
   type(array_trees_t), intent(inout) :: cache_arrays
   integer, intent(inout) :: partref_idx, link_idx
   real(rk), intent(out) :: partial_combinations, total_combinations
   type(chaintree_node_t), pointer :: child_chain
   type(chain_node_t), pointer :: link
   type(partref_node_t), pointer :: partref
   integer :: chain_idx, current_link_idx, child_count
   real(rk) :: child_combinations, child_product
   integer :: items2_count

   if (.not. associated(chain)) return

   ! If this is a leaf level (no children), both values are 1
   if (chain%num_children == 0) then
      partial_combinations = 1
      total_combinations = 1
   else
      ! Initialize accumulators for non-leaf nodes
      partial_combinations = 0
      total_combinations = 1
   end if

   ! Convert this chain (existing conversion logic)
   chain_idx = chain%global_idx

   cache_arrays%assigntree(chain_idx)%atoms1_size = chain%atoms1_size
   cache_arrays%assigntree(chain_idx)%atoms2_size = chain%atoms2_size
   cache_arrays%assigntree(chain_idx)%num_links = chain%num_links
   cache_arrays%assigntree(chain_idx)%num_children = chain%num_children

   ! Cross-tree reference
   if (associated(chain%split_part)) then
      cache_arrays%assigntree(chain_idx)%split_part_idx = chain%split_part%global_idx
   else
      cache_arrays%assigntree(chain_idx)%split_part_idx = 0
   end if

   ! Chain relationships
   if (associated(chain%parent_chain)) then
      cache_arrays%assigntree(chain_idx)%parent_chain_idx = chain%parent_chain%global_idx
   else
      cache_arrays%assigntree(chain_idx)%parent_chain_idx = 0
   end if

   if (associated(chain%first_child_chain)) then
      cache_arrays%assigntree(chain_idx)%first_child_idx = chain%first_child_chain%global_idx
   else
      cache_arrays%assigntree(chain_idx)%first_child_idx = 0
   end if

   if (associated(chain%last_child_chain)) then
      cache_arrays%assigntree(chain_idx)%last_child_idx = chain%last_child_chain%global_idx
   else
      cache_arrays%assigntree(chain_idx)%last_child_idx = 0
   end if

   if (associated(chain%next_sibling_chain)) then
      cache_arrays%assigntree(chain_idx)%next_sibling_idx = chain%next_sibling_chain%global_idx
   else
      cache_arrays%assigntree(chain_idx)%next_sibling_idx = 0
   end if

   ! OPTIMIZATION: Populate direct child access array for faster traversal
   if (chain%num_children > 0) then
      allocate(cache_arrays%assigntree(chain_idx)%child_indices(chain%num_children))
      child_chain => chain%first_child_chain
      child_count = 0
      do while (associated(child_chain))
         child_count = child_count + 1
         cache_arrays%assigntree(chain_idx)%child_indices(child_count) = child_chain%global_idx
         child_chain => child_chain%next_sibling_chain
      end do
   end if

   ! Set link offset (offset = start_idx - 1)
   cache_arrays%assigntree(chain_idx)%link_offset = link_idx

   ! Convert links in this chain using offset approach
   link => chain%first_link
   do while (associated(link))
      link_idx = link_idx + 1
      current_link_idx = link_idx

      cache_arrays%chain(current_link_idx)%num_parts = link%num_parts
      cache_arrays%chain(current_link_idx)%parent_chain_idx = chain%global_idx

      ! Set partref offset (offset = start_idx - 1)
      cache_arrays%chain(current_link_idx)%partref_offset = partref_idx

      ! Convert partrefs to pure array format
      partref => link%first_partref
      do while (associated(partref))
         partref_idx = partref_idx + 1
         cache_arrays%partref_entries(partref_idx) = partref%part%global_idx
         partref => partref%nextref
      end do

      link => link%next_link
   end do

   ! Recursively convert child chains and accumulate statistics
   child_chain => chain%first_child_chain
   do while (associated(child_chain))
      call convert_chains_recurse(child_chain, cache_arrays, partref_idx, link_idx, &
                                  child_combinations, child_product)

      ! Count statistics if this is not a leaf
      if (chain%num_children > 0) then
         items2_count = child_chain%split_part%num_items2

         ! Update combinations (sum): first item1 with each item2
         partial_combinations = partial_combinations + (items2_count * child_combinations)

         ! Update product (multiply): split sizes only
         total_combinations = total_combinations * (items2_count * child_product)
      end if

      child_chain => child_chain%next_sibling_chain
   end do
end subroutine

subroutine print_tree_items_array(cache_arrays)
   type(array_trees_t), intent(in) :: cache_arrays

   if (cache_arrays%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "      Part Items"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Root part is always at index 1, print its children recursively
   call print_items_recursive_array(cache_arrays, 1)
   write(stderr, *)
end subroutine

recursive subroutine print_items_recursive_array(cache_arrays, part_idx)
   type(array_trees_t), intent(in) :: cache_arrays
   integer, intent(in) :: part_idx
   integer :: child_idx, i

   ! Use direct array access instead of linked traversal for better performance
   do i = 1, cache_arrays%partree(part_idx)%num_children
      child_idx = cache_arrays%partree(part_idx)%child_indices(i)

      ! Print the child items with part index prefix
      write(stderr, '(A,I0,A)', advance='no') "Part ", child_idx, ':'
      call print_part_items_array(cache_arrays, child_idx)

      ! Recursively print this child's children
      call print_items_recursive_array(cache_arrays, child_idx)
   end do
end subroutine

subroutine print_part_items_array(cache_arrays, part_idx)
   type(array_trees_t), intent(in) :: cache_arrays
   integer, intent(in) :: part_idx
   integer :: i

   ! Print items1 using offset-based access
   do i = 1, cache_arrays%partree(part_idx)%items1_count
      write(stderr, '(1X,I0)', advance='no') cache_arrays%atomidcs1(cache_arrays%partree(part_idx)%items1_offset + i)
   end do

   write(stderr, '(A)', advance='no') ' /'

   ! Print items2 using offset-based access
   do i = 1, cache_arrays%partree(part_idx)%items2_count
      write(stderr, '(1X,I0)', advance='no') cache_arrays%atomidcs2(cache_arrays%partree(part_idx)%items2_offset + i)
   end do

   write(stderr, *)
end subroutine

! Array-based tree printing procedures

subroutine print_part_tree_array(cache_arrays)
   type(array_trees_t), intent(in) :: cache_arrays
   logical, dimension(:), allocatable :: is_last_child

   if (cache_arrays%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

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
   call print_part_recursive_array(cache_arrays, 1, 0, is_last_child)
   write(stderr, *)

   deallocate(is_last_child)
end subroutine

subroutine print_part_signatures_array(cache_arrays)
   type(array_trees_t), intent(in) :: cache_arrays

   if (cache_arrays%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "  Part Signatures"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Print signatures for all parts (root is always at index 1)
   call print_signatures_recursive_array(cache_arrays, 1)
   write(stderr, *)
end subroutine

recursive subroutine print_signatures_recursive_array(cache_arrays, part_idx)
   type(array_trees_t), intent(in) :: cache_arrays
   integer, intent(in) :: part_idx
   integer :: child_idx, i, j

   if (part_idx == 0) return

   ! Process all children using direct array access
   do i = 1, cache_arrays%partree(part_idx)%num_children
      child_idx = cache_arrays%partree(part_idx)%child_indices(i)

      ! Print the child signature with frequencies and total length
      write(stderr,'(A,I0,A,I0,A)',advance='no') 'Part ', child_idx, ' (len=', &
         cache_arrays%partree(child_idx)%signature_size, '):'

      ! Print unique signature values with frequencies
      if (cache_arrays%partree(child_idx)%signature_unique_count > 0) then
         write(stderr, '(A)', advance='no') ' ['
         do j = 1, cache_arrays%partree(child_idx)%signature_unique_count
            if (j > 1) write(stderr, '(A)', advance='no') ', '
            write(stderr, '(I0,A,I0)', advance='no') &
               cache_arrays%partree(child_idx)%signature_values(j), '×', &
               cache_arrays%partree(child_idx)%signature_frequencies(j)
         end do
         write(stderr, '(A)') ']'
      else
         write(stderr, '(A)') ' []'
      end if

      ! Recursively print this child's children
      call print_signatures_recursive_array(cache_arrays, child_idx)
   end do
end subroutine

subroutine print_leaf_items_array(cache_arrays)
   type(array_trees_t), intent(in) :: cache_arrays

   if (cache_arrays%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "   Leaf Items"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Print leaf items recursively (root is always at index 1)
   call print_leaf_items_recursive_array(cache_arrays, 1)
   write(stderr, *)
end subroutine

recursive subroutine print_leaf_items_recursive_array(cache_arrays, part_idx)
   type(array_trees_t), intent(in) :: cache_arrays
   integer, intent(in) :: part_idx
   integer :: child_idx, i

   if (part_idx == 0) return

   ! Process all children using direct array access
   do i = 1, cache_arrays%partree(part_idx)%num_children
      child_idx = cache_arrays%partree(part_idx)%child_indices(i)

      ! Only print items if this is a leaf part (no children)
      if (cache_arrays%partree(child_idx)%num_children == 0) then
         write(stderr, '(A,I0,A)', advance='no') 'Part ', child_idx, ':'
         call print_part_items_array(cache_arrays, child_idx)
      end if

      ! Recursively traverse this child's children to find more leaves
      call print_leaf_items_recursive_array(cache_arrays, child_idx)
   end do
end subroutine

subroutine print_chain_details_array(cache_arrays)
   type(array_trees_t), intent(in) :: cache_arrays
   integer :: i, link_idx, j, k, part_idx

   write(stderr, '(A)') repeat("=", 40)
   write(stderr, '(A)') "        Chain Details"
   write(stderr, '(A)') repeat("=", 40)
   write(stderr, *)

   do i = 1, cache_arrays%total_chains
      write(stderr, '(A,I0,A,I0,A,I0,A)') 'Chain ', i, ': ', &
         cache_arrays%assigntree(i)%num_links, ' links, ', &
         cache_arrays%assigntree(i)%num_children, ' children'

      if (cache_arrays%assigntree(i)%split_part_idx > 0) then
         write(stderr, '(A,I0)') '  Split part: ', cache_arrays%assigntree(i)%split_part_idx
      end if

      ! Show links in this chain using offset-based access
      do j = 1, cache_arrays%assigntree(i)%num_links
         link_idx = cache_arrays%assigntree(i)%link_offset + j
         write(stderr, '(A,I0,A,I0,A)', advance='no') '  Link ', link_idx, &
            ' (', cache_arrays%chain(link_idx)%num_parts, ' parts): '

         ! Show parts in this link using offset-based access
         do k = 1, cache_arrays%chain(link_idx)%num_parts
            part_idx = cache_arrays%partref_entries(cache_arrays%chain(link_idx)%partref_offset + k)
            write(stderr, '(I0)', advance='no') part_idx
            if (k < cache_arrays%chain(link_idx)%num_parts) write(stderr, '(A)', advance='no') ', '
         end do
         write(stderr, *)
      end do

      if (i < cache_arrays%total_chains) write(stderr, *)
   end do

   write(stderr, *)
end subroutine

subroutine print_first_level_items_array(cache_arrays)
   type(array_trees_t), intent(in) :: cache_arrays
   integer :: child_idx, i

   if (cache_arrays%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "  First Level Items"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   ! Print items for all parts at first level using direct array access
   do i = 1, cache_arrays%partree(1)%num_children
      child_idx = cache_arrays%partree(1)%child_indices(i)
      write(stderr, '(A,I0,A)', advance='no') 'Part ', child_idx, ':'
      call print_part_items_array(cache_arrays, child_idx)
   end do

   write(stderr, *)
end subroutine

recursive subroutine print_part_recursive_array(cache_arrays, part_idx, depth, is_last_child)
   type(array_trees_t), intent(in) :: cache_arrays
   integer, intent(in) :: part_idx, depth
   logical, dimension(:), intent(inout) :: is_last_child
   integer :: child_idx, i, j

   if (part_idx == 0) return

   ! Process all children using direct array access
   do i = 1, cache_arrays%partree(part_idx)%num_children
      child_idx = cache_arrays%partree(part_idx)%child_indices(i)

      ! Check if this is the last child
      is_last_child(depth + 1) = (i == cache_arrays%partree(part_idx)%num_children)

      ! Print prefix components directly
      do j = 1, depth
         if (is_last_child(j)) then
            write(stderr, '(A)', advance='no') "   "
         else
            write(stderr, '(A)', advance='no') "|  "
         end if
      end do

      ! Add branch characters
      if (is_last_child(depth + 1)) then
         write(stderr, '(A)', advance='no') "`--"
      else
         write(stderr, '(A)', advance='no') "|--"
      end if

      ! Print part index with item counts
      write(stderr, '(A,I0,A,I0,A)') '* (', &
         cache_arrays%partree(child_idx)%items1_count, '/', &
         cache_arrays%partree(child_idx)%items2_count, ')'

      ! Recursively print this child's children
      call print_part_recursive_array(cache_arrays, child_idx, depth + 1, is_last_child)
   end do
end subroutine

subroutine print_chain_tree_array(atomtypes, cache_arrays)
   type(partition_t), intent(in) :: atomtypes
   type(array_trees_t), intent(in) :: cache_arrays
   logical, dimension(:), allocatable :: is_last_child

   if (cache_arrays%total_chains == 0) then
      write(stdout, '(A)') "Assignment tree is empty"
      return
   end if

   ! Allocate tracking array for tree lines (max depth 100)
   allocate(is_last_child(100))
   is_last_child = .false.

   ! Print children recursively (root is always at index 1)
   write(stdout, '(A)') '*'
   call print_chain_recursive_array(atomtypes, cache_arrays, 1, 0, is_last_child)
   write(stdout, *)

   ! Print assignment statistics
   if (cache_arrays%total_combinations > 2**24) then
      write(stdout, '(A,ES8.2)') "Total combinations: ", cache_arrays%total_combinations
   else
      write(stdout, '(A,I0)') "Total combinations: ", int(cache_arrays%total_combinations)
   end if
   if (cache_arrays%partial_combinations > 2**24) then
      write(stdout, '(A,ES8.2)') "Sum of partial combinations: ", cache_arrays%partial_combinations
   else
      write(stdout, '(A,I0)') "Sum of partial combinations: ", int(cache_arrays%partial_combinations)
   end if
   write(stdout, *)

   deallocate(is_last_child)
end subroutine

recursive subroutine print_chain_recursive_array(atomtypes, cache_arrays, chain_idx, depth, is_last_child)
   type(partition_t), intent(in) :: atomtypes
   type(array_trees_t), intent(in) :: cache_arrays
   integer, intent(in) :: chain_idx, depth
   logical, dimension(:), intent(inout) :: is_last_child
   integer :: child_idx, split_part_idx, first_atom_idx
   integer :: i, j

   if (chain_idx == 0) return

   ! Process all children using direct array access instead of linked traversal
   do i = 1, cache_arrays%assigntree(chain_idx)%num_children
      child_idx = cache_arrays%assigntree(chain_idx)%child_indices(i)

      ! Check if this is the last child
      is_last_child(depth + 1) = (i == cache_arrays%assigntree(chain_idx)%num_children)

      ! Print prefix components directly
      do j = 1, depth
         if (is_last_child(j)) then
            write(stdout, '(A)', advance='no') "   "
         else
            write(stdout, '(A)', advance='no') "|  "
         end if
      end do

      ! Add branch characters
      if (is_last_child(depth + 1)) then
         write(stdout, '(A)', advance='no') "'--"
      else
         write(stdout, '(A)', advance='no') "|--"
      end if

      ! Print the split part index with item counts
      split_part_idx = cache_arrays%assigntree(child_idx)%split_part_idx
      first_atom_idx = cache_arrays%atomidcs1(cache_arrays%partree(split_part_idx)%items1_offset+1)
      write(stdout, '(A,"*",I0)') &
         trim(atomic_symbols(atomtypes%parts(atomtypes%itemdir1(first_atom_idx))%elnum)), &
         cache_arrays%partree(split_part_idx)%items1_count

      ! Recursively print this child's children
      call print_chain_recursive_array(atomtypes, cache_arrays, child_idx, depth + 1, is_last_child)
   end do
end subroutine

end module
