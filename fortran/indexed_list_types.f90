! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

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

module indexed_list_types
! Array representation of the partition and assignment trees. They are
! built as linked lists (linked_list_types), which suits their dynamic
! growth, and converted once to fixed-size arrays for the searches, which
! improves cache locality and avoids pointer chasing.
use parameters
use common_types
use linked_list_types
use chemdata
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

! Part of the HNA partition tree
type, public :: partree_item_t
   integer(ik) :: depth
   integer(ik) :: n_children
   ! Relationships (0 = null)
   integer(ik) :: parent_part_idx
   integer(ik) :: first_child_idx
   integer(ik) :: last_child_idx
   integer(ik) :: next_sibling_idx
   integer(ik), allocatable :: child_indices(:)
   ! Atoms of the part: atomidcs1(items1_offset+1 : items1_offset+n_items1)
   ! and likewise for molecule 2
   integer(ik) :: items1_offset, n_items1
   integer(ik) :: items2_offset, n_items2
   ! Signature as distinct values with their frequencies
   integer(ik) :: signature_values(MAX_COORDNUM)
   integer(ik) :: signature_frequencies(MAX_COORDNUM)
   integer(ik) :: signature_unique_count  ! number of distinct values
   integer(ik) :: signature_size          ! sum of frequencies
end type

! Link: one refinement level of an assignment tree node
type, public :: chain_item_t
   integer(ik) :: n_parts
   integer(ik) :: parent_chain_idx    ! node that owns this link
   ! Parts refined at this level:
   ! partref_entries(partref_offset+1 : partref_offset+n_parts)
   integer(ik) :: partref_offset
end type

! Assignment tree node
type, public :: assigntree_item_t
   integer(ik) :: n_atoms1, n_atoms2
   integer(ik) :: n_links, n_children
   ! Part individualized at this node (0 = none, for the root)
   integer(ik) :: split_part_idx
   ! Relationships (0 = null)
   integer(ik) :: parent_chain_idx
   integer(ik) :: first_child_idx
   integer(ik) :: last_child_idx
   integer(ik) :: next_sibling_idx
   integer(ik), allocatable :: child_indices(:)
   ! Links of the node: chain(link_offset+1 : link_offset+n_links)
   integer(ik) :: link_offset
end type

! Partition tree, assignment tree and adjacency lists in array form
type, public :: array_trees_t
   ! Atoms of all parts, in contiguous segments per part
   integer(ik), allocatable :: atomidcs1(:)
   integer(ik), allocatable :: atomidcs2(:)
   type(chain_item_t), allocatable :: chain(:)
   type(partree_item_t), allocatable :: partree(:)
   type(assigntree_item_t), allocatable :: assignment_tree(:)
   ! Vertex directories [atom, link]: part of each atom at each link
   ! (0 = not placed at that link). Each link is a column, so the
   ! consecutive links of a node are cleared as one contiguous block.
   integer(ik), allocatable :: itemdir1_entries(:,:)
   integer(ik), allocatable :: itemdir2_entries(:,:)
   ! Part indices listed by the links
   integer(ik), allocatable :: partref_entries(:)
   ! Adjacency lists [atom, neighbor] with coordination numbers and bond types
   integer(ik), allocatable :: adjcs1_cn(:)
   integer(ik), allocatable :: adjcs2_cn(:)
   integer(ik), allocatable :: adjcs1_list(:,:)
   integer(ik), allocatable :: adjcs2_list(:,:)
   integer(ik), allocatable :: adjcs1_bondtype(:,:)
   integer(ik), allocatable :: adjcs2_bondtype(:,:)
   ! Sizes
   integer(ik) :: n_atoms1, n_atoms2
   integer(ik) :: total_items1, total_items2, total_parts
   integer(ik) :: total_links, total_chains
   integer(ik) :: total_partref_entries
   ! Assignment tree statistics (Table 1 of the paper): number of complete
   ! assignments (product of branch possibilities) and sum of partial
   ! combinations evaluated by the local search
   real(real64) :: partial_combinations
   real(real64) :: total_combinations
end type

contains

subroutine cache_adjacency_lists(adjcs1, adjcs2, cache_arrays)
! Copy the adjacency lists of both molecules into cache_arrays
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik) :: i, n_atoms1, n_atoms2

   n_atoms1 = size(adjcs1)
   n_atoms2 = size(adjcs2)

   allocate(cache_arrays%adjcs1_cn(n_atoms1))
   allocate(cache_arrays%adjcs2_cn(n_atoms2))
   allocate(cache_arrays%adjcs1_list(n_atoms1, MAX_COORDNUM))
   allocate(cache_arrays%adjcs2_list(n_atoms2, MAX_COORDNUM))
   allocate(cache_arrays%adjcs1_bondtype(n_atoms1, MAX_COORDNUM))
   allocate(cache_arrays%adjcs2_bondtype(n_atoms2, MAX_COORDNUM))

   cache_arrays%adjcs1_list = 0
   cache_arrays%adjcs2_list = 0
   cache_arrays%adjcs1_bondtype = NO_BOND
   cache_arrays%adjcs2_bondtype = NO_BOND

   do i = 1, n_atoms1
      cache_arrays%adjcs1_cn(i) = adjcs1(i)%cn
      cache_arrays%adjcs1_list(i, 1:adjcs1(i)%cn) = adjcs1(i)%list(1:adjcs1(i)%cn)
      cache_arrays%adjcs1_bondtype(i, 1:adjcs1(i)%cn) = adjcs1(i)%bondtype(1:adjcs1(i)%cn)
   end do

   do i = 1, n_atoms2
      cache_arrays%adjcs2_cn(i) = adjcs2(i)%cn
      cache_arrays%adjcs2_list(i, 1:adjcs2(i)%cn) = adjcs2(i)%list(1:adjcs2(i)%cn)
      cache_arrays%adjcs2_bondtype(i, 1:adjcs2(i)%cn) = adjcs2(i)%bondtype(1:adjcs2(i)%cn)
   end do
end subroutine

subroutine cache_partition_tree(partition_tree, cache_arrays)
! Convert the linked partition tree into cache_arrays
   type(partition_node_t), pointer, intent(in) :: partition_tree
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik) :: item1_idx, item2_idx

   cache_arrays%total_parts = partition_tree%total_parts
   cache_arrays%total_items1 = partition_tree%total_items1
   cache_arrays%total_items2 = partition_tree%total_items2

   allocate(cache_arrays%atomidcs1(cache_arrays%total_items1))
   allocate(cache_arrays%atomidcs2(cache_arrays%total_items2))
   allocate(cache_arrays%partree(cache_arrays%total_parts))

   cache_arrays%atomidcs1 = 0
   cache_arrays%atomidcs2 = 0

   item1_idx = 0
   item2_idx = 0
   call convert_parts_recurse(partition_tree, cache_arrays, item1_idx, item2_idx)
end subroutine

subroutine cache_assignment_tree(assignment_tree, cache_arrays)
! Convert the linked assignment tree into cache_arrays and compute its
! combination statistics. The vertex directories start empty.
   type(chaintree_node_t), pointer, intent(in) :: assignment_tree
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik) :: partref_idx, link_idx

   cache_arrays%n_atoms1 = assignment_tree%n_atoms1
   cache_arrays%n_atoms2 = assignment_tree%n_atoms2
   cache_arrays%total_chains = assignment_tree%total_chains
   cache_arrays%total_links = assignment_tree%total_links
   cache_arrays%total_partref_entries = assignment_tree%total_partrefs

   allocate(cache_arrays%chain(cache_arrays%total_links))
   allocate(cache_arrays%assignment_tree(cache_arrays%total_chains))
   allocate(cache_arrays%partref_entries(cache_arrays%total_partref_entries))
   allocate(cache_arrays%itemdir1_entries(cache_arrays%n_atoms1, cache_arrays%total_links))
   allocate(cache_arrays%itemdir2_entries(cache_arrays%n_atoms2, cache_arrays%total_links))

   cache_arrays%partref_entries = 0
   cache_arrays%itemdir1_entries = 0
   cache_arrays%itemdir2_entries = 0

   partref_idx = 0
   link_idx = 0
   call convert_chains_recurse(assignment_tree, cache_arrays, partref_idx, link_idx, &
         cache_arrays%partial_combinations, cache_arrays%total_combinations)
end subroutine

subroutine convert_signature(part, cache_arrays, part_idx)
! Store the signature of part as distinct values with their frequencies
   type(partition_node_t), pointer, intent(in) :: part
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(in) :: part_idx
   integer(ik) :: temp_values(MAX_COORDNUM)
   integer(ik) :: temp_count, i, j, value
   logical(lk) :: found

   ! Edge codes (see adjacency::edge_code) of the neighbors with a part
   temp_count = size(part%signature)
   temp_values(:temp_count) = part%signature

   cache_arrays%partree(part_idx)%signature_size = temp_count

   cache_arrays%partree(part_idx)%signature_unique_count = 0
   do i = 1, temp_count
      value = temp_values(i)
      found = .FALSE.

      do j = 1, cache_arrays%partree(part_idx)%signature_unique_count
         if (cache_arrays%partree(part_idx)%signature_values(j) == value) then
            cache_arrays%partree(part_idx)%signature_frequencies(j) = &
               cache_arrays%partree(part_idx)%signature_frequencies(j) + 1
            found = .TRUE.
            exit
         end if
      end do

      if (.not. found) then
         cache_arrays%partree(part_idx)%signature_unique_count = &
            cache_arrays%partree(part_idx)%signature_unique_count + 1
         cache_arrays%partree(part_idx)%signature_values(cache_arrays%partree(part_idx)%signature_unique_count) = value
         cache_arrays%partree(part_idx)%signature_frequencies(cache_arrays%partree(part_idx)%signature_unique_count) = 1
      end if
   end do

   ! Zero the unused entries
   cache_arrays%partree(part_idx)%signature_values(cache_arrays%partree(part_idx)%signature_unique_count + 1:MAX_COORDNUM) = 0
   cache_arrays%partree(part_idx)%signature_frequencies(cache_arrays%partree(part_idx)%signature_unique_count + 1:MAX_COORDNUM) = 0
end subroutine

recursive subroutine convert_parts_recurse(part, cache_arrays, item1_idx, item2_idx)
! Store part and its subtree at their global indices. item1_idx and
! item2_idx are the last used positions of atomidcs1 and atomidcs2.
   type(partition_node_t), pointer, intent(in) :: part
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(inout) :: item1_idx, item2_idx
   type(partition_node_t), pointer :: child_part
   type(item_node_t), pointer :: item
   integer(ik) :: part_idx, i, n_children

   if (.not. associated(part)) return

   part_idx = part%global_idx

   cache_arrays%partree(part_idx)%depth = part%depth
   cache_arrays%partree(part_idx)%n_children = part%n_children

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

   ! Atom segments start after the last used positions
   cache_arrays%partree(part_idx)%items1_offset = item1_idx
   cache_arrays%partree(part_idx)%n_items1 = part%n_items1

   cache_arrays%partree(part_idx)%items2_offset = item2_idx
   cache_arrays%partree(part_idx)%n_items2 = part%n_items2

   call convert_signature(part, cache_arrays, part_idx)

   item => part%first_item1
   i = 0
   do while (associated(item))
      i = i + 1
      item1_idx = item1_idx + 1
      cache_arrays%atomidcs1(item1_idx) = item%idx
      item => item%next_item
   end do

   item => part%first_item2
   i = 0
   do while (associated(item))
      i = i + 1
      item2_idx = item2_idx + 1
      cache_arrays%atomidcs2(item2_idx) = item%idx
      item => item%next_item
   end do

   if (part%n_children > 0) then
      allocate(cache_arrays%partree(part_idx)%child_indices(part%n_children))
      child_part => part%first_child_part
      n_children = 0
      do while (associated(child_part))
         n_children = n_children + 1
         cache_arrays%partree(part_idx)%child_indices(n_children) = child_part%global_idx
         child_part => child_part%next_sibling_part
      end do
   end if

   child_part => part%first_child_part
   do while (associated(child_part))
      call convert_parts_recurse(child_part, cache_arrays, item1_idx, item2_idx)
      child_part => child_part%next_sibling_part
   end do
end subroutine

recursive subroutine convert_chains_recurse(chain, cache_arrays, partref_idx, link_idx, &
                                            partial_combinations, total_combinations)
! Store the node chain and its subtree at their global indices, and return
! the combination statistics of the subtree: total_combinations is the
! number of complete assignments (product over the children of the pair
! choices times their own totals) and partial_combinations the number of
! leaves visited by the local search (sum over the children of the same
! terms).
   type(chaintree_node_t), pointer, intent(in) :: chain
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(inout) :: partref_idx, link_idx
   real(rk), intent(out) :: partial_combinations, total_combinations
   type(chaintree_node_t), pointer :: child_chain
   type(chain_node_t), pointer :: link
   type(partref_node_t), pointer :: partref
   integer(ik) :: chain_idx, current_link_idx, n_children
   real(rk) :: child_combinations, child_product
   integer(ik) :: n_items2

   if (.not. associated(chain)) return

   if (chain%n_children == 0) then
      partial_combinations = 1
      total_combinations = 1
   else
      partial_combinations = 0
      total_combinations = 1
   end if

   chain_idx = chain%global_idx

   cache_arrays%assignment_tree(chain_idx)%n_atoms1 = chain%n_atoms1
   cache_arrays%assignment_tree(chain_idx)%n_atoms2 = chain%n_atoms2
   cache_arrays%assignment_tree(chain_idx)%n_links = chain%n_links
   cache_arrays%assignment_tree(chain_idx)%n_children = chain%n_children

   if (associated(chain%split_part)) then
      cache_arrays%assignment_tree(chain_idx)%split_part_idx = chain%split_part%global_idx
   else
      cache_arrays%assignment_tree(chain_idx)%split_part_idx = 0
   end if

   if (associated(chain%parent_chain)) then
      cache_arrays%assignment_tree(chain_idx)%parent_chain_idx = chain%parent_chain%global_idx
   else
      cache_arrays%assignment_tree(chain_idx)%parent_chain_idx = 0
   end if

   if (associated(chain%first_child_chain)) then
      cache_arrays%assignment_tree(chain_idx)%first_child_idx = chain%first_child_chain%global_idx
   else
      cache_arrays%assignment_tree(chain_idx)%first_child_idx = 0
   end if

   if (associated(chain%last_child_chain)) then
      cache_arrays%assignment_tree(chain_idx)%last_child_idx = chain%last_child_chain%global_idx
   else
      cache_arrays%assignment_tree(chain_idx)%last_child_idx = 0
   end if

   if (associated(chain%next_sibling_chain)) then
      cache_arrays%assignment_tree(chain_idx)%next_sibling_idx = chain%next_sibling_chain%global_idx
   else
      cache_arrays%assignment_tree(chain_idx)%next_sibling_idx = 0
   end if

   if (chain%n_children > 0) then
      allocate(cache_arrays%assignment_tree(chain_idx)%child_indices(chain%n_children))
      child_chain => chain%first_child_chain
      n_children = 0
      do while (associated(child_chain))
         n_children = n_children + 1
         cache_arrays%assignment_tree(chain_idx)%child_indices(n_children) = child_chain%global_idx
         child_chain => child_chain%next_sibling_chain
      end do
   end if

   ! Link segment starts after the last used position
   cache_arrays%assignment_tree(chain_idx)%link_offset = link_idx

   link => chain%first_link
   do while (associated(link))
      link_idx = link_idx + 1
      current_link_idx = link_idx

      cache_arrays%chain(current_link_idx)%n_parts = link%n_parts
      cache_arrays%chain(current_link_idx)%parent_chain_idx = chain%global_idx

      cache_arrays%chain(current_link_idx)%partref_offset = partref_idx

      partref => link%first_partref
      do while (associated(partref))
         partref_idx = partref_idx + 1
         cache_arrays%partref_entries(partref_idx) = partref%part%global_idx
         partref => partref%nextref
      end do

      link => link%next_link
   end do

   child_chain => chain%first_child_chain
   do while (associated(child_chain))
      call convert_chains_recurse(child_chain, cache_arrays, partref_idx, link_idx, &
                                  child_combinations, child_product)

      if (chain%n_children > 0) then
         ! The first atom of molecule 1 is paired with each atom of molecule 2
         n_items2 = child_chain%split_part%n_items2
         partial_combinations = partial_combinations + (n_items2 * child_combinations)
         total_combinations = total_combinations * (n_items2 * child_product)
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

   ! The root part is at index 1
   call print_items_recursive_array(cache_arrays, 1)
   write(stderr, *)
end subroutine

recursive subroutine print_items_recursive_array(cache_arrays, part_idx)
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: part_idx
   integer(ik) :: child_idx, i

   do i = 1, cache_arrays%partree(part_idx)%n_children
      child_idx = cache_arrays%partree(part_idx)%child_indices(i)

      write(stderr, '(A,I0,A)', advance='no') "Part ", child_idx, ':'
      call print_part_items_array(cache_arrays, child_idx)

      call print_items_recursive_array(cache_arrays, child_idx)
   end do
end subroutine

subroutine print_part_items_array(cache_arrays, part_idx)
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: part_idx
   integer(ik) :: i

   do i = 1, cache_arrays%partree(part_idx)%n_items1
      write(stderr, '(1X,I0)', advance='no') cache_arrays%atomidcs1(cache_arrays%partree(part_idx)%items1_offset + i)
   end do

   write(stderr, '(A)', advance='no') ' /'

   do i = 1, cache_arrays%partree(part_idx)%n_items2
      write(stderr, '(1X,I0)', advance='no') cache_arrays%atomidcs2(cache_arrays%partree(part_idx)%items2_offset + i)
   end do

   write(stderr, *)
end subroutine

subroutine print_part_tree_array(cache_arrays)
   type(array_trees_t), intent(in) :: cache_arrays
   logical(lk), dimension(:), allocatable :: is_last_child

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
   is_last_child = .FALSE.

   write(stderr, '(A)') 'ROOT'
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

   call print_signatures_recursive_array(cache_arrays, 1)
   write(stderr, *)
end subroutine

recursive subroutine print_signatures_recursive_array(cache_arrays, part_idx)
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: part_idx
   integer(ik) :: child_idx, i, j

   if (part_idx == 0) return

   do i = 1, cache_arrays%partree(part_idx)%n_children
      child_idx = cache_arrays%partree(part_idx)%child_indices(i)

      write(stderr,'(A,I0,A,I0,A)',advance='no') 'Part ', child_idx, ' (len=', &
         cache_arrays%partree(child_idx)%signature_size, '):'

      ! Distinct values as part_index:bond_type x frequency
      if (cache_arrays%partree(child_idx)%signature_unique_count > 0) then
         write(stderr, '(A)', advance='no') ' ['
         do j = 1, cache_arrays%partree(child_idx)%signature_unique_count
            if (j > 1) write(stderr, '(A)', advance='no') ', '
            write(stderr, '(I0,A,I0,A,I0)', advance='no') &
               cache_arrays%partree(child_idx)%signature_values(j)/BOND_TYPE_RADIX, ':', &
               modulo(cache_arrays%partree(child_idx)%signature_values(j), BOND_TYPE_RADIX), '×', &
               cache_arrays%partree(child_idx)%signature_frequencies(j)
         end do
         write(stderr, '(A)') ']'
      else
         write(stderr, '(A)') ' []'
      end if

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

   call print_leaf_items_recursive_array(cache_arrays, 1)
   write(stderr, *)
end subroutine

recursive subroutine print_leaf_items_recursive_array(cache_arrays, part_idx)
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: part_idx
   integer(ik) :: child_idx, i

   if (part_idx == 0) return

   do i = 1, cache_arrays%partree(part_idx)%n_children
      child_idx = cache_arrays%partree(part_idx)%child_indices(i)

      if (cache_arrays%partree(child_idx)%n_children == 0) then
         write(stderr, '(A,I0,A)', advance='no') 'Part ', child_idx, ':'
         call print_part_items_array(cache_arrays, child_idx)
      end if

      call print_leaf_items_recursive_array(cache_arrays, child_idx)
   end do
end subroutine

subroutine print_chain_details_array(cache_arrays)
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik) :: i, link_idx, j, k, part_idx

   write(stderr, '(A)') repeat("=", 40)
   write(stderr, '(A)') "        Chain Details"
   write(stderr, '(A)') repeat("=", 40)
   write(stderr, *)

   do i = 1, cache_arrays%total_chains
      write(stderr, '(A,I0,A,I0,A,I0,A)') 'Chain ', i, ': ', &
         cache_arrays%assignment_tree(i)%n_links, ' links, ', &
         cache_arrays%assignment_tree(i)%n_children, ' children'

      if (cache_arrays%assignment_tree(i)%split_part_idx > 0) then
         write(stderr, '(A,I0)') '  Split part: ', cache_arrays%assignment_tree(i)%split_part_idx
      end if

      do j = 1, cache_arrays%assignment_tree(i)%n_links
         link_idx = cache_arrays%assignment_tree(i)%link_offset + j
         write(stderr, '(A,I0,A,I0,A)', advance='no') '  Link ', link_idx, &
            ' (', cache_arrays%chain(link_idx)%n_parts, ' parts): '

         do k = 1, cache_arrays%chain(link_idx)%n_parts
            part_idx = cache_arrays%partref_entries(cache_arrays%chain(link_idx)%partref_offset + k)
            write(stderr, '(I0)', advance='no') part_idx
            if (k < cache_arrays%chain(link_idx)%n_parts) write(stderr, '(A)', advance='no') ', '
         end do
         write(stderr, *)
      end do

      if (i < cache_arrays%total_chains) write(stderr, *)
   end do

   write(stderr, *)
end subroutine

subroutine print_first_level_items_array(cache_arrays)
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik) :: child_idx, i

   if (cache_arrays%total_parts == 0) then
      write(stderr, '(A)') "Part tree is empty"
      return
   end if

   write(stderr, '(A)') repeat("=", 25)
   write(stderr, '(A)') "  First Level Items"
   write(stderr, '(A)') repeat("=", 25)
   write(stderr, *)

   do i = 1, cache_arrays%partree(1)%n_children
      child_idx = cache_arrays%partree(1)%child_indices(i)
      write(stderr, '(A,I0,A)', advance='no') 'Part ', child_idx, ':'
      call print_part_items_array(cache_arrays, child_idx)
   end do

   write(stderr, *)
end subroutine

recursive subroutine print_part_recursive_array(cache_arrays, part_idx, depth, is_last_child)
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: part_idx, depth
   logical(lk), dimension(:), intent(inout) :: is_last_child
   integer(ik) :: child_idx, i, j

   if (part_idx == 0) return

   do i = 1, cache_arrays%partree(part_idx)%n_children
      child_idx = cache_arrays%partree(part_idx)%child_indices(i)

      ! Tree drawing prefix
      is_last_child(depth + 1) = (i == cache_arrays%partree(part_idx)%n_children)
      do j = 1, depth
         if (is_last_child(j)) then
            write(stderr, '(A)', advance='no') "   "
         else
            write(stderr, '(A)', advance='no') "|  "
         end if
      end do

      if (is_last_child(depth + 1)) then
         write(stderr, '(A)', advance='no') "`--"
      else
         write(stderr, '(A)', advance='no') "|--"
      end if

      write(stderr, '(A,I0,A,I0,A)') '* (', &
         cache_arrays%partree(child_idx)%n_items1, '/', &
         cache_arrays%partree(child_idx)%n_items2, ')'

      call print_part_recursive_array(cache_arrays, child_idx, depth + 1, is_last_child)
   end do
end subroutine

subroutine print_chain_tree_array(atomtypes, cache_arrays)
! Print the assignment tree (as in Figure 2 of the paper, one node per
! split part labeled element*size) and its combination statistics
   type(partition_t), intent(in) :: atomtypes
   type(array_trees_t), intent(in) :: cache_arrays
   logical(lk), dimension(:), allocatable :: is_last_child

   if (cache_arrays%total_chains == 0) then
      write(stdout, '(A)') "Assignment tree is empty"
      return
   end if

   ! Allocate tracking array for tree lines (max depth 100)
   allocate(is_last_child(100))
   is_last_child = .FALSE.

   write(stdout, '(A)') '*'
   call print_chain_recursive_array(atomtypes, cache_arrays, 1, 0, is_last_child)
   write(stdout, *)

   if (cache_arrays%total_combinations < 1.0E+12) then
      write(stdout, '(A,I0)') "Total combinations: ", int(cache_arrays%total_combinations, kind=int64)
   else
      write(stdout, '(A,ES10.4)') "Total combinations: ", cache_arrays%total_combinations
   end if
   if (cache_arrays%partial_combinations < 1.0E+12) then
      write(stdout, '(A,I0)') "Sum of partial combinations: ", int(cache_arrays%partial_combinations, kind=int64)
   else
      write(stdout, '(A,ES10.4)') "Sum of partial combinations: ", cache_arrays%partial_combinations
   end if
   write(stdout, *)

   deallocate(is_last_child)
end subroutine

recursive subroutine print_chain_recursive_array(atomtypes, cache_arrays, chain_idx, depth, is_last_child)
   type(partition_t), intent(in) :: atomtypes
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: chain_idx, depth
   logical(lk), dimension(:), intent(inout) :: is_last_child
   integer(ik) :: child_idx, split_part_idx, first_atom_idx
   integer(ik) :: i, j

   if (chain_idx == 0) return

   do i = 1, cache_arrays%assignment_tree(chain_idx)%n_children
      child_idx = cache_arrays%assignment_tree(chain_idx)%child_indices(i)

      ! Tree drawing prefix
      is_last_child(depth + 1) = (i == cache_arrays%assignment_tree(chain_idx)%n_children)
      do j = 1, depth
         if (is_last_child(j)) then
            write(stdout, '(A)', advance='no') "   "
         else
            write(stdout, '(A)', advance='no') "|  "
         end if
      end do

      if (is_last_child(depth + 1)) then
         write(stdout, '(A)', advance='no') "'--"
      else
         write(stdout, '(A)', advance='no') "|--"
      end if

      ! Element and size of the split part
      split_part_idx = cache_arrays%assignment_tree(child_idx)%split_part_idx
      first_atom_idx = cache_arrays%atomidcs1(cache_arrays%partree(split_part_idx)%items1_offset+1)
      write(stdout, '(A,"*",I0)') &
         trim(atomic_symbols(atomtypes%parts(atomtypes%itemdir1(first_atom_idx))%elnum)), &
         cache_arrays%partree(split_part_idx)%n_items1

      call print_chain_recursive_array(atomtypes, cache_arrays, child_idx, depth + 1, is_last_child)
   end do
end subroutine

end module
