module mna_recompute_arrays
use parameters
use array_trees
implicit none

! Maximum possible number of children for a part
integer, parameter :: MAX_CHILDREN = 2

! Module-level signature workspace to eliminate allocations
integer :: signature_array(MAX_COORD)
integer :: signature_length

! DFS exploration variables
integer :: combination_count

contains

function signature_equivalence_array(array_trees, part_idx) result(equiv)
   ! OPTIMIZED: Fast path for length-1 signatures (most common case)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   logical :: equiv
   integer :: signature_frequencies, i, j

   if (signature_length /= array_trees%parts(part_idx)%signature_length) then
      equiv = .false.
      return
   end if

   ! FAST PATH: Direct comparison for length-1 signatures (most common)
   if (array_trees%parts(part_idx)%signature_length == 1) then
      equiv = (signature_array(1) == array_trees%parts(part_idx)%signature_values(1))
      return
   end if

   ! GENERIC PATH
   do i = 1, array_trees%parts(part_idx)%signature_unique_count
      signature_frequencies = 0

      ! Count matches in target signature
      do j = 1, signature_length
         if (signature_array(j) == array_trees%parts(part_idx)%signature_values(i)) then
            signature_frequencies = signature_frequencies + 1
         end if
      end do

      if (signature_frequencies /= array_trees%parts(part_idx)%signature_frequencies(i)) then
         equiv = .false.
         return
      end if
   end do

   equiv = .true.
end function

function calculate_leaf_squared_distance_contribution(coords1, coords2, array_trees, part_idx) result(squared_distance)
   ! Calculate squared distance contribution for a leaf part (items1_count == items2_count == 1)
   ! UPDATED: Now returns distance instead of accumulating globally
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)  ! coords(dimension, atom_index)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   real(rk) :: squared_distance
   integer :: item1_idx, item2_idx
   real(rk) :: dx, dy, dz

   ! Verify this is a leaf part with exactly one item from each molecule
   if (array_trees%parts(part_idx)%items1_count /= 1 .or. &
       array_trees%parts(part_idx)%items2_count /= 1) then
      write(stderr, '(A,I0,A,I0,A,I0)') 'Warning: Part ', part_idx, &
         ' is not a proper leaf (items1=', array_trees%parts(part_idx)%items1_count, &
         ', items2=', array_trees%parts(part_idx)%items2_count, ')'
      squared_distance = 0.0_rk
      return
   end if

   ! Get the assigned items
   item1_idx = array_trees%item1_values(array_trees%parts(part_idx)%items1_offset + 1)
   item2_idx = array_trees%item2_values(array_trees%parts(part_idx)%items2_offset + 1)

   ! Calculate squared distance between assigned atoms using coordinate matrices
   dx = coords1(1, item1_idx) - coords2(1, item2_idx)
   dy = coords1(2, item1_idx) - coords2(2, item2_idx)
   dz = coords1(3, item1_idx) - coords2(3, item2_idx)

   squared_distance = dx*dx + dy*dy + dz*dz
end function

function check_leaf_parts_for_squared_distance(coords1, coords2, array_trees, part_idx) result(total_distance)
   ! Check if child parts are leaves and calculate squared distance contributions
   ! UPDATED: Now returns total distance instead of accumulating globally
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   real(rk) :: total_distance
   integer :: i, child_idx

   total_distance = 0.0_rk

   ! Check all children of this part
   do i = 1, array_trees%parts(part_idx)%num_children
      child_idx = array_trees%parts(part_idx)%child_indices(i)

      ! If this child is a leaf (no children), calculate squared distance contribution
      if (array_trees%parts(child_idx)%num_children == 0) then
         total_distance = total_distance + calculate_leaf_squared_distance_contribution(coords1, coords2, array_trees, child_idx)
      end if
   end do
end function

function find_child_part_array(array_trees, parent_idx) result(child_relative_idx)
   ! OPTIMIZED: Assumes exactly 2 children - if signature doesn't match first, it must match second
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: parent_idx
   integer :: child_relative_idx, first_child_idx

   ! Get first child index directly
   first_child_idx = array_trees%parts(parent_idx)%child_indices(1)

   ! Check if signature matches first child
   if (signature_equivalence_array(array_trees, first_child_idx)) then
      child_relative_idx = 1
   else
      ! Must match second child (assumption: exactly 2 children)
      child_relative_idx = 2
   end if
end function

subroutine resplit_part_mna_array(coords1, coords2, array_trees, part_idx, read_link_idx, write_link_idx, total_distance)
! ULTRA-OPTIMIZED: Array-based version with direct 2D adjacency access for maximum performance
! UPDATED: Now accumulates distance contribution from newly created leaf parts into total_distance
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)  ! coords(dimension, atom_index)
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx, read_link_idx, write_link_idx
   real(rk), intent(inout) :: total_distance
   integer :: items1_trackers(MAX_CHILDREN), items2_trackers(MAX_CHILDREN)
   integer :: i, j, target_relative_idx, target_part_idx, item_value, target_idx, part_ref_idx, adj_atom
   integer :: items1_offset, items1_count, items2_offset, items2_count
   integer :: num_children

   ! Extract commonly used values for readability
   items1_offset = array_trees%parts(part_idx)%items1_offset
   items1_count = array_trees%parts(part_idx)%items1_count
   items2_offset = array_trees%parts(part_idx)%items2_offset
   items2_count = array_trees%parts(part_idx)%items2_count
   num_children = array_trees%parts(part_idx)%num_children

   ! INITIALIZATION: Reset redistribution counters using relative indices (1 to num_children)
   do i = 1, num_children
      items1_trackers(i) = 0
      items2_trackers(i) = 0
   end do

   ! Process first molecule items with ULTRA-OPTIMIZED signature generation using direct 2D adjacency
   do i = 1, items1_count
      item_value = array_trees%item1_values(items1_offset + i)

      ! ULTRA-OPTIMIZED: Generate compact signature using direct 2D adjacency access
      signature_length = 0
      do j = 1, array_trees%adj_counts1(item_value)
         adj_atom = array_trees%adj_lists1(item_value, j)
         part_ref_idx = array_trees%itemdir1_entries(read_link_idx, adj_atom)
         if (part_ref_idx /= 0) then
            signature_length = signature_length + 1
            signature_array(signature_length) = part_ref_idx
         end if
      end do

      ! Get relative index directly - no search needed
      target_relative_idx = find_child_part_array(array_trees, part_idx)

      ! Get absolute part index from relative index
      target_part_idx = array_trees%parts(part_idx)%child_indices(target_relative_idx)

      ! Add item to target child using relative index directly
      items1_trackers(target_relative_idx) = items1_trackers(target_relative_idx) + 1
      target_idx = array_trees%parts(target_part_idx)%items1_offset + items1_trackers(target_relative_idx)
      array_trees%item1_values(target_idx) = item_value

      ! Update itemdir using 2D array - no offset calculation needed
      array_trees%itemdir1_entries(write_link_idx, item_value) = target_part_idx
   end do

   ! Process second molecule items with ULTRA-OPTIMIZED signature generation using direct 2D adjacency
   do i = 1, items2_count
      item_value = array_trees%item2_values(items2_offset + i)

      ! ULTRA-OPTIMIZED: Generate compact signature using direct 2D adjacency access
      signature_length = 0
      do j = 1, array_trees%adj_counts2(item_value)
         adj_atom = array_trees%adj_lists2(item_value, j)
         part_ref_idx = array_trees%itemdir2_entries(read_link_idx, adj_atom)
         if (part_ref_idx /= 0) then
            signature_length = signature_length + 1
            signature_array(signature_length) = part_ref_idx
         end if
      end do

      ! Get relative index directly - no search needed
      target_relative_idx = find_child_part_array(array_trees, part_idx)

      ! Get absolute part index from relative index
      target_part_idx = array_trees%parts(part_idx)%child_indices(target_relative_idx)

      ! Add item to target child using relative index directly
      items2_trackers(target_relative_idx) = items2_trackers(target_relative_idx) + 1
      target_idx = array_trees%parts(target_part_idx)%items2_offset + items2_trackers(target_relative_idx)
      array_trees%item2_values(target_idx) = item_value

      ! Update itemdir using 2D array - no offset calculation needed
      array_trees%itemdir2_entries(write_link_idx, item_value) = target_part_idx
   end do

   ! Check if the split created any leaf parts and accumulate distance contribution
   total_distance = total_distance + check_leaf_parts_for_squared_distance(coords1, coords2, array_trees, part_idx)
end subroutine

subroutine recompute_nextlevel_mnas_array(coords1, coords2, array_trees, link_idx, total_distance)
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: link_idx
   real(rk), intent(inout) :: total_distance
   integer :: next_link_idx, i, part_idx
   integer :: num_parts, partref_offset

   next_link_idx = link_idx + 1
   num_parts = array_trees%links(link_idx)%num_parts
   partref_offset = array_trees%links(link_idx)%partref_offset

   do i = 1, num_parts
      part_idx = array_trees%partref_entries(partref_offset + i)
      call resplit_part_mna_array(coords1, coords2, array_trees, part_idx, link_idx, next_link_idx, total_distance)
   end do
end subroutine

subroutine assign_and_recompute_mnas_array(coords1, coords2, array_trees, split_part_idx, child_branch_idx, &
      first_link_idx, chosen_item2_idx, total_distance)
   ! Combined procedure: DFS assignment + MNA recomputation
   ! Makes assignment (first item1 with chosen_item2_idx-th item2) then recomputes MNAs for the branch
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: split_part_idx, child_branch_idx, first_link_idx, chosen_item2_idx
   real(rk), intent(inout) :: total_distance
   integer :: child_part1, child_part2, first_item1, chosen_item2, i, item_value, target_idx
   integer :: items1_offset, items1_count, items2_offset, items2_count
   integer :: link_idx, num_links, link_offset

   ! === PART 1: DFS ASSIGNMENT ===

   ! Extract commonly used offsets and values
   items1_offset = array_trees%parts(split_part_idx)%items1_offset
   items1_count = array_trees%parts(split_part_idx)%items1_count
   items2_offset = array_trees%parts(split_part_idx)%items2_offset
   items2_count = array_trees%parts(split_part_idx)%items2_count

   ! Get child parts using direct array access
   child_part1 = array_trees%parts(split_part_idx)%child_indices(1)
   child_part2 = array_trees%parts(split_part_idx)%child_indices(2)

   ! Get first item1 and chosen item2
   first_item1 = array_trees%item1_values(items1_offset + 1)
   chosen_item2 = array_trees%item2_values(items2_offset + chosen_item2_idx)

   ! Assign chosen items to first child (direct placement)
   array_trees%item1_values(array_trees%parts(child_part1)%items1_offset + 1) = first_item1
   array_trees%item2_values(array_trees%parts(child_part1)%items2_offset + 1) = chosen_item2

   ! Update itemdir using 2D arrays - no offset calculation needed
   array_trees%itemdir1_entries(first_link_idx, first_item1) = child_part1
   array_trees%itemdir2_entries(first_link_idx, chosen_item2) = child_part1

   ! Copy remaining items1 to second child (skip the first item)
   target_idx = array_trees%parts(child_part2)%items1_offset
   do i = 2, items1_count
      item_value = array_trees%item1_values(items1_offset + i)
      target_idx = target_idx + 1
      array_trees%item1_values(target_idx) = item_value
      array_trees%itemdir1_entries(first_link_idx, item_value) = child_part2
   end do

   ! Copy remaining items2 to second child (skip the chosen item)
   target_idx = array_trees%parts(child_part2)%items2_offset
   do i = 1, items2_count
      if (i /= chosen_item2_idx) then
         item_value = array_trees%item2_values(items2_offset + i)
         target_idx = target_idx + 1
         array_trees%item2_values(target_idx) = item_value
         array_trees%itemdir2_entries(first_link_idx, item_value) = child_part2
      end if
   end do

   ! Accumulate distance from assignment (leaf parts created by the split)
   total_distance = total_distance + check_leaf_parts_for_squared_distance(coords1, coords2, array_trees, split_part_idx)

   ! === PART 2: MNA RECOMPUTATION ===

   num_links = array_trees%chains(child_branch_idx)%num_links
   link_offset = array_trees%chains(child_branch_idx)%link_offset

   do i = 1, num_links
      link_idx = link_offset + i
      call recompute_nextlevel_mnas_array(coords1, coords2, array_trees, link_idx, total_distance)
   end do
end subroutine

recursive subroutine redistribute_items_dfs_recursive(coords1, coords2, array_trees, branch_idx, total_distance)
   ! DFS exploration of all assignment possibilities - accumulates sum of best distances from all child branches
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: branch_idx
   real(rk), intent(inout) :: total_distance

   integer :: child_branch_idx, first_link_idx, split_part_idx, i, items2_count, j
   integer, allocatable :: temp_itemdir1_entries(:,:), temp_itemdir2_entries(:,:)
   real(rk) :: combination_distance
   real(rk) :: best_branch_distance

   ! Check if this is a leaf level (no more child branches)
   if (array_trees%chains(branch_idx)%num_children == 0) then
      combination_count = combination_count + 1
      return
   end if

   ! Process each child branch independently
   do i = 1, array_trees%chains(branch_idx)%num_children
      child_branch_idx = array_trees%chains(branch_idx)%child_indices(i)
      first_link_idx = array_trees%chains(child_branch_idx)%link_offset + 1
      split_part_idx = array_trees%chains(child_branch_idx)%split_part_idx

      items2_count = array_trees%parts(split_part_idx)%items2_count
      best_branch_distance = huge(1.0_rk)  ! Best distance for this specific branch

      ! Try pairing first item1 with each item2 to find best assignment for this branch
      do j = 1, items2_count
         ! Save current state for backtracking (automatic allocation)
         temp_itemdir1_entries = array_trees%itemdir1_entries
         temp_itemdir2_entries = array_trees%itemdir2_entries
         combination_distance = 0.0_rk  ! Start fresh for this combination

         ! Make assignment and recompute MNAs in one combined operation
         call assign_and_recompute_mnas_array(coords1, coords2, array_trees, split_part_idx, child_branch_idx, &
            first_link_idx, j, combination_distance)

         ! Recursively explore subtree and accumulate distance
         call redistribute_items_dfs_recursive(coords1, coords2, array_trees, child_branch_idx, combination_distance)

         ! Update best distance for this branch if this assignment is better
         if (combination_distance < best_branch_distance) then
            best_branch_distance = combination_distance
         end if

         ! Backtrack: restore state for next iteration (automatic allocation)
         array_trees%itemdir1_entries = temp_itemdir1_entries
         array_trees%itemdir2_entries = temp_itemdir2_entries
      end do

      ! Add the best distance from this branch to the total
      total_distance = total_distance + best_branch_distance
   end do
end subroutine

subroutine redistribute_items_array(coords1, coords2, array_trees, branch_idx, final_total_distance)
   ! DFS exploration wrapper - finds optimal assignment among all possibilities
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: branch_idx
   real(rk), intent(out), optional :: final_total_distance

   real(rk) :: optimal_distance

   ! Initialize DFS exploration variables
   combination_count = 0
   optimal_distance = 0.0_rk

   ! Clear all itemdir entries using 2D array operations
   array_trees%itemdir1_entries = 0
   array_trees%itemdir2_entries = 0

   write(stderr, '(A)') "=== Starting DFS exploration of all assignment possibilities ==="

   ! Perform DFS exploration to find optimal assignment
   call redistribute_items_dfs_recursive(coords1, coords2, array_trees, branch_idx, optimal_distance)

   ! Report final results
   write(stderr, '(A)') repeat("=", 60)
   write(stderr, '(A,I0)') "Assignment combinations probed: ", combination_count
   write(stderr, '(A,F10.4)') "Optimal total squared distance: ", optimal_distance
   write(stderr, '(A)') repeat("=", 60)

   ! Return optimal distance if requested
   if (present(final_total_distance)) then
      final_total_distance = optimal_distance
   end if
end subroutine

end module
