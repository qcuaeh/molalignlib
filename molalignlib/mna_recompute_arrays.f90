module mna_recompute_arrays
use parameters
use molecule
use array_trees
implicit none

! Maximum possible number of children for a part
integer, parameter :: MAX_CHILDREN = 2

! Module-level signature workspace to eliminate allocations
integer :: signature_array(MAX_COORD)
integer :: signature_length

! Module-level total squared distance tracking variables
real(rk) :: total_squared_distance
integer :: total_assigned_pairs

! DFS exploration variables
real(rk) :: best_total_distance
integer :: exploration_count
integer, allocatable :: best_assignment_state(:)  ! Store best itemdir state

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

subroutine calculate_leaf_squared_distance_contribution(mol1, mol2, array_trees, part_idx)
   ! Calculate squared distance contribution for a leaf part (items1_count == items2_count == 1)
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   integer :: item1_idx, item2_idx
   real(rk) :: squared_distance
   real(rk) :: dx, dy, dz

   ! Verify this is a leaf part with exactly one item from each molecule
   if (array_trees%parts(part_idx)%items1_count /= 1 .or. &
       array_trees%parts(part_idx)%items2_count /= 1) then
      write(stderr, '(A,I0,A,I0,A,I0)') 'Warning: Part ', part_idx, &
         ' is not a proper leaf (items1=', array_trees%parts(part_idx)%items1_count, &
         ', items2=', array_trees%parts(part_idx)%items2_count, ')'
      return
   end if

   ! Get the assigned items
   item1_idx = array_trees%item1_values(array_trees%parts(part_idx)%items1_offset + 1)
   item2_idx = array_trees%item2_values(array_trees%parts(part_idx)%items2_offset + 1)

   ! Calculate squared distance between assigned atoms
   dx = mol1%atoms(item1_idx)%coords(1) - mol2%atoms(item2_idx)%coords(1)
   dy = mol1%atoms(item1_idx)%coords(2) - mol2%atoms(item2_idx)%coords(2)
   dz = mol1%atoms(item1_idx)%coords(3) - mol2%atoms(item2_idx)%coords(3)

   squared_distance = dx*dx + dy*dy + dz*dz

   ! Update progressive total squared distance calculation
   total_squared_distance = total_squared_distance + squared_distance
   total_assigned_pairs = total_assigned_pairs + 1
end subroutine

subroutine check_leaf_parts_for_squared_distance(mol1, mol2, array_trees, part_idx)
   ! Check if child parts are leaves and calculate squared distance contributions
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   integer :: i, child_idx

   ! Check all children of this part
   do i = 1, array_trees%parts(part_idx)%num_children
      child_idx = array_trees%parts(part_idx)%child_indices(i)

      ! If this child is a leaf (no children), calculate squared distance contribution
      if (array_trees%parts(child_idx)%num_children == 0) then
         call calculate_leaf_squared_distance_contribution(mol1, mol2, array_trees, child_idx)
      end if
   end do
end subroutine

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

subroutine resplit_part_mna_array(mol1, mol2, array_trees, part_idx, read_link_idx, write_link_idx)
! Array-based version of resplit_part_mna with inlined signature generation
! UPDATED: Now uses relative indices directly without searching - more efficient
! UPDATED: Added RMSD calculation for newly created leaf parts
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx, read_link_idx, write_link_idx
   integer :: items1_trackers(MAX_CHILDREN), items2_trackers(MAX_CHILDREN)
   integer :: i, j, target_relative_idx, target_part_idx, item_value, target_idx, part_ref_idx
   integer :: items1_offset, items1_count, items2_offset, items2_count
   integer :: itemdir1_offset, itemdir2_offset, read_itemdir1_offset, read_itemdir2_offset
   integer :: num_children

   ! Extract commonly used values for readability
   items1_offset = array_trees%parts(part_idx)%items1_offset
   items1_count = array_trees%parts(part_idx)%items1_count
   items2_offset = array_trees%parts(part_idx)%items2_offset
   items2_count = array_trees%parts(part_idx)%items2_count
   itemdir1_offset = array_trees%links(write_link_idx)%itemdir1_offset
   itemdir2_offset = array_trees%links(write_link_idx)%itemdir2_offset
   read_itemdir1_offset = array_trees%links(read_link_idx)%itemdir1_offset
   read_itemdir2_offset = array_trees%links(read_link_idx)%itemdir2_offset
   num_children = array_trees%parts(part_idx)%num_children

   ! INITIALIZATION: Reset redistribution counters using relative indices (1 to num_children)
   do i = 1, num_children
      items1_trackers(i) = 0
      items2_trackers(i) = 0
   end do

   ! Process first molecule items with inlined signature generation
   do i = 1, items1_count
      item_value = array_trees%item1_values(items1_offset + i)

      ! INLINED: Generate compact signature from itemdir1
      signature_length = 0
      do j = 1, size(mol1%atoms(item_value)%adjlist)
         part_ref_idx = array_trees%itemdir_entries(read_itemdir1_offset + mol1%atoms(item_value)%adjlist(j))
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

      ! Update itemdir
      array_trees%itemdir_entries(itemdir1_offset + item_value) = target_part_idx
   end do

   ! Process second molecule items with inlined signature generation
   do i = 1, items2_count
      item_value = array_trees%item2_values(items2_offset + i)

      ! INLINED: Generate compact signature from itemdir2
      signature_length = 0
      do j = 1, size(mol2%atoms(item_value)%adjlist)
         part_ref_idx = array_trees%itemdir_entries(read_itemdir2_offset + mol2%atoms(item_value)%adjlist(j))
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

      ! Update itemdir
      array_trees%itemdir_entries(itemdir2_offset + item_value) = target_part_idx
   end do

   ! Check if the split created any leaf parts and calculate squared distance contributions
   call check_leaf_parts_for_squared_distance(mol1, mol2, array_trees, part_idx)
end subroutine

subroutine recompute_nextlevel_mnas_array(mol1, mol2, array_trees, link_idx)
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: link_idx
   integer :: next_link_idx, i, part_idx
   integer :: num_parts, partref_offset

   next_link_idx = link_idx + 1
   num_parts = array_trees%links(link_idx)%num_parts
   partref_offset = array_trees%links(link_idx)%partref_offset

   do i = 1, num_parts
      part_idx = array_trees%partref_entries(partref_offset + i)
      call resplit_part_mna_array(mol1, mol2, array_trees, part_idx, link_idx, next_link_idx)
   end do
end subroutine

subroutine recompute_consistent_mnas_array(mol1, mol2, array_trees, branch_idx)
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: branch_idx
   integer :: i, link_idx
   integer :: num_links, link_offset

   num_links = array_trees%chains(branch_idx)%num_links
   link_offset = array_trees%chains(branch_idx)%link_offset

   do i = 1, num_links
      link_idx = link_offset + i
      call recompute_nextlevel_mnas_array(mol1, mol2, array_trees, link_idx)
   end do
end subroutine

subroutine resplit_part_dfs_assignment(mol1, mol2, array_trees, part_idx, write_link_idx, chosen_item2_idx)
   ! DFS assignment version: assigns first item1 with chosen_item2_idx-th item2
   ! UPDATED: Added RMSD calculation for newly created leaf parts
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx, write_link_idx, chosen_item2_idx
   integer :: child_part1, child_part2, first_item1, chosen_item2, i, item_value, target_idx
   integer :: items1_offset, items1_count, items2_offset, items2_count
   integer :: itemdir1_offset, itemdir2_offset

   ! Extract commonly used offsets and values
   items1_offset = array_trees%parts(part_idx)%items1_offset
   items1_count = array_trees%parts(part_idx)%items1_count
   items2_offset = array_trees%parts(part_idx)%items2_offset
   items2_count = array_trees%parts(part_idx)%items2_count
   itemdir1_offset = array_trees%links(write_link_idx)%itemdir1_offset
   itemdir2_offset = array_trees%links(write_link_idx)%itemdir2_offset

   ! Get child parts using direct array access
   child_part1 = array_trees%parts(part_idx)%child_indices(1)
   child_part2 = array_trees%parts(part_idx)%child_indices(2)

   ! Get first item1 and chosen item2
   first_item1 = array_trees%item1_values(items1_offset + 1)
   chosen_item2 = array_trees%item2_values(items2_offset + chosen_item2_idx)

   ! Assign chosen items to first child (direct placement)
   array_trees%item1_values(array_trees%parts(child_part1)%items1_offset + 1) = first_item1
   array_trees%item2_values(array_trees%parts(child_part1)%items2_offset + 1) = chosen_item2
   array_trees%itemdir_entries(itemdir1_offset + first_item1) = child_part1
   array_trees%itemdir_entries(itemdir2_offset + chosen_item2) = child_part1

   ! Copy remaining items1 to second child (skip the first item)
   target_idx = array_trees%parts(child_part2)%items1_offset
   do i = 2, items1_count
      item_value = array_trees%item1_values(items1_offset + i)
      target_idx = target_idx + 1
      array_trees%item1_values(target_idx) = item_value
      array_trees%itemdir_entries(itemdir1_offset + item_value) = child_part2
   end do

   ! Copy remaining items2 to second child (skip the chosen item)
   target_idx = array_trees%parts(child_part2)%items2_offset
   do i = 1, items2_count
      if (i /= chosen_item2_idx) then
         item_value = array_trees%item2_values(items2_offset + i)
         target_idx = target_idx + 1
         array_trees%item2_values(target_idx) = item_value
         array_trees%itemdir_entries(itemdir2_offset + item_value) = child_part2
      end if
   end do

   ! Check if the split created any leaf parts and calculate squared distance contributions
   call check_leaf_parts_for_squared_distance(mol1, mol2, array_trees, part_idx)
end subroutine

subroutine save_itemdir_state(array_trees, saved_state)
   ! Save current itemdir state for backtracking
   type(array_trees_t), intent(in) :: array_trees
   integer, allocatable, intent(out) :: saved_state(:)

   if (allocated(saved_state)) deallocate(saved_state)
   allocate(saved_state(size(array_trees%itemdir_entries)))
   saved_state = array_trees%itemdir_entries
end subroutine

subroutine restore_itemdir_state(array_trees, saved_state)
   ! Restore itemdir state for backtracking
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: saved_state(:)

   array_trees%itemdir_entries = saved_state
end subroutine

recursive subroutine redistribute_items_dfs_recursive(mol1, mol2, array_trees, branch_idx)
   ! DFS exploration of all assignment possibilities
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: branch_idx

   integer :: child_branch_idx, first_link_idx, split_part_idx, i, items2_count, j
   integer, allocatable :: saved_itemdir_state(:)
   real(rk) :: current_distance

   ! Process each child branch
   do i = 1, array_trees%chains(branch_idx)%num_children
      child_branch_idx = array_trees%chains(branch_idx)%child_indices(i)
      first_link_idx = array_trees%chains(child_branch_idx)%link_offset + 1
      split_part_idx = array_trees%chains(child_branch_idx)%split_part_idx

      items2_count = array_trees%parts(split_part_idx)%items2_count

      ! Try pairing first item1 with each item2
      do j = 1, items2_count
         ! Save current state for backtracking
         call save_itemdir_state(array_trees, saved_itemdir_state)

         ! Make assignment: first item1 with j-th item2
         call resplit_part_dfs_assignment(mol1, mol2, array_trees, split_part_idx, first_link_idx, j)

         ! Recompute MNAs for this assignment
         call recompute_consistent_mnas_array(mol1, mol2, array_trees, child_branch_idx)

         ! Recursively explore subtree
         call redistribute_items_dfs_recursive(mol1, mol2, array_trees, child_branch_idx)

         ! Check if this is a leaf level (no more child branches)
         if (array_trees%chains(child_branch_idx)%num_children == 0) then
            exploration_count = exploration_count + 1
            current_distance = total_squared_distance

            if (current_distance < best_total_distance) then
               best_total_distance = current_distance
               ! Save best assignment state
               if (allocated(best_assignment_state)) deallocate(best_assignment_state)
               call save_itemdir_state(array_trees, best_assignment_state)

               write(stderr, '(A,I0,A,F10.4)') "Exploration ", exploration_count, &
                  ", new best distance: ", best_total_distance
            end if
         end if

         ! Backtrack: restore state for next iteration
         call restore_itemdir_state(array_trees, saved_itemdir_state)
         total_squared_distance = 0.0_rk
         total_assigned_pairs = 0

         deallocate(saved_itemdir_state)
      end do
   end do
end subroutine

subroutine redistribute_items_array(mol1, mol2, array_trees, branch_idx, final_total_squared_distance)
   ! DFS exploration wrapper - finds optimal assignment among all possibilities
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: branch_idx
   real(rk), intent(out), optional :: final_total_squared_distance

   ! Initialize DFS exploration variables
   best_total_distance = huge(1.0_rk)  ! Start with worst possible distance
   exploration_count = 0
   total_squared_distance = 0.0_rk
   total_assigned_pairs = 0

   ! Clear all itemdir entries
   array_trees%itemdir_entries = 0

   write(stderr, '(A)') "=== Starting DFS exploration of all assignment possibilities ==="

   ! Perform DFS exploration to find optimal assignment
   call redistribute_items_dfs_recursive(mol1, mol2, array_trees, branch_idx)

   ! Restore the best assignment found
   if (allocated(best_assignment_state)) then
      call restore_itemdir_state(array_trees, best_assignment_state)

      ! Recalculate final distance with best assignment
      total_squared_distance = 0.0_rk
      total_assigned_pairs = 0
      call calculate_final_distance_recursive(mol1, mol2, array_trees, branch_idx)
   end if

   ! Report final results
   write(stderr, '(A)') repeat("=", 60)
   write(stderr, '(A,I0)') "Total assignment combinations explored: ", exploration_count
   write(stderr, '(A,I0)') "Total assigned pairs in optimal solution: ", total_assigned_pairs
   write(stderr, '(A,F10.4)') "Optimal total squared distance: ", best_total_distance
   write(stderr, '(A)') repeat("=", 60)

   ! Return optimal distance if requested
   if (present(final_total_squared_distance)) then
      final_total_squared_distance = best_total_distance
   end if
end subroutine

recursive subroutine calculate_final_distance_recursive(mol1, mol2, array_trees, branch_idx)
   ! Recalculate total distance for the final optimal assignment
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: branch_idx
   integer :: i, child_branch_idx, split_part_idx, j, child_part_idx

   ! Process each child branch
   do i = 1, array_trees%chains(branch_idx)%num_children
      child_branch_idx = array_trees%chains(branch_idx)%child_indices(i)
      split_part_idx = array_trees%chains(child_branch_idx)%split_part_idx

      ! Check leaf parts in this split part
      do j = 1, array_trees%parts(split_part_idx)%num_children
         child_part_idx = array_trees%parts(split_part_idx)%child_indices(j)
         if (array_trees%parts(child_part_idx)%num_children == 0) then
            call calculate_leaf_squared_distance_contribution(mol1, mol2, array_trees, child_part_idx)
         end if
      end do

      ! Recursively process child branches
      call calculate_final_distance_recursive(mol1, mol2, array_trees, child_branch_idx)
   end do
end subroutine

! Legacy functions for compatibility - no longer used in DFS mode
subroutine resplit_part_first_array(mol1, mol2, array_trees, part_idx, write_link_idx)
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx, write_link_idx

   ! Call DFS assignment with first item2 (index 1)
   call resplit_part_dfs_assignment(mol1, mol2, array_trees, part_idx, write_link_idx, 1)
end subroutine

subroutine resplit_part_random_array(mol1, mol2, array_trees, part_idx, write_link_idx)
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx, write_link_idx
   integer :: random_item2_idx, items2_count
   real :: random_real

   items2_count = array_trees%parts(part_idx)%items2_count

   ! Generate random index for item2
   call random_number(random_real)
   random_item2_idx = int(random_real * items2_count) + 1

   ! Call DFS assignment with random item2
   call resplit_part_dfs_assignment(mol1, mol2, array_trees, part_idx, write_link_idx, random_item2_idx)
end subroutine

end module
