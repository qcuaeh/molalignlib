module assigntree_distribute
use parameters
use random
use derived_types
use lcrs_frame
use permutation
use spatial_transforms
implicit none
private
public distribute_items_dfs
public distribute_items_random

! Maximum possible number of children for a part
integer, parameter :: MAX_CHILDREN = 8

! Module-level signature workspace to eliminate allocations
integer :: signature_array(MAX_COORD)
integer :: signature_length

! DFS exploration variables
integer :: combination_count

contains

function signature_equivalence_array(assign_frame, part_idx) result(equiv)
   ! OPTIMIZED: Fast path for length-1 signatures (most common case)
   type(array_trees_t), intent(in) :: assign_frame
   integer, intent(in) :: part_idx
   logical :: equiv
   integer :: signature_frequencies, i, j

   if (signature_length /= assign_frame%partree(part_idx)%signature_length) then
      equiv = .false.
      return
   end if

   ! FAST PATH: Direct comparison for length-1 signatures (most common)
   if (assign_frame%partree(part_idx)%signature_length == 1) then
      equiv = (signature_array(1) == assign_frame%partree(part_idx)%signature_values(1))
      return
   end if

   ! GENERIC PATH
   do i = 1, assign_frame%partree(part_idx)%signature_unique_count
      signature_frequencies = 0

      ! Count matches in target signature
      do j = 1, signature_length
         if (signature_array(j) == assign_frame%partree(part_idx)%signature_values(i)) then
            signature_frequencies = signature_frequencies + 1
         end if
      end do

      if (signature_frequencies /= assign_frame%partree(part_idx)%signature_frequencies(i)) then
         equiv = .false.
         return
      end if
   end do

   equiv = .true.
end function

function find_child_part_array(assign_frame, parent_idx) result(child_relative_idx)
   ! OPTIMIZED: Assumes exactly 2 children - if signature doesn't match first, it must match second
   type(array_trees_t), intent(in) :: assign_frame
   integer, intent(in) :: parent_idx
   integer :: child_relative_idx
   integer :: num_children, child_idx, i

   ! Get first child index directly
   num_children = assign_frame%partree(parent_idx)%num_children

   ! Check if signature matches first or second child
   do i = 1, num_children
      child_idx = assign_frame%partree(parent_idx)%child_indices(i)
      if (signature_equivalence_array(assign_frame, child_idx)) then
         child_relative_idx = i
         return
      end if
   end do

   error stop 'Part signature does not match any child'
end function

subroutine collect_leaf_assignments(assign_frame, part_idx, subperm)
   ! Collect assignment pairs from leaf parts into assignment
   type(array_trees_t), intent(in) :: assign_frame
   integer, intent(in) :: part_idx
   type(subperm_t), intent(inout) :: subperm
   integer :: i, child_idx, item1_idx, item2_idx

   ! Check all children of this part
   do i = 1, assign_frame%partree(part_idx)%num_children
      child_idx = assign_frame%partree(part_idx)%child_indices(i)

      ! If this child is a leaf (no children), collect its assignment
      if (assign_frame%partree(child_idx)%num_children == 0) then
         ! Verify this is a proper leaf part with exactly one item from each molecule
         if (assign_frame%partree(child_idx)%items1_count == 1 .and. &
             assign_frame%partree(child_idx)%items2_count == 1) then

            ! Get the assigned items
            item1_idx = assign_frame%item1_values(assign_frame%partree(child_idx)%items1_offset + 1)
            item2_idx = assign_frame%item2_values(assign_frame%partree(child_idx)%items2_offset + 1)

            ! Add item pair to assignment
            call subperm_add(subperm, item1_idx, item2_idx)
         end if
      end if
   end do
end subroutine

subroutine resplit_part_mna(assign_frame, part_idx, read_link_idx, write_link_idx, subperm)
! ULTRA-OPTIMIZED: Array-based version with direct 2D adjacency access for maximum performance
! UPDATED: Now collects assignment pairs from newly created leaf parts into subperm
   type(array_trees_t), intent(inout) :: assign_frame
   integer, intent(in) :: part_idx, read_link_idx, write_link_idx
   type(subperm_t), intent(inout) :: subperm
   integer, dimension(MAX_CHILDREN) :: items1_trackers, items2_trackers
   integer :: i, j, target_relative_idx, target_part_idx, item_value, target_idx, part_ref_idx, adj_atom
   integer :: items1_offset, items1_count, items2_offset, items2_count
   integer :: num_children

   ! Extract commonly used values for readability
   items1_offset = assign_frame%partree(part_idx)%items1_offset
   items1_count = assign_frame%partree(part_idx)%items1_count
   items2_offset = assign_frame%partree(part_idx)%items2_offset
   items2_count = assign_frame%partree(part_idx)%items2_count
   num_children = assign_frame%partree(part_idx)%num_children

   ! INITIALIZATION: Reset trackers using relative indices (1 to num_children)
   do i = 1, num_children
      items1_trackers(i) = 0
      items2_trackers(i) = 0
   end do

   ! Process first molecule items with ULTRA-OPTIMIZED signature generation using direct 2D adjacency
   do i = 1, items1_count
      item_value = assign_frame%item1_values(items1_offset + i)

      ! ULTRA-OPTIMIZED: Generate compact signature using direct 2D adjacency access
      signature_length = 0
      do j = 1, assign_frame%adj_counts1(item_value)
         adj_atom = assign_frame%adj_lists1(item_value, j)
         part_ref_idx = assign_frame%itemdir1_entries(read_link_idx, adj_atom)
         if (part_ref_idx /= 0) then
            signature_length = signature_length + 1
            signature_array(signature_length) = part_ref_idx
         end if
      end do

      ! Get relative index directly - no search needed
      target_relative_idx = find_child_part_array(assign_frame, part_idx)

      ! Get absolute part index from relative index
      target_part_idx = assign_frame%partree(part_idx)%child_indices(target_relative_idx)

      ! Add item to target child using relative index directly
      items1_trackers(target_relative_idx) = items1_trackers(target_relative_idx) + 1
      target_idx = assign_frame%partree(target_part_idx)%items1_offset + items1_trackers(target_relative_idx)
      assign_frame%item1_values(target_idx) = item_value

      ! Update itemdir using 2D array - no offset calculation needed
      assign_frame%itemdir1_entries(write_link_idx, item_value) = target_part_idx
   end do

   ! Process second molecule items with ULTRA-OPTIMIZED signature generation using direct 2D adjacency
   do i = 1, items2_count
      item_value = assign_frame%item2_values(items2_offset + i)

      ! ULTRA-OPTIMIZED: Generate compact signature using direct 2D adjacency access
      signature_length = 0
      do j = 1, assign_frame%adj_counts2(item_value)
         adj_atom = assign_frame%adj_lists2(item_value, j)
         part_ref_idx = assign_frame%itemdir2_entries(read_link_idx, adj_atom)
         if (part_ref_idx /= 0) then
            signature_length = signature_length + 1
            signature_array(signature_length) = part_ref_idx
         end if
      end do

      ! Get relative index directly - no search needed
      target_relative_idx = find_child_part_array(assign_frame, part_idx)

      ! Get absolute part index from relative index
      target_part_idx = assign_frame%partree(part_idx)%child_indices(target_relative_idx)

      ! Add item to target child using relative index directly
      items2_trackers(target_relative_idx) = items2_trackers(target_relative_idx) + 1
      target_idx = assign_frame%partree(target_part_idx)%items2_offset + items2_trackers(target_relative_idx)
      assign_frame%item2_values(target_idx) = item_value

      ! Update itemdir using 2D array - no offset calculation needed
      assign_frame%itemdir2_entries(write_link_idx, item_value) = target_part_idx
   end do

   ! Collect assignment pairs from newly created leaf parts into subperm
   call collect_leaf_assignments(assign_frame, part_idx, subperm)
end subroutine

subroutine recompute_mna_partition(assign_frame, link_idx, subperm)
   type(array_trees_t), intent(inout) :: assign_frame
   integer, intent(in) :: link_idx
   type(subperm_t), intent(inout) :: subperm
   integer :: next_link_idx, i, part_idx
   integer :: num_parts, partref_offset

   next_link_idx = link_idx + 1
   num_parts = assign_frame%chain(link_idx)%num_parts
   partref_offset = assign_frame%chain(link_idx)%partref_offset

   do i = 1, num_parts
      part_idx = assign_frame%partref_entries(partref_offset + i)
      call resplit_part_mna(assign_frame, part_idx, link_idx, next_link_idx, subperm)
   end do
end subroutine

subroutine distribute_part_items(assign_frame, split_part_idx, child_branch_idx, &
      first_link_idx, chosen_item1_idx, chosen_item2_idx, subperm)
   ! Combined procedure: assignment + MNA recomputation
   ! Makes assignment (chosen_item1_idx-th item1 with chosen_item2_idx-th item2) then recomputes MNAs for the branch
   ! UPDATED: Now accepts both item1 and item2 indices to match original random assignment behavior
   type(array_trees_t), intent(inout) :: assign_frame
   integer, intent(in) :: split_part_idx, child_branch_idx, first_link_idx, chosen_item1_idx, chosen_item2_idx
   type(subperm_t), intent(inout) :: subperm
   integer :: child_part1, child_part2, chosen_item1, chosen_item2, i, item_value, target_idx
   integer :: items1_offset, items1_count, items2_offset, items2_count
   integer :: link_idx, num_links, link_offset

   ! === PART 1: PAIR ASSIGNMENT ===

   ! Extract commonly used offsets and values
   items1_offset = assign_frame%partree(split_part_idx)%items1_offset
   items1_count = assign_frame%partree(split_part_idx)%items1_count
   items2_offset = assign_frame%partree(split_part_idx)%items2_offset
   items2_count = assign_frame%partree(split_part_idx)%items2_count

   ! Get child parts using direct array access
   child_part1 = assign_frame%partree(split_part_idx)%child_indices(1)
   child_part2 = assign_frame%partree(split_part_idx)%child_indices(2)

   ! Get chosen items based on provided indices
   chosen_item1 = assign_frame%item1_values(items1_offset + chosen_item1_idx)
   chosen_item2 = assign_frame%item2_values(items2_offset + chosen_item2_idx)

   ! Assign chosen items to first child (direct placement)
   assign_frame%item1_values(assign_frame%partree(child_part1)%items1_offset + 1) = chosen_item1
   assign_frame%item2_values(assign_frame%partree(child_part1)%items2_offset + 1) = chosen_item2

   ! Update itemdir using 2D arrays - no offset calculation needed
   assign_frame%itemdir1_entries(first_link_idx, chosen_item1) = child_part1
   assign_frame%itemdir2_entries(first_link_idx, chosen_item2) = child_part1

   ! Copy remaining items1 to second child (skip the chosen item)
   target_idx = assign_frame%partree(child_part2)%items1_offset
   do i = 1, items1_count
      if (i /= chosen_item1_idx) then
         item_value = assign_frame%item1_values(items1_offset + i)
         target_idx = target_idx + 1
         assign_frame%item1_values(target_idx) = item_value
         assign_frame%itemdir1_entries(first_link_idx, item_value) = child_part2
      end if
   end do

   ! Copy remaining items2 to second child (skip the chosen item)
   target_idx = assign_frame%partree(child_part2)%items2_offset
   do i = 1, items2_count
      if (i /= chosen_item2_idx) then
         item_value = assign_frame%item2_values(items2_offset + i)
         target_idx = target_idx + 1
         assign_frame%item2_values(target_idx) = item_value
         assign_frame%itemdir2_entries(first_link_idx, item_value) = child_part2
      end if
   end do

   ! Collect assignment pairs from assignment (leaf parts created by the split)
   call collect_leaf_assignments(assign_frame, split_part_idx, subperm)

   ! === PART 2: SCNA RECOMPUTATION ===

   num_links = assign_frame%assigntree(child_branch_idx)%num_links
   link_offset = assign_frame%assigntree(child_branch_idx)%link_offset

   do i = 1, num_links
      link_idx = link_offset + i
      call recompute_mna_partition(assign_frame, link_idx, subperm)
   end do
end subroutine

recursive subroutine distribute_items_dfs_recursive(coords1, coords2, assign_frame, branch_idx, treeperm)
   ! DFS exploration of all assignment possibilities - finds permutation that minimizes total distance
   ! OPTIMIZED: Reduces allocations by reusing arrays within branch scope, but maintains isolation between branches
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)  ! coordinates needed for distance calculation
   type(array_trees_t), intent(inout) :: assign_frame
   integer, intent(in) :: branch_idx
   type(subperm_t), intent(inout) :: treeperm

   integer :: child_branch_idx, first_link_idx, split_part_idx, i, items2_count, j
   integer :: link_idx, branch_link_offset, branch_num_links
   type(subperm_t) :: best_branchperm, branchperm
   real(rk) :: branch_distance, best_branch_distance
   integer :: num_atoms

   num_atoms = assign_frame%num_atoms1

   ! Check if this is a leaf level (no more child branches)
   if (assign_frame%assigntree(branch_idx)%num_children == 0) then
      combination_count = combination_count + 1
      return
   end if

   ! Process each child branch independently
   do i = 1, assign_frame%assigntree(branch_idx)%num_children
      child_branch_idx = assign_frame%assigntree(branch_idx)%child_indices(i)
      first_link_idx = assign_frame%assigntree(child_branch_idx)%link_offset + 1
      split_part_idx = assign_frame%assigntree(child_branch_idx)%split_part_idx
      branch_link_offset = assign_frame%assigntree(child_branch_idx)%link_offset
      branch_num_links = assign_frame%assigntree(child_branch_idx)%num_links

      items2_count = assign_frame%partree(split_part_idx)%items2_count
      best_branch_distance = huge(1.0_rk)  ! Best distance for this specific branch

      call subperm_init(best_branchperm, num_atoms)
      call subperm_init(branchperm, num_atoms)

      ! Try pairing first item1 with each item2 to find best assignment for this branch
      do j = 1, items2_count
         ! Reset branch assignment for this iteration
         branchperm%size = 0

         ! Make assignment and recompute MNAs in one combined operation (using first item1, index=1)
         call distribute_part_items(assign_frame, split_part_idx, child_branch_idx, first_link_idx, 1, j, &
            branchperm)

         ! Recursively explore subtree and collect child permutation
         call distribute_items_dfs_recursive(coords1, coords2, assign_frame, child_branch_idx, branchperm)

         ! Calculate partial distance for this branch
         branch_distance = total_sqdist(branchperm, coords1, coords2)

         ! Update best distance for this branch if this assignment is better
         if (branch_distance < best_branch_distance) then
            best_branch_distance = branch_distance
            best_branchperm = branchperm
         end if

         ! Reset state for next iteration - only reset links used by this branch
         do link_idx = branch_link_offset + 1, branch_link_offset + branch_num_links
            assign_frame%itemdir1_entries(link_idx, :) = 0
            assign_frame%itemdir2_entries(link_idx, :) = 0
         end do
      end do

      ! Update the optimal assignment with the best assigment from this branch
      call subperm_merge(treeperm, best_branchperm)
   end do
end subroutine

recursive subroutine distribute_items_random_recursive(coords1, coords2, assign_frame, branch_idx, treeperm)
   ! Random exploration - generates one assignment randomly using same traversal order as random module
   ! Similar to distribute_items_dfs_recursive but picks one random assignment instead of exploring all
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: assign_frame
   integer, intent(in) :: branch_idx
   type(subperm_t), intent(inout) :: treeperm

   integer :: rand_idx1, rand_idx2
   integer :: link_idx, branch_link_offset, branch_num_links
   integer :: child_branch_idx, first_link_idx, split_part_idx
   integer :: items1_count, items2_count, i

   ! Check if this is a leaf level (no more child branches)
   if (assign_frame%assigntree(branch_idx)%num_children == 0) then
      return
   end if

   ! Process each child branch using same traversal order as random module
   do i = 1, assign_frame%assigntree(branch_idx)%num_children
      child_branch_idx = assign_frame%assigntree(branch_idx)%child_indices(i)
      first_link_idx = assign_frame%assigntree(child_branch_idx)%link_offset + 1
      split_part_idx = assign_frame%assigntree(child_branch_idx)%split_part_idx
      branch_link_offset = assign_frame%assigntree(child_branch_idx)%link_offset
      branch_num_links = assign_frame%assigntree(child_branch_idx)%num_links

      items1_count = assign_frame%partree(split_part_idx)%items1_count
      items2_count = assign_frame%partree(split_part_idx)%items2_count

      ! Generate random choices for both item1 and item2 indices (matching original random module)
      rand_idx1 = random_uniform_integer(1, items1_count)
      rand_idx2 = random_uniform_integer(1, items2_count)

      ! Make random assignment and recompute MNAs using existing procedure
      call distribute_part_items(assign_frame, split_part_idx, child_branch_idx, first_link_idx, &
         rand_idx1, rand_idx2, treeperm)

      ! Recursively explore child branch
      call distribute_items_random_recursive(coords1, coords2, assign_frame, child_branch_idx, treeperm)

      ! Reset state for next iteration - only reset links used by this branch
      do link_idx = branch_link_offset + 1, branch_link_offset + branch_num_links
         assign_frame%itemdir1_entries(link_idx, :) = 0
         assign_frame%itemdir2_entries(link_idx, :) = 0
      end do
   end do
end subroutine

subroutine distribute_items_random(coords1, coords2, assign_frame, treeperm)
   ! Random exploration wrapper - generates one random assignment
   ! Similar to distribute_items_dfs but generates random assignment instead of optimal
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: assign_frame
   type(subperm_t), intent(out) :: treeperm
   ! Local variables
   real(rk) :: total_distance

!   call random_init(.true., .true.)

   ! Initialize random assignment
   call subperm_init(treeperm, assign_frame%num_atoms1)

   ! Initialize assignment with preassigned pairs
   call collect_leaf_assignments(assign_frame, 1, treeperm)

   ! Perform random exploration to generate one assignment (starting from root chain at index 1)
   call distribute_items_random_recursive(coords1, coords2, assign_frame, 1, treeperm)

   ! Calculate distance from the random permutation array
   total_distance = total_sqdist(treeperm, coords1, coords2)
   write(stderr, '(A)') repeat("=", 60)
   write(stderr, '(A,I0,A,I0,A)') "Atoms assigned: ", treeperm%size, " out of ", &
         assign_frame%num_atoms1, " total atoms"
   write(stderr, '(A,F10.4)') "Random assignment total squared distance: ", total_distance
   write(stderr, '(A)') repeat("=", 60)
!   call check_subperm(treeperm)
end subroutine

subroutine distribute_items_dfs(coords1, coords2, assign_frame, treeperm)
   ! DFS exploration wrapper - finds optimal assignment among all possibilities
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: assign_frame
   type(subperm_t), intent(out) :: treeperm
   ! Local variables
   integer :: num_atoms

   num_atoms = assign_frame%num_atoms1

   ! Initialize optimal assignment
   call subperm_init(treeperm, num_atoms)

   ! Initialize optimal assignment with preassigned pairs
   call collect_leaf_assignments(assign_frame, 1, treeperm)

   ! Initialize DFS exploration variables
   combination_count = 0

   ! Perform DFS exploration to find optimal assignment (starting from root chain at index 1)
   call distribute_items_dfs_recursive(coords1, coords2, assign_frame, 1, treeperm)

!block
!   integer :: assigned_count
!   real(rk) :: total_distance
!   assigned_count = treeperm%size
!   total_distance = total_sqdist(treeperm, coords1, coords2)
!   write(stderr, '(A)') repeat("=", 60)
!   write(stderr, '(A,I0)') "Assignment combinations probed: ", combination_count
!   write(stderr, '(A,I0,A,I0,A)') "Atoms assigned: ", assigned_count, " out of ", num_atoms, " total atoms"
!   write(stderr, '(A,F10.4)') "Optimal total squared distance: ", total_distance
!   write(stderr, '(A)') repeat("=", 60)
!end block
end subroutine

end module
