module assigntree_distribute
use parameters
use random
use derived_types
use lcrs_tree_arrays
use permutation
use spatial_transforms
implicit none
private
public distribute_items_branch
public distribute_items_tree
public distribute_items_random

! Maximum possible number of children for a part
integer, parameter :: MAX_CHILDREN = 8

! Module-level signature workspace to eliminate allocations
integer :: signature_array(MAX_COORD)
integer :: signature_length

! DFS exploration variables
integer :: combination_count

contains

function signature_equivalence_array(assign_arrays, part_idx) result(equiv)
   ! OPTIMIZED: Fast path for length-1 signatures (most common case)
   type(array_trees_t), intent(in) :: assign_arrays
   integer, intent(in) :: part_idx
   logical :: equiv
   integer :: signature_frequencies, i, j

   if (signature_length /= assign_arrays%partree(part_idx)%signature_length) then
      equiv = .false.
      return
   end if

   ! FAST PATH: Direct comparison for length-1 signatures (most common)
   if (assign_arrays%partree(part_idx)%signature_length == 1) then
      equiv = (signature_array(1) == assign_arrays%partree(part_idx)%signature_values(1))
      return
   end if

   ! GENERIC PATH
   do i = 1, assign_arrays%partree(part_idx)%signature_unique_count
      signature_frequencies = 0

      ! Count matches in target signature
      do j = 1, signature_length
         if (signature_array(j) == assign_arrays%partree(part_idx)%signature_values(i)) then
            signature_frequencies = signature_frequencies + 1
         end if
      end do

      if (signature_frequencies /= assign_arrays%partree(part_idx)%signature_frequencies(i)) then
         equiv = .false.
         return
      end if
   end do

   equiv = .true.
end function

function find_child_part_array(assign_arrays, parent_idx) result(child_relative_idx)
   type(array_trees_t), intent(in) :: assign_arrays
   integer, intent(in) :: parent_idx
   integer :: child_relative_idx
   integer :: num_children, child_idx, i

   ! Get first child index directly
   num_children = assign_arrays%partree(parent_idx)%num_children

   ! Find the child that matches the signature
   do i = 1, num_children
      child_idx = assign_arrays%partree(parent_idx)%child_indices(i)
      if (signature_equivalence_array(assign_arrays, child_idx)) then
         child_relative_idx = i
         return
      end if
   end do

   ! If the signature did not match any child it is a bug
   error stop 'Part signature did not match any child'
end function

subroutine collect_leaf_assignments(assign_arrays, part_idx, subperm)
   ! Collect assignment pairs from leaf parts into assignment
   type(array_trees_t), intent(in) :: assign_arrays
   integer, intent(in) :: part_idx
   type(subperm_t), intent(inout) :: subperm
   integer :: i, child_idx, item1_idx, item2_idx

   ! Check all children of this part
   do i = 1, assign_arrays%partree(part_idx)%num_children
      child_idx = assign_arrays%partree(part_idx)%child_indices(i)

      ! If this child is a leaf (no children), collect its assignment
      if (assign_arrays%partree(child_idx)%num_children == 0) then
         ! Verify this is a proper leaf part with exactly one item from each molecule
         if (assign_arrays%partree(child_idx)%items1_count == 1 .and. &
             assign_arrays%partree(child_idx)%items2_count == 1) then

            ! Get the assigned items
            item1_idx = assign_arrays%vertidcs1(assign_arrays%partree(child_idx)%items1_offset + 1)
            item2_idx = assign_arrays%vertidcs2(assign_arrays%partree(child_idx)%items2_offset + 1)

            ! Add item pair to assignment
            call subperm_add(subperm, item1_idx, item2_idx)
         end if
      end if
   end do
end subroutine

subroutine resplit_part_mna(assign_arrays, part_idx, read_link_idx, write_link_idx, subperm)
! ULTRA-OPTIMIZED: Array-based version with direct 2D adjacency access for maximum performance
! UPDATED: Now collects assignment pairs from newly created leaf parts into subperm
   type(array_trees_t), intent(inout) :: assign_arrays
   integer, intent(in) :: part_idx, read_link_idx, write_link_idx
   type(subperm_t), intent(inout) :: subperm
   integer, dimension(MAX_CHILDREN) :: items1_trackers, items2_trackers
   integer :: i, j, target_relative_idx, target_part_idx, item_vertidx, target_idx, part_ref_idx, adj_atom
   integer :: items1_offset, items1_count, items2_offset, items2_count
   integer :: num_children

   ! Extract commonly used values for readability
   items1_offset = assign_arrays%partree(part_idx)%items1_offset
   items1_count = assign_arrays%partree(part_idx)%items1_count
   items2_offset = assign_arrays%partree(part_idx)%items2_offset
   items2_count = assign_arrays%partree(part_idx)%items2_count
   num_children = assign_arrays%partree(part_idx)%num_children

   ! INITIALIZATION: Reset trackers using relative indices (1 to num_children)
   do i = 1, num_children
      items1_trackers(i) = 0
      items2_trackers(i) = 0
   end do

   ! Process first molecule items with ULTRA-OPTIMIZED signature generation using direct 2D adjacency
   do i = 1, items1_count
      item_vertidx = assign_arrays%vertidcs1(items1_offset + i)

      ! ULTRA-OPTIMIZED: Generate compact signature using direct 2D adjacency access
      signature_length = 0
      do j = 1, assign_arrays%adj_counts1(item_vertidx)
         adj_atom = assign_arrays%adj_lists1(item_vertidx, j)
         part_ref_idx = assign_arrays%itemdir1_entries(read_link_idx, adj_atom)
         if (part_ref_idx /= 0) then
            signature_length = signature_length + 1
            signature_array(signature_length) = part_ref_idx
         end if
      end do

      ! Get relative index directly - no search needed
      target_relative_idx = find_child_part_array(assign_arrays, part_idx)

      ! Get absolute part index from relative index
      target_part_idx = assign_arrays%partree(part_idx)%child_indices(target_relative_idx)

      ! Add item to target child using relative index directly
      items1_trackers(target_relative_idx) = items1_trackers(target_relative_idx) + 1
      target_idx = assign_arrays%partree(target_part_idx)%items1_offset + items1_trackers(target_relative_idx)
      assign_arrays%vertidcs1(target_idx) = item_vertidx

      ! Update itemdir using 2D array - no offset calculation needed
      assign_arrays%itemdir1_entries(write_link_idx, item_vertidx) = target_part_idx
   end do

   ! Process second molecule items with ULTRA-OPTIMIZED signature generation using direct 2D adjacency
   do i = 1, items2_count
      item_vertidx = assign_arrays%vertidcs2(items2_offset + i)

      ! ULTRA-OPTIMIZED: Generate compact signature using direct 2D adjacency access
      signature_length = 0
      do j = 1, assign_arrays%adj_counts2(item_vertidx)
         adj_atom = assign_arrays%adj_lists2(item_vertidx, j)
         part_ref_idx = assign_arrays%itemdir2_entries(read_link_idx, adj_atom)
         if (part_ref_idx /= 0) then
            signature_length = signature_length + 1
            signature_array(signature_length) = part_ref_idx
         end if
      end do

      ! Get relative index directly - no search needed
      target_relative_idx = find_child_part_array(assign_arrays, part_idx)

      ! Get absolute part index from relative index
      target_part_idx = assign_arrays%partree(part_idx)%child_indices(target_relative_idx)

      ! Add item to target child using relative index directly
      items2_trackers(target_relative_idx) = items2_trackers(target_relative_idx) + 1
      target_idx = assign_arrays%partree(target_part_idx)%items2_offset + items2_trackers(target_relative_idx)
      assign_arrays%vertidcs2(target_idx) = item_vertidx

      ! Update itemdir using 2D array - no offset calculation needed
      assign_arrays%itemdir2_entries(write_link_idx, item_vertidx) = target_part_idx
   end do

   ! Collect assignment pairs from newly created leaf parts into subperm
   call collect_leaf_assignments(assign_arrays, part_idx, subperm)
end subroutine

subroutine recompute_mna_partition(assign_arrays, link_idx, subperm)
   type(array_trees_t), intent(inout) :: assign_arrays
   integer, intent(in) :: link_idx
   type(subperm_t), intent(inout) :: subperm
   integer :: next_link_idx, i, part_idx
   integer :: num_parts, partref_offset

   next_link_idx = link_idx + 1
   num_parts = assign_arrays%chain(link_idx)%num_parts
   partref_offset = assign_arrays%chain(link_idx)%partref_offset

   do i = 1, num_parts
      part_idx = assign_arrays%partref_entries(partref_offset + i)
      call resplit_part_mna(assign_arrays, part_idx, link_idx, next_link_idx, subperm)
   end do
end subroutine

subroutine distribute_part_items(assign_arrays, split_part_idx, child_branch_idx, &
      first_link_idx, chosen_item1_idx, chosen_item2_idx, subperm)
   ! Combined procedure: assignment + MNA recomputation
   ! Makes assignment (chosen_item1_idx-th item1 with chosen_item2_idx-th item2) then recomputes MNAs for the branch
   ! UPDATED: Now accepts both item1 and item2 indices to match original random assignment behavior
   type(array_trees_t), intent(inout) :: assign_arrays
   integer, intent(in) :: split_part_idx, child_branch_idx, first_link_idx, chosen_item1_idx, chosen_item2_idx
   type(subperm_t), intent(inout) :: subperm
   integer :: child_part1, child_part2, chosen_item1, chosen_item2, i, item_vertidx, target_idx
   integer :: items1_offset, items1_count, items2_offset, items2_count
   integer :: link_idx, num_links, link_offset

   ! === PART 1: PAIR ASSIGNMENT ===

   ! Extract commonly used offsets and values
   items1_offset = assign_arrays%partree(split_part_idx)%items1_offset
   items1_count = assign_arrays%partree(split_part_idx)%items1_count
   items2_offset = assign_arrays%partree(split_part_idx)%items2_offset
   items2_count = assign_arrays%partree(split_part_idx)%items2_count

   ! Get child parts using direct array access
   child_part1 = assign_arrays%partree(split_part_idx)%child_indices(1)
   child_part2 = assign_arrays%partree(split_part_idx)%child_indices(2)

   ! Get chosen items based on provided indices
   chosen_item1 = assign_arrays%vertidcs1(items1_offset + chosen_item1_idx)
   chosen_item2 = assign_arrays%vertidcs2(items2_offset + chosen_item2_idx)

   ! Assign chosen items to first child (direct placement)
   assign_arrays%vertidcs1(assign_arrays%partree(child_part1)%items1_offset + 1) = chosen_item1
   assign_arrays%vertidcs2(assign_arrays%partree(child_part1)%items2_offset + 1) = chosen_item2

   ! Update itemdir using 2D arrays - no offset calculation needed
   assign_arrays%itemdir1_entries(first_link_idx, chosen_item1) = child_part1
   assign_arrays%itemdir2_entries(first_link_idx, chosen_item2) = child_part1

   ! Copy remaining items1 to second child (skip the chosen item)
   target_idx = assign_arrays%partree(child_part2)%items1_offset
   do i = 1, items1_count
      if (i /= chosen_item1_idx) then
         item_vertidx = assign_arrays%vertidcs1(items1_offset + i)
         target_idx = target_idx + 1
         assign_arrays%vertidcs1(target_idx) = item_vertidx
         assign_arrays%itemdir1_entries(first_link_idx, item_vertidx) = child_part2
      end if
   end do

   ! Copy remaining items2 to second child (skip the chosen item)
   target_idx = assign_arrays%partree(child_part2)%items2_offset
   do i = 1, items2_count
      if (i /= chosen_item2_idx) then
         item_vertidx = assign_arrays%vertidcs2(items2_offset + i)
         target_idx = target_idx + 1
         assign_arrays%vertidcs2(target_idx) = item_vertidx
         assign_arrays%itemdir2_entries(first_link_idx, item_vertidx) = child_part2
      end if
   end do

   ! Collect assignment pairs from assignment (leaf parts created by the split)
   call collect_leaf_assignments(assign_arrays, split_part_idx, subperm)

   ! === PART 2: SCNA RECOMPUTATION ===

   num_links = assign_arrays%assigntree(child_branch_idx)%num_links
   link_offset = assign_arrays%assigntree(child_branch_idx)%link_offset

   do i = 1, num_links
      link_idx = link_offset + i
      call recompute_mna_partition(assign_arrays, link_idx, subperm)
   end do
end subroutine

recursive subroutine recurse_distribute_items_random(coords1, coords2, assign_arrays, branch_idx, best_perm)
   ! Random exploration - generates one assignment randomly using same traversal order as random module
   ! Similar to recurse_distribute_items_branch but picks one random assignment instead of exploring all
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: assign_arrays
   integer, intent(in) :: branch_idx
   type(subperm_t), intent(inout) :: best_perm

   integer :: rand_idx1, rand_idx2
   integer :: link_idx, branch_link_offset, branch_num_links
   integer :: child_branch_idx, first_link_idx, split_part_idx
   integer :: items1_count, items2_count, i

   ! Check if this is a leaf level (no more child branches)
   if (assign_arrays%assigntree(branch_idx)%num_children == 0) then
      return
   end if

   ! Process each child branch using same traversal order as random module
   do i = 1, assign_arrays%assigntree(branch_idx)%num_children
      child_branch_idx = assign_arrays%assigntree(branch_idx)%child_indices(i)
      first_link_idx = assign_arrays%assigntree(child_branch_idx)%link_offset + 1
      split_part_idx = assign_arrays%assigntree(child_branch_idx)%split_part_idx
      branch_link_offset = assign_arrays%assigntree(child_branch_idx)%link_offset
      branch_num_links = assign_arrays%assigntree(child_branch_idx)%num_links

      items1_count = assign_arrays%partree(split_part_idx)%items1_count
      items2_count = assign_arrays%partree(split_part_idx)%items2_count

      ! Generate random choices for both item1 and item2 indices (matching original random module)
      rand_idx1 = random_uniform_integer(1, items1_count)
      rand_idx2 = random_uniform_integer(1, items2_count)

      ! Make random assignment and recompute MNAs using existing procedure
      call distribute_part_items(assign_arrays, split_part_idx, child_branch_idx, first_link_idx, &
         rand_idx1, rand_idx2, best_perm)

      ! Recursively explore child branch
      call recurse_distribute_items_random(coords1, coords2, assign_arrays, child_branch_idx, best_perm)

      ! Reset state for next iteration - only reset links used by this branch
      do link_idx = branch_link_offset + 1, branch_link_offset + branch_num_links
         assign_arrays%itemdir1_entries(link_idx, :) = 0
         assign_arrays%itemdir2_entries(link_idx, :) = 0
      end do
   end do
end subroutine

subroutine distribute_items_random(coords1, coords2, assign_arrays, best_perm)
   ! Random exploration wrapper - generates one random assignment
   ! Similar to distribute_items_branch but generates random assignment instead of optimal
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: assign_arrays
   type(subperm_t), intent(out) :: best_perm
   ! Local variables
   real(rk) :: total_dist

!   call random_init(.true., .true.)

   ! Initialize random assignment
   call subperm_init(best_perm, assign_arrays%num_atoms1)

   ! Initialize assignment with preassigned pairs
   call collect_leaf_assignments(assign_arrays, 1, best_perm)

   ! Perform random exploration to generate one assignment (starting from root chain at index 1)
   call recurse_distribute_items_random(coords1, coords2, assign_arrays, 1, best_perm)

   ! Calculate distance from the random permutation array
   total_dist = total_sqdist(best_perm, coords1, coords2)
   write(stderr, '(A)') repeat("=", 60)
   write(stderr, '(A,I0,A,I0,A)') "Atoms assigned: ", best_perm%count, " out of ", &
         assign_arrays%num_atoms1, " total atoms"
   write(stderr, '(A,F10.4)') "Random assignment total squared distance: ", total_dist
   write(stderr, '(A)') repeat("=", 60)
   call check_subperm(best_perm)
end subroutine

subroutine collect_split_parts(assign_arrays, split_parts, num_split_parts)
   ! Collect all parts that need assignments (have >= 2 items in both molecules)
   type(array_trees_t), intent(in) :: assign_arrays
   integer, allocatable, intent(out) :: split_parts(:)
   integer, intent(out) :: num_split_parts
   integer :: i, count

   ! First pass: count split parts
   count = 0
   do i = 1, assign_arrays%total_chains
      if (assign_arrays%assigntree(i)%split_part_idx > 0) then
         if (assign_arrays%partree(assign_arrays%assigntree(i)%split_part_idx)%items1_count >= 2 .and. &
             assign_arrays%partree(assign_arrays%assigntree(i)%split_part_idx)%items2_count >= 2) then
            count = count + 1
         end if
      end if
   end do

   num_split_parts = count
   allocate(split_parts(num_split_parts))

   ! Second pass: collect split part indices
   count = 0
   do i = 1, assign_arrays%total_chains
      if (assign_arrays%assigntree(i)%split_part_idx > 0) then
         if (assign_arrays%partree(assign_arrays%assigntree(i)%split_part_idx)%items1_count >= 2 .and. &
             assign_arrays%partree(assign_arrays%assigntree(i)%split_part_idx)%items2_count >= 2) then
            count = count + 1
            split_parts(count) = assign_arrays%assigntree(i)%split_part_idx
         end if
      end if
   end do
end subroutine

recursive subroutine recurse_distribute_items_tree(coords1, coords2, assign_arrays, split_parts, &
                                num_split_parts, current_split_idx, this_perm, best_perm, min_dist)
   ! Recursively try all assignments for all split parts
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: assign_arrays
   integer, intent(in) :: split_parts(:), num_split_parts, current_split_idx
   type(subperm_t), intent(inout) :: this_perm, best_perm
   real(rk), intent(inout) :: min_dist

   integer :: split_part_idx, child_branch_idx, first_link_idx
   integer :: items2_count, i, j
   integer :: link_idx, branch_link_offset, branch_num_links
   type(subperm_t) :: saved_treeperm
   real(rk) :: total_dist

   ! Base case: we've made assignments for all split parts - evaluate complete permutation
   if (current_split_idx > num_split_parts) then
      combination_count = combination_count + 1
      total_dist = least_total_sqdist(this_perm, coords1, coords2)
      
      ! Update global best if this permutation is better
      if (total_dist < min_dist) then
         min_dist = total_dist
         best_perm = this_perm
      end if
      return
   end if

   ! Get the current split part
   split_part_idx = split_parts(current_split_idx)
   
   ! Find the child branch that uses this split part
   child_branch_idx = 0
   do i = 1, assign_arrays%total_chains
      if (assign_arrays%assigntree(i)%split_part_idx == split_part_idx) then
         child_branch_idx = i
         exit
      end if
   end do

   if (child_branch_idx == 0) then
      error stop 'Could not find child branch for split part'
   end if

   first_link_idx = assign_arrays%assigntree(child_branch_idx)%link_offset + 1
   branch_link_offset = assign_arrays%assigntree(child_branch_idx)%link_offset
   branch_num_links = assign_arrays%assigntree(child_branch_idx)%num_links

   items2_count = assign_arrays%partree(split_part_idx)%items2_count

   ! Try all possible assignments for this split part (only vary item2, keep item1 at index 1)
   do j = 1, items2_count
      ! Save current permutation state
      saved_treeperm = this_perm

      ! Make assignment for this split part (always use first item1, index=1)
      call distribute_part_items(assign_arrays, split_part_idx, child_branch_idx, first_link_idx, 1, j, &
         this_perm)

      ! Recursively try assignments for remaining split parts
      call recurse_distribute_items_tree(coords1, coords2, assign_arrays, split_parts, num_split_parts, &
         current_split_idx + 1, this_perm, best_perm, min_dist)

      ! Restore permutation state
      this_perm = saved_treeperm

      ! Reset itemdir state for this branch
      do link_idx = branch_link_offset + 1, branch_link_offset + branch_num_links
         assign_arrays%itemdir1_entries(link_idx, :) = 0
         assign_arrays%itemdir2_entries(link_idx, :) = 0
      end do
   end do
end subroutine

subroutine distribute_items_tree(coords1, coords2, assign_arrays, best_perm)
   ! DFS exploration of all permutations
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: assign_arrays
   type(subperm_t), intent(out) :: best_perm
   ! Local variables
   type(subperm_t) :: this_perm
   real(rk) :: min_dist
   integer, allocatable :: split_parts(:)
   integer :: num_split_parts, num_atoms

   num_atoms = assign_arrays%num_atoms1

   ! Initialize permutations
   call subperm_init(this_perm, num_atoms)
   call subperm_init(best_perm, num_atoms)

   ! Initialize both permutations with preassigned pairs
   call collect_leaf_assignments(assign_arrays, 1, this_perm)
   call collect_leaf_assignments(assign_arrays, 1, best_perm)

   ! Initialize global minimum distance
   min_dist = huge(1.0_rk)

   ! Initialize DFS exploration variables
   combination_count = 0

   ! Collect all split parts
   call collect_split_parts(assign_arrays, split_parts, num_split_parts)

   ! Try all assignments for all split parts
   call recurse_distribute_items_tree(coords1, coords2, assign_arrays, split_parts, num_split_parts, 1, &
      this_perm, best_perm, min_dist)

   ! Clean up
   deallocate(split_parts)

   ! Output results
!   write(stderr, '(A)') repeat("=", 60)
!   write(stderr, '(A,I0)') "Split parts found: ", num_split_parts
!   write(stderr, '(A,I0)') "Total permutations explored: ", combination_count
!   write(stderr, '(A,I0,A,I0,A)') "Atoms assigned: ", best_perm%count, " out of ", num_atoms, " total atoms"
!   write(stderr, '(A,F10.4)') "Minimum total squared distance: ", min_dist
!   write(stderr, '(A)') repeat("=", 60)
!   call check_subperm(best_perm)
end subroutine

recursive subroutine recurse_distribute_items_branch(coords1, coords2, assign_arrays, branch_idx, best_perm)
   ! DFS exploration of all assignment possibilities - finds permutation that minimizes total distance
   ! OPTIMIZED: Reduces allocations by reusing arrays within branch scope, but maintains isolation between branches
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)  ! coordinates needed for distance calculation
   type(array_trees_t), intent(inout) :: assign_arrays
   integer, intent(in) :: branch_idx
   type(subperm_t), intent(inout) :: best_perm

   integer :: child_branch_idx, first_link_idx, split_part_idx, i, items2_count, j
   integer :: link_idx, branch_link_offset, branch_num_links
   type(subperm_t) :: best_branch_perm, branch_perm
   real(rk) :: branch_dist, min_branch_dist
   integer :: num_atoms

   num_atoms = assign_arrays%num_atoms1

   ! Check if this is a leaf level (no more child branches)
   if (assign_arrays%assigntree(branch_idx)%num_children == 0) then
      combination_count = combination_count + 1
      return
   end if

   ! Process each child branch independently
   do i = 1, assign_arrays%assigntree(branch_idx)%num_children
      child_branch_idx = assign_arrays%assigntree(branch_idx)%child_indices(i)
      first_link_idx = assign_arrays%assigntree(child_branch_idx)%link_offset + 1
      split_part_idx = assign_arrays%assigntree(child_branch_idx)%split_part_idx
      branch_link_offset = assign_arrays%assigntree(child_branch_idx)%link_offset
      branch_num_links = assign_arrays%assigntree(child_branch_idx)%num_links

      items2_count = assign_arrays%partree(split_part_idx)%items2_count
      min_branch_dist = huge(1.0_rk)  ! Best distance for this specific branch

      call subperm_init(best_branch_perm, num_atoms)
      call subperm_init(branch_perm, num_atoms)

      ! Try pairing first item1 with each item2 to find best assignment for this branch
      do j = 1, items2_count
         ! Reset branch assignment for this iteration
         branch_perm%count = 0

         ! Make assignment and recompute MNAs in one combined operation (using first item1, index=1)
         call distribute_part_items(assign_arrays, split_part_idx, child_branch_idx, first_link_idx, 1, j, &
            branch_perm)

         ! Recursively explore subtree and collect child permutation
         call recurse_distribute_items_branch(coords1, coords2, assign_arrays, child_branch_idx, branch_perm)

         ! Calculate partial distance for this branch
         branch_dist = total_sqdist(branch_perm, coords1, coords2)

         ! Update best distance for this branch if this assignment is better
         if (branch_dist < min_branch_dist) then
            min_branch_dist = branch_dist
            best_branch_perm = branch_perm
         end if

         ! Reset state for next iteration - only reset links used by this branch
         do link_idx = branch_link_offset + 1, branch_link_offset + branch_num_links
            assign_arrays%itemdir1_entries(link_idx, :) = 0
            assign_arrays%itemdir2_entries(link_idx, :) = 0
         end do
      end do

      ! Update the optimal assignment with the best assigment from this branch
      call subperm_merge(best_perm, best_branch_perm)
   end do
end subroutine

subroutine distribute_items_branch(coords1, coords2, assign_arrays, best_perm)
   ! DFS exploration wrapper - finds optimal assignment among all possibilities
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: assign_arrays
   type(subperm_t), intent(out) :: best_perm
   ! Local variables
   integer :: num_atoms

   num_atoms = assign_arrays%num_atoms1

   ! Initialize optimal assignment
   call subperm_init(best_perm, num_atoms)

   ! Initialize optimal assignment with preassigned pairs
   call collect_leaf_assignments(assign_arrays, 1, best_perm)

   ! Initialize DFS exploration variables
   combination_count = 0

   ! Perform DFS exploration to find optimal assignment (starting from root chain at index 1)
   call recurse_distribute_items_branch(coords1, coords2, assign_arrays, 1, best_perm)

!block
!   integer :: assigned_count
!   real(rk) :: total_dist
!   assigned_count = best_perm%count
!   total_dist = total_sqdist(best_perm, coords1, coords2)
!   write(stderr, '(A)') repeat("=", 60)
!   write(stderr, '(A,I0)') "Assignment combinations probed: ", combination_count
!   write(stderr, '(A,I0,A,I0,A)') "Atoms assigned: ", assigned_count, " out of ", num_atoms, " total atoms"
!   write(stderr, '(A,F10.4)') "Minimum total squared distance: ", total_dist
!   write(stderr, '(A)') repeat("=", 60)
!   call check_subperm(best_perm)
!end block
end subroutine

end module
