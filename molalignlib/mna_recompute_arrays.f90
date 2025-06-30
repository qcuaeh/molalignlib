module mna_recompute_arrays
use parameters
use permutation
use spatial
use array_trees
implicit none

! Maximum possible number of children for a part
integer, parameter :: MAX_CHILDREN = 2

! Module-level signature workspace to eliminate allocations
integer :: signature_array(MAX_COORD)
integer :: signature_length

! DFS exploration variables
integer :: combination_count

! OPTIMIZATION: Reduced allocations by allocating once per branch instead of per combination

! NEW: Derived type to track assignments without using zero sentinels
type :: assignment_t
   integer, allocatable :: permutation(:)      ! permutation(i) = j means atom i -> atom j
   integer, allocatable :: assigned_indices(:) ! indices of assigned entries in permutation
   integer :: num_assigned                     ! number of assigned entries
end type

contains

subroutine init_assignment(assignment, num_atoms)
   type(assignment_t), intent(out) :: assignment
   integer, intent(in) :: num_atoms

   allocate(assignment%permutation(num_atoms))
   allocate(assignment%assigned_indices(num_atoms))
   assignment%num_assigned = 0
end subroutine

function signature_equivalence_array(array_trees, part_idx) result(equiv)
   ! OPTIMIZED: Fast path for length-1 signatures (most common case)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   logical :: equiv
   integer :: signature_frequencies, i, j

   if (signature_length /= array_trees%partree(part_idx)%signature_length) then
      equiv = .false.
      return
   end if

   ! FAST PATH: Direct comparison for length-1 signatures (most common)
   if (array_trees%partree(part_idx)%signature_length == 1) then
      equiv = (signature_array(1) == array_trees%partree(part_idx)%signature_values(1))
      return
   end if

   ! GENERIC PATH
   do i = 1, array_trees%partree(part_idx)%signature_unique_count
      signature_frequencies = 0

      ! Count matches in target signature
      do j = 1, signature_length
         if (signature_array(j) == array_trees%partree(part_idx)%signature_values(i)) then
            signature_frequencies = signature_frequencies + 1
         end if
      end do

      if (signature_frequencies /= array_trees%partree(part_idx)%signature_frequencies(i)) then
         equiv = .false.
         return
      end if
   end do

   equiv = .true.
end function

subroutine collect_leaf_assignments(array_trees, part_idx, assignment)
   ! Collect assignment pairs from leaf parts into assignment
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   type(assignment_t), intent(inout) :: assignment
   integer :: i, child_idx, item1_idx, item2_idx

   ! Check all children of this part
   do i = 1, array_trees%partree(part_idx)%num_children
      child_idx = array_trees%partree(part_idx)%child_indices(i)

      ! If this child is a leaf (no children), collect its assignment
      if (array_trees%partree(child_idx)%num_children == 0) then
         ! Verify this is a proper leaf part with exactly one item from each molecule
         if (array_trees%partree(child_idx)%items1_count == 1 .and. &
             array_trees%partree(child_idx)%items2_count == 1) then

            ! Get the assigned items
            item1_idx = array_trees%item1_values(array_trees%partree(child_idx)%items1_offset + 1)
            item2_idx = array_trees%item2_values(array_trees%partree(child_idx)%items2_offset + 1)

            ! Add item pair to assignment
            assignment%permutation(item1_idx) = item2_idx
            assignment%num_assigned = assignment%num_assigned + 1
            assignment%assigned_indices(assignment%num_assigned) = item1_idx
         end if
      end if
   end do
end subroutine

function find_child_part_array(array_trees, parent_idx) result(child_relative_idx)
   ! OPTIMIZED: Assumes exactly 2 children - if signature doesn't match first, it must match second
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: parent_idx
   integer :: child_relative_idx, first_child_idx, second_child_idx

   ! Get first child index directly
   first_child_idx = array_trees%partree(parent_idx)%child_indices(1)
   second_child_idx = array_trees%partree(parent_idx)%child_indices(2)

   ! Check if signature matches first or second child
   if (signature_equivalence_array(array_trees, first_child_idx)) then
      child_relative_idx = 1
   else if (signature_equivalence_array(array_trees, second_child_idx)) then
      child_relative_idx = 2
   else
      error stop 'Part signature does not match any child'
   end if
end function

subroutine resplit_part_mna(array_trees, part_idx, read_link_idx, write_link_idx, assignment)
! ULTRA-OPTIMIZED: Array-based version with direct 2D adjacency access for maximum performance
! UPDATED: Now collects assignment pairs from newly created leaf parts into assignment
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx, read_link_idx, write_link_idx
   type(assignment_t), intent(inout) :: assignment
   integer, dimension(MAX_CHILDREN) :: items1_trackers, items2_trackers
   integer :: i, j, target_relative_idx, target_part_idx, item_value, target_idx, part_ref_idx, adj_atom
   integer :: items1_offset, items1_count, items2_offset, items2_count
   integer :: num_children

   ! Extract commonly used values for readability
   items1_offset = array_trees%partree(part_idx)%items1_offset
   items1_count = array_trees%partree(part_idx)%items1_count
   items2_offset = array_trees%partree(part_idx)%items2_offset
   items2_count = array_trees%partree(part_idx)%items2_count
   num_children = array_trees%partree(part_idx)%num_children

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
      target_part_idx = array_trees%partree(part_idx)%child_indices(target_relative_idx)

      ! Add item to target child using relative index directly
      items1_trackers(target_relative_idx) = items1_trackers(target_relative_idx) + 1
      target_idx = array_trees%partree(target_part_idx)%items1_offset + items1_trackers(target_relative_idx)
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
      target_part_idx = array_trees%partree(part_idx)%child_indices(target_relative_idx)

      ! Add item to target child using relative index directly
      items2_trackers(target_relative_idx) = items2_trackers(target_relative_idx) + 1
      target_idx = array_trees%partree(target_part_idx)%items2_offset + items2_trackers(target_relative_idx)
      array_trees%item2_values(target_idx) = item_value

      ! Update itemdir using 2D array - no offset calculation needed
      array_trees%itemdir2_entries(write_link_idx, item_value) = target_part_idx
   end do

   ! Collect assignment pairs from newly created leaf parts into assignment
   call collect_leaf_assignments(array_trees, part_idx, assignment)
end subroutine

subroutine recompute_nextlevel_mnas(array_trees, link_idx, assignment)
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: link_idx
   type(assignment_t), intent(inout) :: assignment
   integer :: next_link_idx, i, part_idx
   integer :: num_parts, partref_offset

   next_link_idx = link_idx + 1
   num_parts = array_trees%chain(link_idx)%num_parts
   partref_offset = array_trees%chain(link_idx)%partref_offset

   do i = 1, num_parts
      part_idx = array_trees%partref_entries(partref_offset + i)
      call resplit_part_mna(array_trees, part_idx, link_idx, next_link_idx, assignment)
   end do
end subroutine

subroutine assign_and_recompute_mnas(array_trees, split_part_idx, child_branch_idx, &
      first_link_idx, chosen_item1_idx, chosen_item2_idx, assignment)
   ! Combined procedure: assignment + MNA recomputation
   ! Makes assignment (chosen_item1_idx-th item1 with chosen_item2_idx-th item2) then recomputes MNAs for the branch
   ! UPDATED: Now accepts both item1 and item2 indices to match original random assignment behavior
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: split_part_idx, child_branch_idx, first_link_idx, chosen_item1_idx, chosen_item2_idx
   type(assignment_t), intent(inout) :: assignment
   integer :: child_part1, child_part2, chosen_item1, chosen_item2, i, item_value, target_idx
   integer :: items1_offset, items1_count, items2_offset, items2_count
   integer :: link_idx, num_links, link_offset

   ! === PART 1: ASSIGNMENT ===

   ! Extract commonly used offsets and values
   items1_offset = array_trees%partree(split_part_idx)%items1_offset
   items1_count = array_trees%partree(split_part_idx)%items1_count
   items2_offset = array_trees%partree(split_part_idx)%items2_offset
   items2_count = array_trees%partree(split_part_idx)%items2_count

   ! Get child parts using direct array access
   child_part1 = array_trees%partree(split_part_idx)%child_indices(1)
   child_part2 = array_trees%partree(split_part_idx)%child_indices(2)

   ! Get chosen items based on provided indices
   chosen_item1 = array_trees%item1_values(items1_offset + chosen_item1_idx)
   chosen_item2 = array_trees%item2_values(items2_offset + chosen_item2_idx)

   ! Assign chosen items to first child (direct placement)
   array_trees%item1_values(array_trees%partree(child_part1)%items1_offset + 1) = chosen_item1
   array_trees%item2_values(array_trees%partree(child_part1)%items2_offset + 1) = chosen_item2

   ! Update itemdir using 2D arrays - no offset calculation needed
   array_trees%itemdir1_entries(first_link_idx, chosen_item1) = child_part1
   array_trees%itemdir2_entries(first_link_idx, chosen_item2) = child_part1

   ! Copy remaining items1 to second child (skip the chosen item)
   target_idx = array_trees%partree(child_part2)%items1_offset
   do i = 1, items1_count
      if (i /= chosen_item1_idx) then
         item_value = array_trees%item1_values(items1_offset + i)
         target_idx = target_idx + 1
         array_trees%item1_values(target_idx) = item_value
         array_trees%itemdir1_entries(first_link_idx, item_value) = child_part2
      end if
   end do

   ! Copy remaining items2 to second child (skip the chosen item)
   target_idx = array_trees%partree(child_part2)%items2_offset
   do i = 1, items2_count
      if (i /= chosen_item2_idx) then
         item_value = array_trees%item2_values(items2_offset + i)
         target_idx = target_idx + 1
         array_trees%item2_values(target_idx) = item_value
         array_trees%itemdir2_entries(first_link_idx, item_value) = child_part2
      end if
   end do

   ! Collect assignment pairs from assignment (leaf parts created by the split)
   call collect_leaf_assignments(array_trees, split_part_idx, assignment)

   ! === PART 2: MNA RECOMPUTATION ===

   num_links = array_trees%assigntree(child_branch_idx)%num_links
   link_offset = array_trees%assigntree(child_branch_idx)%link_offset

   do i = 1, num_links
      link_idx = link_offset + i
      call recompute_nextlevel_mnas(array_trees, link_idx, assignment)
   end do
end subroutine

recursive subroutine redistribute_items_dfs_recursive(coords1, coords2, array_trees, branch_idx, optimal_assignment)
   ! DFS exploration of all assignment possibilities - finds permutation that minimizes total distance
   ! OPTIMIZED: Reduces allocations by reusing arrays within branch scope, but maintains isolation between branches
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)  ! coordinates needed for distance calculation
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: branch_idx
   type(assignment_t), intent(inout) :: optimal_assignment

   integer :: child_branch_idx, first_link_idx, split_part_idx, i, items2_count, j
   integer :: link_idx, branch_link_offset, branch_num_links
   type(assignment_t) :: best_branch_assignment, branch_assignment
   real(rk) :: branch_distance, best_branch_distance
   integer :: num_atoms

   num_atoms = array_trees%num_atoms1

   ! Check if this is a leaf level (no more child branches)
   if (array_trees%assigntree(branch_idx)%num_children == 0) then
      combination_count = combination_count + 1
      return
   end if

   ! Process each child branch independently
   do i = 1, array_trees%assigntree(branch_idx)%num_children
      child_branch_idx = array_trees%assigntree(branch_idx)%child_indices(i)
      first_link_idx = array_trees%assigntree(child_branch_idx)%link_offset + 1
      split_part_idx = array_trees%assigntree(child_branch_idx)%split_part_idx
      branch_link_offset = array_trees%assigntree(child_branch_idx)%link_offset
      branch_num_links = array_trees%assigntree(child_branch_idx)%num_links

      items2_count = array_trees%partree(split_part_idx)%items2_count
      best_branch_distance = huge(1.0_rk)  ! Best distance for this specific branch

      call init_assignment(best_branch_assignment, num_atoms)
      call init_assignment(branch_assignment, num_atoms)

      ! Try pairing first item1 with each item2 to find best assignment for this branch
      do j = 1, items2_count
         ! Reset branch assignment for this iteration
         branch_assignment%num_assigned = 0

         ! Make assignment and recompute MNAs in one combined operation (using first item1, index=1)
         call assign_and_recompute_mnas(array_trees, split_part_idx, child_branch_idx, first_link_idx, 1, j, &
            branch_assignment)

         ! Recursively explore subtree and collect child permutation
         call redistribute_items_dfs_recursive(coords1, coords2, array_trees, child_branch_idx, branch_assignment)

         ! Calculate partial distance for this branch
         branch_distance = totsqdist(branch_assignment%assigned_indices(1:branch_assignment%num_assigned), &
            branch_assignment%permutation, coords1, coords2)

         ! Update best distance for this branch if this assignment is better
         if (branch_distance < best_branch_distance) then
            best_branch_distance = branch_distance
            best_branch_assignment = branch_assignment
         end if

         ! Reset state for next iteration - only reset links used by this branch
         do link_idx = branch_link_offset + 1, branch_link_offset + branch_num_links
            array_trees%itemdir1_entries(link_idx, :) = 0
            array_trees%itemdir2_entries(link_idx, :) = 0
         end do
      end do

      ! Update the optimal assignment with the best assigment from this branch
      call update_assignment(optimal_assignment, best_branch_assignment)
   end do
end subroutine

recursive subroutine redistribute_items_random_recursive(coords1, coords2, array_trees, branch_idx, assignment)
   ! Random exploration - generates one assignment randomly using same traversal order as random module
   ! Similar to redistribute_items_dfs_recursive but picks one random assignment instead of exploring all
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: branch_idx
   type(assignment_t), intent(inout) :: assignment

   integer :: child_branch_idx, first_link_idx, split_part_idx, i, items1_count, items2_count
   integer :: random_choice1, random_choice2
   integer :: link_idx, branch_link_offset, branch_num_links
   real :: random_real

   ! Check if this is a leaf level (no more child branches)
   if (array_trees%assigntree(branch_idx)%num_children == 0) then
      return
   end if

   ! Process each child branch using same traversal order as random module
   do i = 1, array_trees%assigntree(branch_idx)%num_children
      child_branch_idx = array_trees%assigntree(branch_idx)%child_indices(i)
      first_link_idx = array_trees%assigntree(child_branch_idx)%link_offset + 1
      split_part_idx = array_trees%assigntree(child_branch_idx)%split_part_idx
      branch_link_offset = array_trees%assigntree(child_branch_idx)%link_offset
      branch_num_links = array_trees%assigntree(child_branch_idx)%num_links

      items1_count = array_trees%partree(split_part_idx)%items1_count
      items2_count = array_trees%partree(split_part_idx)%items2_count

      ! Generate random choices for both item1 and item2 indices (matching original random module)
      call random_number(random_real)
      random_choice1 = int(random_real * items1_count) + 1
      call random_number(random_real)
      random_choice2 = int(random_real * items2_count) + 1

      ! Make random assignment and recompute MNAs using existing procedure
      call assign_and_recompute_mnas(array_trees, split_part_idx, child_branch_idx, first_link_idx, &
         random_choice1, random_choice2, assignment)

      ! Recursively explore child branch
      call redistribute_items_random_recursive(coords1, coords2, array_trees, child_branch_idx, assignment)

      ! Reset state for next iteration - only reset links used by this branch
      do link_idx = branch_link_offset + 1, branch_link_offset + branch_num_links
         array_trees%itemdir1_entries(link_idx, :) = 0
         array_trees%itemdir2_entries(link_idx, :) = 0
      end do
   end do
end subroutine

subroutine update_assignment(target, source)
   ! Merge source assignment into target assignment
   type(assignment_t), intent(inout) :: target
   type(assignment_t), intent(in) :: source
   integer :: i, atom1_idx, atom2_idx

   do i = 1, source%num_assigned
      atom1_idx = source%assigned_indices(i)
      atom2_idx = source%permutation(atom1_idx)

      ! Check for conflicts in existing assignments
      if (any(target%assigned_indices(1:target%num_assigned) == atom1_idx)) then
         write(stderr, '(A,I0,A)') &
            "ERROR: Attempting to overwrite assignment at position ", atom1_idx
         error stop "Assignment update conflict"
      end if

      ! Add assignment directly
      target%num_assigned = target%num_assigned + 1
      target%assigned_indices(target%num_assigned) = atom1_idx
      target%permutation(atom1_idx) = atom2_idx
   end do
end subroutine

subroutine redistribute_items_dfs(coords1, coords2, array_trees, optimal_permutation, total_distance)
   ! DFS exploration wrapper - finds optimal assignment among all possibilities
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: array_trees
   integer, allocatable, intent(out) :: optimal_permutation(:)
   real(rk), intent(out) :: total_distance

   type(assignment_t) :: optimal_assignment
   integer :: num_atoms, assigned_count

   num_atoms = array_trees%num_atoms1

   ! Initialize optimal assignment
   call init_assignment(optimal_assignment, num_atoms)

   ! Initialize optimal assignment with preassigned pairs
   call collect_leaf_assignments(array_trees, 1, optimal_assignment)

   ! Initialize DFS exploration variables
   combination_count = 0

   write(stderr, '(A)') "=== Starting DFS exploration of all assignment possibilities ==="

   ! Perform DFS exploration to find optimal assignment (starting from root chain at index 1)
   call redistribute_items_dfs_recursive(coords1, coords2, array_trees, 1, optimal_assignment)

   ! Copy final permutation array from assignment
   optimal_permutation = optimal_assignment%permutation

   ! Calculate distance from the optimal permutation array for verification
   total_distance = totsqdist(optimal_assignment%assigned_indices(1:optimal_assignment%num_assigned), &
      optimal_assignment%permutation, coords1, coords2)

   ! Count assigned atoms
   assigned_count = optimal_assignment%num_assigned

   ! Report final results
   write(stderr, '(A)') repeat("=", 60)
   write(stderr, '(A,I0)') "Assignment combinations probed: ", combination_count
   write(stderr, '(A,I0,A,I0,A)') "Atoms assigned: ", assigned_count, " out of ", num_atoms, " total atoms"
   write(stderr, '(A,F10.4)') "Optimal total squared distance: ", total_distance
   write(stderr, '(A)') repeat("=", 60)

   ! Validate permutation consistency
   call validate_perm(optimal_permutation)
end subroutine

subroutine redistribute_items_random(coords1, coords2, array_trees, random_permutation, total_distance)
   ! Random exploration wrapper - generates one random assignment
   ! Similar to redistribute_items_dfs but generates random assignment instead of optimal
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: array_trees
   integer, allocatable, intent(out) :: random_permutation(:)
   real(rk), intent(out) :: total_distance

   type(assignment_t) :: random_assignment
   integer :: num_atoms, assigned_count

   num_atoms = array_trees%num_atoms1

   ! Initialize random assignment
   call init_assignment(random_assignment, num_atoms)

   ! Initialize assignment with preassigned pairs
   call collect_leaf_assignments(array_trees, 1, random_assignment)

   write(stderr, '(A)') "=== Starting random assignment generation ==="

   ! Perform random exploration to generate one assignment (starting from root chain at index 1)
   call redistribute_items_random_recursive(coords1, coords2, array_trees, 1, random_assignment)

   ! Copy final permutation array from assignment
   random_permutation = random_assignment%permutation

   ! Calculate distance from the random permutation array
   total_distance = totsqdist(random_assignment%assigned_indices(1:random_assignment%num_assigned), &
      random_assignment%permutation, coords1, coords2)

   ! Count assigned atoms
   assigned_count = random_assignment%num_assigned

   ! Report final results
   write(stderr, '(A)') repeat("=", 60)
   write(stderr, '(A,I0,A,I0,A)') "Atoms assigned: ", assigned_count, " out of ", num_atoms, " total atoms"
   write(stderr, '(A,F10.4)') "Random assignment total squared distance: ", total_distance
   write(stderr, '(A)') repeat("=", 60)

   ! Validate permutation consistency
   call validate_perm(random_permutation)
end subroutine

end module
