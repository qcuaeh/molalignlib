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

module assignment_conformer
use parameters
use random
use common_types
use permutation
use euclidean
use indexed_list_types
implicit none
private
public assign_atoms_greedy
public assign_atoms_global
public assign_atoms_local
public assign_atoms_local_pruned

! Maximum number of part children
integer(ik), parameter :: MAX_CHILD = 10

! DFS exploration variables
integer(ik) :: n_combinations

contains

function signature_equivalence_array(cache_arrays, part_idx, signature_size, signature) result(equiv)
   ! OPTIMIZED: Fast path for length-1 signatures (most common case)
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: part_idx
   integer(ik), intent(in) :: signature_size, signature(:)
   logical(lk) :: equiv
   integer(ik) :: signature_frequencies, i, j

   if (signature_size /= cache_arrays%partree(part_idx)%signature_size) then
      equiv = .FALSE.
      return
   end if

   ! FAST PATH: Direct comparison for length-1 signatures (most common)
   if (cache_arrays%partree(part_idx)%signature_size == 1) then
      equiv = (signature(1) == cache_arrays%partree(part_idx)%signature_values(1))
      return
   end if

   ! GENERIC PATH
   do i = 1, cache_arrays%partree(part_idx)%signature_unique_count
      signature_frequencies = 0

      ! Count matches in target signature
      do j = 1, signature_size
         if (signature(j) == cache_arrays%partree(part_idx)%signature_values(i)) then
            signature_frequencies = signature_frequencies + 1
         end if
      end do

      if (signature_frequencies /= cache_arrays%partree(part_idx)%signature_frequencies(i)) then
         equiv = .FALSE.
         return
      end if
   end do

   equiv = .TRUE.
end function

function find_child_part_array(cache_arrays, parent_idx, signature_size, signature) result(child_relative_idx)
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: parent_idx
   integer(ik), intent(in) :: signature_size, signature(:)
   integer(ik) :: child_relative_idx
   integer(ik) :: n_children, child_idx, i

   ! Get first child index directly
   n_children = cache_arrays%partree(parent_idx)%n_children

   ! Find the child that matches the signature
   do i = 1, n_children
      child_idx = cache_arrays%partree(parent_idx)%child_indices(i)
      if (signature_equivalence_array(cache_arrays, child_idx, signature_size, signature)) then
         child_relative_idx = i
         return
      end if
   end do

   ! If the signature did not match any child it is a bug
   error stop 'Part signature did not match any child'
end function

subroutine get_mapping(atom_map, mapping1)
! Export a completed atom_map as a permutation array. The assignment tree
! only contains included atoms, so once the search is over every atom has
! been assigned. Only defined entries are read.
   type(partmap_t), intent(in) :: atom_map
   integer(ik), dimension(:), allocatable, intent(out) :: mapping1
   integer(ik) :: i, i1

   allocate (mapping1(size(atom_map%mapping)))
   mapping1 = 0
   do i = 1, atom_map%subset_size
      i1 = atom_map%subset(i)
      mapping1(i1) = atom_map%mapping(i1)
   end do

   if (DEBUG_TESTS) then
      if (.not. is_permutation(mapping1)) then
         error stop 'get_mapping: assignment is incomplete'
      end if
   end if
end subroutine

subroutine collect_leaf_assignments(cache_arrays, part_idx, atom_map)
   ! Collect assignment pairs from leaf parts into assignment
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: part_idx
   type(partmap_t), intent(inout) :: atom_map
   integer(ik) :: i, child_idx, item1_idx, item2_idx

   ! Check all children of this part
   do i = 1, cache_arrays%partree(part_idx)%n_children
      child_idx = cache_arrays%partree(part_idx)%child_indices(i)

      ! If this child is a leaf (no children), collect its assignment
      if (cache_arrays%partree(child_idx)%n_children == 0) then
         ! Verify this is a proper leaf part with exactly one item from each molecule
         if (cache_arrays%partree(child_idx)%n_items1 == 1 .and. &
             cache_arrays%partree(child_idx)%n_items2 == 1) then

            ! Get the assigned items
            item1_idx = cache_arrays%atomidcs1(cache_arrays%partree(child_idx)%items1_offset + 1)
            item2_idx = cache_arrays%atomidcs2(cache_arrays%partree(child_idx)%items2_offset + 1)

            ! Add item pair to assignment
            call submap_add(atom_map, item1_idx, item2_idx)
         end if
      end if
   end do
end subroutine

subroutine update_hna_part(cache_arrays, part_idx, read_link_idx, write_link_idx, atom_map)
! ULTRA-OPTIMIZED: Array-based version with direct 2D adjacency access for maximum performance
! UPDATED: Now collects assignment pairs from newly created leaf parts into atom_map
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(in) :: part_idx, read_link_idx, write_link_idx
   type(partmap_t), intent(inout) :: atom_map
   integer(ik), save :: signature_size, signature(MAX_COORDNUM)
   integer(ik), dimension(MAX_CHILD), save :: items1_trackers, items2_trackers
   integer(ik) :: target_relative_idx, target_part_idx, item_idx, target_idx, part_ref_idx
   integer(ik) :: items1_offset, n_items1, items2_offset, n_items2, n_children
   integer(ik) :: adj_atom
   integer(ik) :: i, j

   ! Extract commonly used values for readability
   items1_offset = cache_arrays%partree(part_idx)%items1_offset
   n_items1 = cache_arrays%partree(part_idx)%n_items1
   items2_offset = cache_arrays%partree(part_idx)%items2_offset
   n_items2 = cache_arrays%partree(part_idx)%n_items2
   n_children = cache_arrays%partree(part_idx)%n_children

   ! INITIALIZATION: Reset trackers using relative indices (1 to n_children)
   do i = 1, n_children
      items1_trackers(i) = 0
      items2_trackers(i) = 0
   end do

   ! Process first molecule items with ULTRA-OPTIMIZED signature generation using direct 2D adjacency
   do i = 1, n_items1
      item_idx = cache_arrays%atomidcs1(items1_offset + i)

      ! ULTRA-OPTIMIZED: Generate compact signature using direct 2D adjacency access
      signature_size = 0
      do j = 1, cache_arrays%adjcs1_cn(item_idx)
         adj_atom = cache_arrays%adjcs1_list(item_idx, j)
         part_ref_idx = cache_arrays%itemdir1_entries(adj_atom, read_link_idx)
         if (part_ref_idx /= 0) then
            signature_size = signature_size + 1
            ! Inlined edge_code(part_ref_idx, bond type) from module adjacency
            signature(signature_size) = part_ref_idx*BOND_TYPE_RADIX &
                                      + cache_arrays%adjcs1_bondtype(item_idx, j)
         end if
      end do

      ! Get relative index directly - no search needed
      target_relative_idx = find_child_part_array(cache_arrays, part_idx, signature_size, signature)

      ! Get absolute part index from relative index
      target_part_idx = cache_arrays%partree(part_idx)%child_indices(target_relative_idx)

      ! Add item to target child using relative index directly
      items1_trackers(target_relative_idx) = items1_trackers(target_relative_idx) + 1
      target_idx = cache_arrays%partree(target_part_idx)%items1_offset &
                 + items1_trackers(target_relative_idx)
      cache_arrays%atomidcs1(target_idx) = item_idx

      ! Update itemdir using 2D array - no offset calculation needed
      cache_arrays%itemdir1_entries(item_idx, write_link_idx) = target_part_idx
   end do

   ! Process second molecule items with ULTRA-OPTIMIZED signature generation using direct 2D adjacency
   do i = 1, n_items2
      item_idx = cache_arrays%atomidcs2(items2_offset + i)

      ! ULTRA-OPTIMIZED: Generate compact signature using direct 2D adjacency access
      signature_size = 0
      do j = 1, cache_arrays%adjcs2_cn(item_idx)
         adj_atom = cache_arrays%adjcs2_list(item_idx, j)
         part_ref_idx = cache_arrays%itemdir2_entries(adj_atom, read_link_idx)
         if (part_ref_idx /= 0) then
            signature_size = signature_size + 1
            ! Inlined edge_code(part_ref_idx, bond type) from module adjacency
            signature(signature_size) = part_ref_idx*BOND_TYPE_RADIX &
                                      + cache_arrays%adjcs2_bondtype(item_idx, j)
         end if
      end do

      ! Get relative index directly - no search needed
      target_relative_idx = find_child_part_array(cache_arrays, part_idx, signature_size, signature)

      ! Get absolute part index from relative index
      target_part_idx = cache_arrays%partree(part_idx)%child_indices(target_relative_idx)

      ! Add item to target child using relative index directly
      items2_trackers(target_relative_idx) = items2_trackers(target_relative_idx) + 1
      target_idx = cache_arrays%partree(target_part_idx)%items2_offset &
                 + items2_trackers(target_relative_idx)
      cache_arrays%atomidcs2(target_idx) = item_idx

      ! Update itemdir using 2D array - no offset calculation needed
      cache_arrays%itemdir2_entries(item_idx, write_link_idx) = target_part_idx
   end do

   ! Collect assignment pairs from newly created leaf parts into atom_map
   call collect_leaf_assignments(cache_arrays, part_idx, atom_map)
end subroutine

subroutine assign_pair_to_children(cache_arrays, split_part_idx, first_link_idx, &
      chosen_item1_idx, chosen_item2_idx, atom_map)
   ! Assigns chosen items to first child and remaining items to second child
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(in) :: split_part_idx, first_link_idx, chosen_item1_idx, chosen_item2_idx
   type(partmap_t), intent(inout) :: atom_map
   integer(ik) :: child_part1, child_part2, chosen_item1, chosen_item2, i, item_idx, target_idx
   integer(ik) :: items1_offset, n_items1, items2_offset, n_items2

   ! Extract commonly used offsets and values
   items1_offset = cache_arrays%partree(split_part_idx)%items1_offset
   n_items1 = cache_arrays%partree(split_part_idx)%n_items1
   items2_offset = cache_arrays%partree(split_part_idx)%items2_offset
   n_items2 = cache_arrays%partree(split_part_idx)%n_items2

   ! Get child parts using direct array access
   child_part1 = cache_arrays%partree(split_part_idx)%child_indices(1)
   child_part2 = cache_arrays%partree(split_part_idx)%child_indices(2)

   ! Get chosen items based on provided indices
   chosen_item1 = cache_arrays%atomidcs1(items1_offset + chosen_item1_idx)
   chosen_item2 = cache_arrays%atomidcs2(items2_offset + chosen_item2_idx)

   ! Assign chosen items to first child (direct placement)
   cache_arrays%atomidcs1(cache_arrays%partree(child_part1)%items1_offset + 1) = chosen_item1
   cache_arrays%atomidcs2(cache_arrays%partree(child_part1)%items2_offset + 1) = chosen_item2

   ! Update itemdir using 2D arrays - no offset calculation needed
   cache_arrays%itemdir1_entries(chosen_item1, first_link_idx) = child_part1
   cache_arrays%itemdir2_entries(chosen_item2, first_link_idx) = child_part1

   ! Copy remaining items1 to second child (skip the chosen item)
   target_idx = cache_arrays%partree(child_part2)%items1_offset
   do i = 1, n_items1
      if (i /= chosen_item1_idx) then
         item_idx = cache_arrays%atomidcs1(items1_offset + i)
         target_idx = target_idx + 1
         cache_arrays%atomidcs1(target_idx) = item_idx
         cache_arrays%itemdir1_entries(item_idx, first_link_idx) = child_part2
      end if
   end do

   ! Copy remaining items2 to second child (skip the chosen item)
   target_idx = cache_arrays%partree(child_part2)%items2_offset
   do i = 1, n_items2
      if (i /= chosen_item2_idx) then
         item_idx = cache_arrays%atomidcs2(items2_offset + i)
         target_idx = target_idx + 1
         cache_arrays%atomidcs2(target_idx) = item_idx
         cache_arrays%itemdir2_entries(item_idx, first_link_idx) = child_part2
      end if
   end do

   ! Collect assignment pairs from assignment (leaf parts created by the split)
   call collect_leaf_assignments(cache_arrays, split_part_idx, atom_map)
end subroutine

subroutine assign_branch_atoms(cache_arrays, split_part_idx, child_branch_idx, &
      first_link_idx, chosen_item1_idx, chosen_item2_idx, atom_map)
   ! Combined procedure: assignment + HNA recomputation
   ! Makes assignment (chosen_item1_idx-th item1 with chosen_item2_idx-th item2)
   ! then recomputes HNAs for the branch
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(in) :: split_part_idx, child_branch_idx, first_link_idx
   integer(ik), intent(in) :: chosen_item1_idx, chosen_item2_idx
   type(partmap_t), intent(inout) :: atom_map
   integer(ik) :: link_idx, next_link_idx, part_idx
   integer(ik) :: n_links, n_parts, link_offset, partref_offset
   integer(ik) :: i, j

   ! === Pair assignment ===
   call assign_pair_to_children(cache_arrays, split_part_idx, first_link_idx, &
         chosen_item1_idx, chosen_item2_idx, atom_map)

   ! === HNA chain recomputation ===
   n_links = cache_arrays%assigntree(child_branch_idx)%n_links
   link_offset = cache_arrays%assigntree(child_branch_idx)%link_offset

   if (DEBUG_TESTS) then
      ! The loop below reads link i and writes link i+1, but the searches only
      ! reset links link_offset+1 .. link_offset+n_links. The last link must
      ! therefore have no parts, or its writes would land one link past the
      ! branch and survive the reset.
      if (n_links > 0) then
         if (cache_arrays%chain(link_offset + n_links)%n_parts /= 0) then
            error stop 'Last link of the branch has parts to update'
         end if
      end if
   end if

   do i = 1, n_links
      link_idx = link_offset + i
      next_link_idx = link_idx + 1
      n_parts = cache_arrays%chain(link_idx)%n_parts
      partref_offset = cache_arrays%chain(link_idx)%partref_offset

      do j = 1, n_parts
         part_idx = cache_arrays%partref_entries(partref_offset + j)
         call update_hna_part(cache_arrays, part_idx, link_idx, next_link_idx, atom_map)
      end do
   end do
end subroutine

subroutine collect_split_parts(cache_arrays, split_parts, n_split_parts)
   ! Collect all parts that need assignments (have > 1 items in both molecules)
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), allocatable, intent(out) :: split_parts(:)
   integer(ik), intent(out) :: n_split_parts
   integer(ik) :: i, count

   ! First pass: count split parts
   count = 0
   do i = 1, cache_arrays%total_chains
      if (cache_arrays%assigntree(i)%split_part_idx > 0) then
         if (cache_arrays%partree(cache_arrays%assigntree(i)%split_part_idx)%n_items1 > 1 .and. &
             cache_arrays%partree(cache_arrays%assigntree(i)%split_part_idx)%n_items2 > 1) then
            count = count + 1
         end if
      end if
   end do

   n_split_parts = count
   allocate(split_parts(n_split_parts))

   ! Second pass: collect split part indices
   count = 0
   do i = 1, cache_arrays%total_chains
      if (cache_arrays%assigntree(i)%split_part_idx > 0) then
         if (cache_arrays%partree(cache_arrays%assigntree(i)%split_part_idx)%n_items1 > 1 .and. &
             cache_arrays%partree(cache_arrays%assigntree(i)%split_part_idx)%n_items2 > 1) then
            count = count + 1
            split_parts(count) = cache_arrays%assigntree(i)%split_part_idx
         end if
      end if
   end do
end subroutine

recursive subroutine recur_assign_atoms_greedy(coords1, coords2, cache_arrays, &
      branch_idx, greedy_map)
   ! Greedy exploration - always picks the closest pair at each split
   ! Similar to recur_assign_atoms_random but chooses minimum distance instead of random
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(in) :: branch_idx
   type(partmap_t), intent(inout) :: greedy_map

   integer(ik) :: i, idx1, idx2
   integer(ik) :: branch_link_offset, branch_num_links
   integer(ik) :: child_branch_idx, first_link_idx, split_part_idx
   integer(ik) :: n_items1, n_items2
   real(rk) :: min_dist, current_dist
   integer(ik) :: greedy_idx1, greedy_idx2
   integer(ik) :: item1_idx, item2_idx
   integer(ik) :: items1_offset, items2_offset

   ! Check if this is a leaf level (no more child branches)
   if (cache_arrays%assigntree(branch_idx)%n_children == 0) then
      return
   end if

   ! Process each child branch using same traversal order as random module
   do i = 1, cache_arrays%assigntree(branch_idx)%n_children
      child_branch_idx = cache_arrays%assigntree(branch_idx)%child_indices(i)
      first_link_idx = cache_arrays%assigntree(child_branch_idx)%link_offset + 1
      split_part_idx = cache_arrays%assigntree(child_branch_idx)%split_part_idx
      branch_link_offset = cache_arrays%assigntree(child_branch_idx)%link_offset
      branch_num_links = cache_arrays%assigntree(child_branch_idx)%n_links

      n_items1 = cache_arrays%partree(split_part_idx)%n_items1
      n_items2 = cache_arrays%partree(split_part_idx)%n_items2
      items1_offset = cache_arrays%partree(split_part_idx)%items1_offset
      items2_offset = cache_arrays%partree(split_part_idx)%items2_offset

      ! Find the closest pair (greedy choice)
      min_dist = huge(min_dist)
      greedy_idx1 = 1
      greedy_idx2 = 1

      do idx1 = 1, n_items1
         item1_idx = cache_arrays%atomidcs1(items1_offset + idx1)
         do idx2 = 1, n_items2
            item2_idx = cache_arrays%atomidcs2(items2_offset + idx2)

            ! Calculate squared distance between atoms
            current_dist = sum((coords1(:, item1_idx) - coords2(:, item2_idx))**2)

            if (current_dist < min_dist) then
               min_dist = current_dist
               greedy_idx1 = idx1
               greedy_idx2 = idx2
            end if
         end do
      end do

      ! Make greedy assignment (closest pair)
      call assign_branch_atoms(cache_arrays, split_part_idx, child_branch_idx, &
            first_link_idx, greedy_idx1, greedy_idx2, greedy_map)

      ! Recursively explore child branch
      call recur_assign_atoms_greedy(coords1, coords2, cache_arrays, child_branch_idx, &
            greedy_map)

      ! Reset state for next iteration - only reset links used by this branch
      cache_arrays%itemdir1_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
      cache_arrays%itemdir2_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
   end do
end subroutine

subroutine assign_atoms_greedy(coords1, coords2, cache_arrays, mapping1, mapdist)
   ! Greedy exploration wrapper with recursive JVC optimization
   ! At each split level, applies JVC to find optimal pairing for that subset
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), dimension(:), allocatable, intent(out) :: mapping1
   real(rk), intent(out) :: mapdist
   ! Local variables
   type(partmap_t) :: greedy_map

   ! Initialize greedy assignment
   call submap_init(greedy_map, cache_arrays%n_atoms1)

   ! Initialize assignment with preassigned pairs
   call collect_leaf_assignments(cache_arrays, 1, greedy_map)

   ! Perform greedy exploration with JVC optimization (starting from root chain at index 1)
   call recur_assign_atoms_greedy(coords1, coords2, cache_arrays, 1, greedy_map)

   ! Calculate total distance
   mapdist = sqdistsum(greedy_map%subset(:greedy_map%subset_size), greedy_map%mapping, coords1, coords2)

   call get_mapping(greedy_map, mapping1)
end subroutine

recursive subroutine recur_assign_atoms_global(coords1, coords2, cache_arrays, &
      split_parts, n_split_parts, current_split_idx, this_map, best_map, min_dist)
   ! Recursively try all assignments for all split parts
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(in) :: split_parts(:), n_split_parts, current_split_idx
   type(partmap_t), intent(inout) :: this_map, best_map
   real(rk), intent(inout) :: min_dist

   integer(ik) :: split_part_idx, child_branch_idx, first_link_idx
   integer(ik) :: n_items2, i, j
   integer(ik) :: branch_link_offset, branch_num_links
   integer(ik) :: saved_size
   real(rk) :: total_dist

   ! Base case: we've made assignments for all split parts - evaluate complete permutation
   if (current_split_idx > n_split_parts) then
      n_combinations = n_combinations + 1

      ! Calculate least distance of this permutation
      total_dist = least_sqdistsum(this_map%subset(:this_map%subset_size), this_map%mapping, coords1, coords2)

      ! Update global best if this permutation is better
      if (total_dist < min_dist) then
         min_dist = total_dist
         best_map%subset_size = 0
         call submap_merge(best_map, this_map)
      end if
      return
   end if

   ! Get the current split part
   split_part_idx = split_parts(current_split_idx)

   ! Find the child branch that uses this split part
   child_branch_idx = 0
   do i = 1, cache_arrays%total_chains
      if (cache_arrays%assigntree(i)%split_part_idx == split_part_idx) then
         child_branch_idx = i
         exit
      end if
   end do

   if (child_branch_idx == 0) then
      error stop 'Could not find child branch for split part'
   end if

   first_link_idx = cache_arrays%assigntree(child_branch_idx)%link_offset + 1
   branch_link_offset = cache_arrays%assigntree(child_branch_idx)%link_offset
   branch_num_links = cache_arrays%assigntree(child_branch_idx)%n_links

   n_items2 = cache_arrays%partree(split_part_idx)%n_items2

   ! Try all possible assignments for this split part (only vary item2, keep item1 at index 1)
   do j = 1, n_items2
      ! Save current permutation state (pairs are only appended, so the
      ! size is enough to roll back)
      saved_size = this_map%subset_size

      ! Make assignment for this split part (always use first item1, index=1)
      call assign_branch_atoms(cache_arrays, split_part_idx, child_branch_idx, &
            first_link_idx, 1, j, this_map)

      ! Recursively try assignments for remaining split parts
      call recur_assign_atoms_global(coords1, coords2, cache_arrays, split_parts, &
            n_split_parts, current_split_idx + 1, this_map, best_map, min_dist)

      ! Restore permutation state
      this_map%subset_size = saved_size

      ! Reset itemdir state for this branch
      cache_arrays%itemdir1_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
      cache_arrays%itemdir2_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
   end do
end subroutine

subroutine assign_atoms_global(coords1, coords2, cache_arrays, mapping1)
   ! DFS exploration of all permutations
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), dimension(:), allocatable, intent(out) :: mapping1
   ! Local variables
   type(partmap_t) :: this_map, best_map
   real(rk) :: min_dist
   integer(ik), allocatable :: split_parts(:)
   integer(ik) :: n_split_parts, n_atoms

   n_atoms = cache_arrays%n_atoms1

   ! Initialize permutations
   call submap_init(this_map, n_atoms)
   call submap_init(best_map, n_atoms)

   ! Initialize both permutations with preassigned pairs
   call collect_leaf_assignments(cache_arrays, 1, this_map)
   call collect_leaf_assignments(cache_arrays, 1, best_map)

   ! Initialize global minimum distance
   min_dist = huge(min_dist)

   ! Initialize DFS exploration variables
   n_combinations = 0

   ! Collect all split parts
   call collect_split_parts(cache_arrays, split_parts, n_split_parts)

   ! Try all assignments for all split parts
   call recur_assign_atoms_global(coords1, coords2, cache_arrays, split_parts, &
         n_split_parts, 1, this_map, best_map, min_dist)

   ! Clean up
   deallocate(split_parts)

   ! Convert partmap_t to permutation array
   call get_mapping(best_map, mapping1)
end subroutine

recursive subroutine recur_assign_atoms_local(coords1, coords2, cache_arrays, &
      branch_idx, best_map, accumulated_dist)
   ! DFS exploration of all assignment possibilities - finds permutation that minimizes total distance
   ! OPTIMIZED: Incremental distance calculation to avoid redundant O(n) sqdistsum calls
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(in) :: branch_idx
   type(partmap_t), intent(inout) :: best_map
   real(rk), intent(inout) :: accumulated_dist  ! NEW: incrementally track distance

   integer(ik) :: child_branch_idx, first_link_idx, split_part_idx, i, n_items2, j
   integer(ik) :: branch_link_offset, branch_num_links
   type(partmap_t) :: best_branch_map, branch_map
   real(rk) :: branch_dist, min_branch_dist
   integer(ik) :: n_atoms

   n_atoms = cache_arrays%n_atoms1

   ! Check if this is a leaf level (no more child branches)
   if (cache_arrays%assigntree(branch_idx)%n_children == 0) then
      n_combinations = n_combinations + 1
      return
   end if

   ! Allocate the branch subperms once per call; they are emptied per child
   call submap_init(best_branch_map, n_atoms)
   call submap_init(branch_map, n_atoms)

   ! Process each child branch independently
   do i = 1, cache_arrays%assigntree(branch_idx)%n_children
      child_branch_idx = cache_arrays%assigntree(branch_idx)%child_indices(i)
      first_link_idx = cache_arrays%assigntree(child_branch_idx)%link_offset + 1
      split_part_idx = cache_arrays%assigntree(child_branch_idx)%split_part_idx
      branch_link_offset = cache_arrays%assigntree(child_branch_idx)%link_offset
      branch_num_links = cache_arrays%assigntree(child_branch_idx)%n_links

      n_items2 = cache_arrays%partree(split_part_idx)%n_items2
      min_branch_dist = huge(min_branch_dist)

      best_branch_map%subset_size = 0
      branch_map%subset_size = 0

      ! Try pairing first item1 with each item2 to find best assignment for this branch
      do j = 1, n_items2
         ! Reset branch assignment and distance for this iteration
         branch_map%subset_size = 0
         branch_dist = 0

         ! Make assignment and recompute HNAs (using first item1, index=1)
         call assign_branch_atoms(cache_arrays, split_part_idx, child_branch_idx, &
               first_link_idx, 1, j, branch_map)

         ! Add new assigned pairs distance to branch distance
         branch_dist = branch_dist + &
               sqdistsum(branch_map%subset(:branch_map%subset_size), branch_map%mapping, coords1, coords2)

         ! Recursively explore subtree - distance is accumulated in branch_dist
         call recur_assign_atoms_local(coords1, coords2, cache_arrays, &
               child_branch_idx, branch_map, branch_dist)

         ! branch_dist now contains total accumulated distance - no recalculation needed!
         if (branch_dist < min_branch_dist) then
            min_branch_dist = branch_dist
            best_branch_map%subset_size = 0
            call submap_merge(best_branch_map, branch_map)
         end if

         ! Reset state for next iteration - only reset links used by this branch
         cache_arrays%itemdir1_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
         cache_arrays%itemdir2_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
      end do

      ! Update the optimal assignment with the best assignment from this branch
      call submap_merge(best_map, best_branch_map)

      ! Accumulate the best distance from this branch
      accumulated_dist = accumulated_dist + min_branch_dist
   end do
end subroutine

subroutine assign_atoms_local(coords1, coords2, cache_arrays, mapping1, total_dist)
   ! DFS exploration wrapper - finds optimal assignment among all possibilities
   ! OPTIMIZED: Uses incremental distance calculation
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), dimension(:), allocatable, intent(out) :: mapping1
   real(rk), intent(out) :: total_dist
   ! Local variables
   type(partmap_t) :: best_map
   integer(ik) :: n_atoms

   n_atoms = cache_arrays%n_atoms1

   ! Initialize optimal assignment
   call submap_init(best_map, n_atoms)

   ! Initialize optimal assignment with preassigned pairs
   call collect_leaf_assignments(cache_arrays, 1, best_map)

   ! Initialize DFS exploration variables
   n_combinations = 0

   ! Initialize distance accumulator with preassigned pairs
   total_dist = sqdistsum(best_map%subset(:best_map%subset_size), best_map%mapping, coords1, coords2)

   ! Perform DFS exploration to find optimal assignment (starting from root chain at index 1)
   call recur_assign_atoms_local(coords1, coords2, cache_arrays, 1, best_map, total_dist)

   ! Convert partmap_t to permutation array
   call get_mapping(best_map, mapping1)
end subroutine

function estimate_unassigned_lower_bound(coords1, coords2, cache_arrays, branch_idx) result(lower_bound)
   ! Estimate lower bound distance for all unassigned atoms using nearest-neighbor approach
   ! For each unassigned atom in molecule 1, finds closest unassigned atom in molecule 2
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: branch_idx
   real(rk) :: dist, min_dist, lower_bound
   integer(ik) :: child_idx, split_part_idx, item1_idx, item2_idx
   integer(ik) :: items1_offset, items2_offset, n_items1, n_items2
   integer(ik) :: i, j, k
   
   lower_bound = 0.0_rk
   
   ! Process all child branches to accumulate unassigned atoms
   do i = 1, cache_arrays%assigntree(branch_idx)%n_children
      child_idx = cache_arrays%assigntree(branch_idx)%child_indices(i)
      split_part_idx = cache_arrays%assigntree(child_idx)%split_part_idx
      
      items1_offset = cache_arrays%partree(split_part_idx)%items1_offset
      n_items1 = cache_arrays%partree(split_part_idx)%n_items1
      items2_offset = cache_arrays%partree(split_part_idx)%items2_offset
      n_items2 = cache_arrays%partree(split_part_idx)%n_items2
      
      ! For each atom in molecule 1 at this split, find nearest atom in molecule 2
      do j = 1, n_items1
         item1_idx = cache_arrays%atomidcs1(items1_offset + j)
         min_dist = huge(min_dist)
         
         ! Find minimum distance to any atom in molecule 2 at this split
         do k = 1, n_items2
            item2_idx = cache_arrays%atomidcs2(items2_offset + k)
            dist = sum((coords1(:, item1_idx) - coords2(:, item2_idx))**2)
            
            if (dist < min_dist) then
               min_dist = dist
            end if
         end do
         
         lower_bound = lower_bound + min_dist
      end do
   end do
end function

recursive subroutine recur_assign_atoms_local_pruned(coords1, coords2, cache_arrays, &
         branch_idx, total_budget, best_map, accumulated_dist, success)
   ! DFS exploration with pruning_atoms - finds permutation that minimizes total distance
   ! OPTIMIZED: Incremental distance calculation + nearest-neighbor lower bound pruning
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(in) :: branch_idx
   real(rk), intent(in) :: total_budget
   type(partmap_t), intent(inout) :: best_map
   real(rk), intent(inout) :: accumulated_dist  ! NEW: incrementally track distance
   logical(lk), intent(inout) :: success

   integer(ik) :: child_branch_idx, first_link_idx, split_part_idx, i, n_items2, j
   integer(ik) :: branch_link_offset, branch_num_links
   type(partmap_t) :: best_branch_map, branch_map
   real(rk) :: branch_dist, min_branch_dist
   logical(lk) :: branch_success, child_success
   real(rk) :: remaining_budget, lower_bound_estimate
   integer(ik) :: n_atoms

   n_atoms = cache_arrays%n_atoms1

   ! Check if this is a leaf level (no more child branches)
   if (cache_arrays%assigntree(branch_idx)%n_children == 0) then
      n_combinations = n_combinations + 1
      success = .TRUE.
      return
   end if

   ! Calculate remaining budget at this level
   remaining_budget = total_budget - accumulated_dist
   
   ! AGGRESSIVE PRUNING: Estimate lower bound for unassigned atoms
   lower_bound_estimate = estimate_unassigned_lower_bound(coords1, coords2, cache_arrays, branch_idx)
   
   ! If even the best-case scenario exceeds budget, prune this branch
   if (lower_bound_estimate >= remaining_budget) then
      success = .FALSE.
      return
   end if

   ! Allocate the branch subperms once per call; they are emptied per child
   call submap_init(best_branch_map, n_atoms)
   call submap_init(branch_map, n_atoms)

   ! Process each child branch independently
   do i = 1, cache_arrays%assigntree(branch_idx)%n_children
      ! Early exit: if no threshold budget remains, remaining branches cannot succeed
      if (remaining_budget < 0) then
         success = .FALSE.
         return
      end if

      child_branch_idx = cache_arrays%assigntree(branch_idx)%child_indices(i)
      first_link_idx = cache_arrays%assigntree(child_branch_idx)%link_offset + 1
      split_part_idx = cache_arrays%assigntree(child_branch_idx)%split_part_idx
      branch_link_offset = cache_arrays%assigntree(child_branch_idx)%link_offset
      branch_num_links = cache_arrays%assigntree(child_branch_idx)%n_links

      n_items2 = cache_arrays%partree(split_part_idx)%n_items2
      min_branch_dist = huge(min_branch_dist)  ! Best distance for this specific branch
      branch_success = .FALSE.

      best_branch_map%subset_size = 0
      branch_map%subset_size = 0

      ! Try pairing first item1 with each item2 to find best assignment for this branch
      do j = 1, n_items2
         ! Reset branch assignment for this iteration
         branch_map%subset_size = 0
         branch_dist = 0

         ! Make assignment and recompute HNAs operation (using first item1, index=1)
         call assign_branch_atoms(cache_arrays, split_part_idx, child_branch_idx, &
               first_link_idx, 1, j, branch_map)

         ! Add new assigned pairs distance to branch distance
         branch_dist = branch_dist + &
               sqdistsum(branch_map%subset(:branch_map%subset_size), branch_map%mapping, coords1, coords2)

         ! PRUNING: Early check - only continue if current partial distance is within remaining threshold
         if (branch_dist < remaining_budget) then
            child_success = .FALSE.

            ! Recursively explore subtree - distance accumulates in branch_dist
            call recur_assign_atoms_local_pruned(coords1, coords2, cache_arrays, &
                  child_branch_idx, total_budget, branch_map, branch_dist, child_success)

            ! Only consider this branch if the recursive call succeeded
            if (child_success) then
               ! branch_dist now contains total accumulated distance - no recalculation needed!
               if (branch_dist < min_branch_dist) then
                  min_branch_dist = branch_dist
                  best_branch_map%subset_size = 0
                  call submap_merge(best_branch_map, branch_map)
                  branch_success = .TRUE.
               end if
            end if
         end if

         ! Reset state for next iteration - only reset links used by this branch
         cache_arrays%itemdir1_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
         cache_arrays%itemdir2_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
      end do

      ! If this branch failed to find a solution, the entire recursion fails
      if (.not. branch_success) then
         success = .FALSE.
         return
      end if

      ! Update the optimal assignment with the best assignment from this branch
      call submap_merge(best_map, best_branch_map)

      ! Accumulate the best distance from this branch
      accumulated_dist = accumulated_dist + min_branch_dist

      ! Update remaining budget (subtract this branch's contribution)
      remaining_budget = remaining_budget - min_branch_dist
   end do

   ! If we get here, check if we're within budget
   success = (remaining_budget >= 0)
end subroutine

subroutine assign_atoms_local_pruned(coords1, coords2, cache_arrays, mapping1, mapdist)
   ! DFS exploration with pruning_atoms threshold - finds assignment within threshold
   ! OPTIMIZED: Uses incremental distance calculation
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), dimension(:), allocatable, intent(out) :: mapping1
   real(rk), intent(inout) :: mapdist
   ! Local variables
   real(rk) :: total_budget
   type(partmap_t) :: best_map
   logical(lk) :: success
   integer(ik) :: n_atoms

   total_budget = mapdist + SQDIST_TOL
   n_atoms = cache_arrays%n_atoms1

   ! Initialize optimal assignment
   call submap_init(best_map, n_atoms)

   ! Initialize assignment with preassigned pairs
   call collect_leaf_assignments(cache_arrays, 1, best_map)

   ! Initialize accumulated distance with preassigned pairs
   mapdist = sqdistsum(best_map%subset(:best_map%subset_size), best_map%mapping, coords1, coords2)

   ! Initialize DFS exploration variables
   success = .FALSE.
   n_combinations = 0

   ! Perform pruned DFS exploration (starting from root chain at index 1)
   call recur_assign_atoms_local_pruned(coords1, coords2, cache_arrays, 1, total_budget, &
         best_map, mapdist, success)

   if (.not. success) error stop 'Assignment failed'

   ! Convert partmap_t to permutation array
   call get_mapping(best_map, mapping1)
end subroutine

end module
