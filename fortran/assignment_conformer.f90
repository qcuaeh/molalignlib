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
! Symmetry-corrected atom assignment between conformers at a fixed relative
! orientation, by traversal of the assignment tree built from the
! self-consistent HNA partition (see refinement::build_assignment_tree and
! J. Chem. Theory Comput., doi:10.1021/acs.jctc.6c00545).
!
! Each assignment tree node (chain) is a branch point that owns a split
! part: a leaf of the HNA partition holding several topologically
! equivalent atom pairs. Assigning one pair of the split part
! (individualization) and replaying the refinement recorded in the links of
! the chain fixes the pairs of every part that becomes a singleton. The
! children of a node are independent subproblems, so their optimal
! assignments are found separately and combined, which evaluates the sum of
! the branch possibilities instead of their product. Since the atoms of a
! part are equivalent, the exact searches always pair the first atom of
! molecule 1 in the split part and only vary its partner in molecule 2.
!
! All procedures work on the array representation of the trees
! (array_trees_t). The per-link vertex directories (itemdir1/2_entries) are
! rewritten while a pair choice is explored and cleared afterwards.
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

! Maximum number of children of a part
integer(ik), parameter :: MAX_CHILD = 10

! Complete assignments (global search) or tree leaves (local searches)
! visited by the last search; diagnostic counter
integer(ik) :: n_combinations

contains

function signature_equivalence_array(cache_arrays, part_idx, signature_size, signature) result(equiv)
! Multiset equality between signature and the cached signature of part
! part_idx, stored as its distinct values with their frequencies
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: part_idx
   integer(ik), intent(in) :: signature_size, signature(:)
   logical(lk) :: equiv
   integer(ik) :: signature_frequencies, i, j

   if (signature_size /= cache_arrays%partree(part_idx)%signature_size) then
      equiv = .FALSE.
      return
   end if

   ! Single-entry signatures (the most common case) are compared directly
   if (cache_arrays%partree(part_idx)%signature_size == 1) then
      equiv = (signature(1) == cache_arrays%partree(part_idx)%signature_values(1))
      return
   end if

   ! Same size and same frequency of every cached value
   do i = 1, cache_arrays%partree(part_idx)%signature_unique_count
      signature_frequencies = 0
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
! Position (1..n_children) of the child of part parent_idx whose signature
! equals signature. The search replays the refinement that created the
! children, so a match must exist.
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: parent_idx
   integer(ik), intent(in) :: signature_size, signature(:)
   integer(ik) :: child_relative_idx
   integer(ik) :: n_children, child_idx, i

   n_children = cache_arrays%partree(parent_idx)%n_children

   do i = 1, n_children
      child_idx = cache_arrays%partree(parent_idx)%child_indices(i)
      if (signature_equivalence_array(cache_arrays, child_idx, signature_size, signature)) then
         child_relative_idx = i
         return
      end if
   end do

   error stop 'Part signature did not match any child'
end function

subroutine get_mapping(atom_map, mapping1)
! Export a complete atom_map as a permutation array. Once a search is over
! every atom of the tree has been assigned.
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
! Add to atom_map the pair of every child of part part_idx that is a
! singleton leaf (one atom of each molecule), i.e. the pairs fixed by the
! last subdivision of that part
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: part_idx
   type(partmap_t), intent(inout) :: atom_map
   integer(ik) :: i, child_idx, item1_idx, item2_idx

   do i = 1, cache_arrays%partree(part_idx)%n_children
      child_idx = cache_arrays%partree(part_idx)%child_indices(i)
      if (cache_arrays%partree(child_idx)%n_children == 0) then
         if (cache_arrays%partree(child_idx)%n_items1 == 1 .and. &
             cache_arrays%partree(child_idx)%n_items2 == 1) then

            item1_idx = cache_arrays%atomidcs1(cache_arrays%partree(child_idx)%items1_offset + 1)
            item2_idx = cache_arrays%atomidcs2(cache_arrays%partree(child_idx)%items2_offset + 1)
            call submap_add(atom_map, item1_idx, item2_idx)
         end if
      end if
   end do
end subroutine

subroutine update_hna_part(cache_arrays, part_idx, read_link_idx, write_link_idx, atom_map)
! Refine part part_idx for the current pair choices: every atom goes to the
! child whose cached signature matches its signature under the vertex
! directory of link read_link_idx, and its new part is written to the
! directory of link write_link_idx. The children and their sizes were fixed
! when the tree was built; only the atoms they hold change. Only neighbors
! with a part in the read directory, i.e. those moved at the previous level,
! enter the signature; the others contribute equally to all atoms of a
! self-consistent part.
! Pairs of new singleton children are added to atom_map.
! The work arrays are saved, so this routine is not thread safe.
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(in) :: part_idx, read_link_idx, write_link_idx
   type(partmap_t), intent(inout) :: atom_map
   integer(ik), save :: signature_size, signature(MAX_COORDNUM)
   integer(ik), dimension(MAX_CHILD), save :: items1_trackers, items2_trackers
   integer(ik) :: target_relative_idx, target_part_idx, item_idx, target_idx, part_ref_idx
   integer(ik) :: items1_offset, n_items1, items2_offset, n_items2, n_children
   integer(ik) :: adj_atom
   integer(ik) :: i, j

   items1_offset = cache_arrays%partree(part_idx)%items1_offset
   n_items1 = cache_arrays%partree(part_idx)%n_items1
   items2_offset = cache_arrays%partree(part_idx)%items2_offset
   n_items2 = cache_arrays%partree(part_idx)%n_items2
   n_children = cache_arrays%partree(part_idx)%n_children

   ! Number of atoms placed so far in each child
   do i = 1, n_children
      items1_trackers(i) = 0
      items2_trackers(i) = 0
   end do

   ! Atoms of molecule 1
   do i = 1, n_items1
      item_idx = cache_arrays%atomidcs1(items1_offset + i)

      ! Signature of the atom
      signature_size = 0
      do j = 1, cache_arrays%adjcs1_cn(item_idx)
         adj_atom = cache_arrays%adjcs1_list(item_idx, j)
         part_ref_idx = cache_arrays%itemdir1_entries(adj_atom, read_link_idx)
         if (part_ref_idx /= 0) then
            signature_size = signature_size + 1
            ! Inlined adjacency::edge_code
            signature(signature_size) = part_ref_idx*BOND_TYPE_RADIX &
                                      + cache_arrays%adjcs1_bondtype(item_idx, j)
         end if
      end do

      ! Move the atom to the child with the same signature
      target_relative_idx = find_child_part_array(cache_arrays, part_idx, signature_size, signature)
      target_part_idx = cache_arrays%partree(part_idx)%child_indices(target_relative_idx)
      items1_trackers(target_relative_idx) = items1_trackers(target_relative_idx) + 1
      target_idx = cache_arrays%partree(target_part_idx)%items1_offset &
                 + items1_trackers(target_relative_idx)
      cache_arrays%atomidcs1(target_idx) = item_idx

      cache_arrays%itemdir1_entries(item_idx, write_link_idx) = target_part_idx
   end do

   ! Atoms of molecule 2
   do i = 1, n_items2
      item_idx = cache_arrays%atomidcs2(items2_offset + i)

      ! Signature of the atom
      signature_size = 0
      do j = 1, cache_arrays%adjcs2_cn(item_idx)
         adj_atom = cache_arrays%adjcs2_list(item_idx, j)
         part_ref_idx = cache_arrays%itemdir2_entries(adj_atom, read_link_idx)
         if (part_ref_idx /= 0) then
            signature_size = signature_size + 1
            ! Inlined adjacency::edge_code
            signature(signature_size) = part_ref_idx*BOND_TYPE_RADIX &
                                      + cache_arrays%adjcs2_bondtype(item_idx, j)
         end if
      end do

      ! Move the atom to the child with the same signature
      target_relative_idx = find_child_part_array(cache_arrays, part_idx, signature_size, signature)
      target_part_idx = cache_arrays%partree(part_idx)%child_indices(target_relative_idx)
      items2_trackers(target_relative_idx) = items2_trackers(target_relative_idx) + 1
      target_idx = cache_arrays%partree(target_part_idx)%items2_offset &
                 + items2_trackers(target_relative_idx)
      cache_arrays%atomidcs2(target_idx) = item_idx

      cache_arrays%itemdir2_entries(item_idx, write_link_idx) = target_part_idx
   end do

   call collect_leaf_assignments(cache_arrays, part_idx, atom_map)
end subroutine

subroutine assign_pair_to_children(cache_arrays, split_part_idx, first_link_idx, &
      chosen_item1_idx, chosen_item2_idx, atom_map)
! Individualize a pair of split part split_part_idx: the atoms at positions
! chosen_item1_idx and chosen_item2_idx of the part go to its first child
! and the remaining atoms to its second child, in the vertex directory of
! link first_link_idx. The new pair is added to atom_map.
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(in) :: split_part_idx, first_link_idx, chosen_item1_idx, chosen_item2_idx
   type(partmap_t), intent(inout) :: atom_map
   integer(ik) :: child_part1, child_part2, chosen_item1, chosen_item2, i, item_idx, target_idx
   integer(ik) :: items1_offset, n_items1, items2_offset, n_items2

   items1_offset = cache_arrays%partree(split_part_idx)%items1_offset
   n_items1 = cache_arrays%partree(split_part_idx)%n_items1
   items2_offset = cache_arrays%partree(split_part_idx)%items2_offset
   n_items2 = cache_arrays%partree(split_part_idx)%n_items2

   child_part1 = cache_arrays%partree(split_part_idx)%child_indices(1)
   child_part2 = cache_arrays%partree(split_part_idx)%child_indices(2)

   chosen_item1 = cache_arrays%atomidcs1(items1_offset + chosen_item1_idx)
   chosen_item2 = cache_arrays%atomidcs2(items2_offset + chosen_item2_idx)

   ! Chosen pair to the first child
   cache_arrays%atomidcs1(cache_arrays%partree(child_part1)%items1_offset + 1) = chosen_item1
   cache_arrays%atomidcs2(cache_arrays%partree(child_part1)%items2_offset + 1) = chosen_item2
   cache_arrays%itemdir1_entries(chosen_item1, first_link_idx) = child_part1
   cache_arrays%itemdir2_entries(chosen_item2, first_link_idx) = child_part1

   ! Remaining atoms to the second child
   target_idx = cache_arrays%partree(child_part2)%items1_offset
   do i = 1, n_items1
      if (i /= chosen_item1_idx) then
         item_idx = cache_arrays%atomidcs1(items1_offset + i)
         target_idx = target_idx + 1
         cache_arrays%atomidcs1(target_idx) = item_idx
         cache_arrays%itemdir1_entries(item_idx, first_link_idx) = child_part2
      end if
   end do

   target_idx = cache_arrays%partree(child_part2)%items2_offset
   do i = 1, n_items2
      if (i /= chosen_item2_idx) then
         item_idx = cache_arrays%atomidcs2(items2_offset + i)
         target_idx = target_idx + 1
         cache_arrays%atomidcs2(target_idx) = item_idx
         cache_arrays%itemdir2_entries(item_idx, first_link_idx) = child_part2
      end if
   end do

   call collect_leaf_assignments(cache_arrays, split_part_idx, atom_map)
end subroutine

subroutine assign_branch_atoms(cache_arrays, split_part_idx, child_branch_idx, &
      first_link_idx, chosen_item1_idx, chosen_item2_idx, atom_map)
! Assign one pair of the split part of node child_branch_idx and propagate
! it through the links of the node (Pseudocode S3, lines 9-12): the parts
! listed in link i are refined reading the vertex directory of link i and
! writing that of link i+1. Every pair fixed on the way is added to
! atom_map.
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(in) :: split_part_idx, child_branch_idx, first_link_idx
   integer(ik), intent(in) :: chosen_item1_idx, chosen_item2_idx
   type(partmap_t), intent(inout) :: atom_map
   integer(ik) :: link_idx, next_link_idx, part_idx
   integer(ik) :: n_links, n_parts, link_offset, partref_offset
   integer(ik) :: i, j

   call assign_pair_to_children(cache_arrays, split_part_idx, first_link_idx, &
         chosen_item1_idx, chosen_item2_idx, atom_map)

   n_links = cache_arrays%assigntree(child_branch_idx)%n_links
   link_offset = cache_arrays%assigntree(child_branch_idx)%link_offset

   if (DEBUG_TESTS) then
      ! The searches only clear links link_offset+1..link_offset+n_links, so
      ! the last link must have no parts, or its writes to the next link
      ! would survive the reset
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
! Split parts of all assignment tree nodes, in node order (parents before
! children), for the exhaustive search
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), allocatable, intent(out) :: split_parts(:)
   integer(ik), intent(out) :: n_split_parts
   integer(ik) :: i, count

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
! Greedy descent of the assignment tree: every split part is assigned its
! closest pair of atoms. The result is a valid assignment whose distance
! is an upper bound for the pruned search.
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

   if (cache_arrays%assigntree(branch_idx)%n_children == 0) then
      return
   end if

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

      ! Closest pair of the split part
      min_dist = huge(min_dist)
      greedy_idx1 = 1
      greedy_idx2 = 1

      do idx1 = 1, n_items1
         item1_idx = cache_arrays%atomidcs1(items1_offset + idx1)
         do idx2 = 1, n_items2
            item2_idx = cache_arrays%atomidcs2(items2_offset + idx2)

            current_dist = sum((coords1(:, item1_idx) - coords2(:, item2_idx))**2)

            if (current_dist < min_dist) then
               min_dist = current_dist
               greedy_idx1 = idx1
               greedy_idx2 = idx2
            end if
         end do
      end do

      call assign_branch_atoms(cache_arrays, split_part_idx, child_branch_idx, &
            first_link_idx, greedy_idx1, greedy_idx2, greedy_map)
      call recur_assign_atoms_greedy(coords1, coords2, cache_arrays, child_branch_idx, &
            greedy_map)

      ! Clear the vertex directories of the links of this branch
      cache_arrays%itemdir1_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
      cache_arrays%itemdir2_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
   end do
end subroutine

subroutine assign_atoms_greedy(coords1, coords2, cache_arrays, mapping1, mapdist)
! Greedy assignment (see recur_assign_atoms_greedy); mapdist is its sum of
! squared distances
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), dimension(:), allocatable, intent(out) :: mapping1
   real(rk), intent(out) :: mapdist
   ! Local variables
   type(partmap_t) :: greedy_map

   call submap_init(greedy_map, cache_arrays%n_atoms1)

   ! Pairs already fixed by the self-consistent partition
   call collect_leaf_assignments(cache_arrays, 1, greedy_map)

   ! Descend from the root node
   call recur_assign_atoms_greedy(coords1, coords2, cache_arrays, 1, greedy_map)

   mapdist = sqdistsum(greedy_map%subset(:greedy_map%subset_size), greedy_map%mapping, coords1, coords2)

   call get_mapping(greedy_map, mapping1)
end subroutine

recursive subroutine recur_assign_atoms_global(coords1, coords2, cache_arrays, &
      split_parts, n_split_parts, current_split_idx, this_map, best_map, min_dist)
! Exhaustive enumeration of complete assignments: every pair choice of split
! part current_split_idx is combined with every choice of the following
! ones, and each complete assignment is scored by its sum of squared
! distances after optimal superposition. The cost is the product of the
! branch possibilities, so this is only used when that product is small
! compared with their sum (adaptive alignment strategy).
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

   ! All split parts assigned: score the complete assignment
   if (current_split_idx > n_split_parts) then
      n_combinations = n_combinations + 1
      total_dist = least_sqdistsum(this_map%subset(:this_map%subset_size), this_map%mapping, coords1, coords2)
      if (total_dist < min_dist) then
         min_dist = total_dist
         best_map%subset_size = 0
         call submap_merge(best_map, this_map)
      end if
      return
   end if

   split_part_idx = split_parts(current_split_idx)

   ! Node that owns this split part
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

   ! Pair the first atom of molecule 1 with each atom of molecule 2
   do j = 1, n_items2
      ! Pairs are only appended, so the size is enough to roll back
      saved_size = this_map%subset_size

      call assign_branch_atoms(cache_arrays, split_part_idx, child_branch_idx, &
            first_link_idx, 1, j, this_map)
      call recur_assign_atoms_global(coords1, coords2, cache_arrays, split_parts, &
            n_split_parts, current_split_idx + 1, this_map, best_map, min_dist)

      this_map%subset_size = saved_size

      ! Clear the vertex directories of the links of this branch
      cache_arrays%itemdir1_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
      cache_arrays%itemdir2_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
   end do
end subroutine

subroutine assign_atoms_global(coords1, coords2, cache_arrays, mapping1)
! Orientation-independent exhaustive search: the valid assignment with the
! lowest RMSD after optimal superposition (see recur_assign_atoms_global)
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), dimension(:), allocatable, intent(out) :: mapping1
   ! Local variables
   type(partmap_t) :: this_map, best_map
   real(rk) :: min_dist
   integer(ik), allocatable :: split_parts(:)
   integer(ik) :: n_split_parts, n_atoms

   n_atoms = cache_arrays%n_atoms1

   call submap_init(this_map, n_atoms)
   call submap_init(best_map, n_atoms)

   ! Pairs already fixed by the self-consistent partition
   call collect_leaf_assignments(cache_arrays, 1, this_map)
   call collect_leaf_assignments(cache_arrays, 1, best_map)

   min_dist = huge(min_dist)
   n_combinations = 0

   call collect_split_parts(cache_arrays, split_parts, n_split_parts)
   call recur_assign_atoms_global(coords1, coords2, cache_arrays, split_parts, &
         n_split_parts, 1, this_map, best_map, min_dist)
   deallocate(split_parts)

   call get_mapping(best_map, mapping1)
end subroutine

recursive subroutine recur_assign_atoms_local(coords1, coords2, cache_arrays, &
      branch_idx, best_map, accumulated_dist)
! Depth-first search for the assignment of the subtree of node branch_idx
! that minimizes the sum of squared distances at fixed orientation
! (Pseudocode S3). The children of a node are independent: for each child,
! every pair choice of its split part is tried and its subtree solved
! recursively, and the best choice is merged into best_map. The squared
! distances of the best choices are added to accumulated_dist.
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(in) :: branch_idx
   type(partmap_t), intent(inout) :: best_map
   real(rk), intent(inout) :: accumulated_dist

   integer(ik) :: child_branch_idx, first_link_idx, split_part_idx, i, n_items2, j
   integer(ik) :: branch_link_offset, branch_num_links
   type(partmap_t) :: best_branch_map, branch_map
   real(rk) :: branch_dist, min_branch_dist
   integer(ik) :: n_atoms

   n_atoms = cache_arrays%n_atoms1

   if (cache_arrays%assigntree(branch_idx)%n_children == 0) then
      n_combinations = n_combinations + 1
      return
   end if

   ! Allocated once per call and emptied for each child
   call submap_init(best_branch_map, n_atoms)
   call submap_init(branch_map, n_atoms)

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

      ! Pair the first atom of molecule 1 with each atom of molecule 2
      do j = 1, n_items2
         branch_map%subset_size = 0
         branch_dist = 0

         ! Individualize the pair and propagate it through the branch
         call assign_branch_atoms(cache_arrays, split_part_idx, child_branch_idx, &
               first_link_idx, 1, j, branch_map)

         ! Distance of the pairs fixed at this level
         branch_dist = branch_dist + &
               sqdistsum(branch_map%subset(:branch_map%subset_size), branch_map%mapping, coords1, coords2)

         ! Add the optimal distance of the subtree
         call recur_assign_atoms_local(coords1, coords2, cache_arrays, &
               child_branch_idx, branch_map, branch_dist)

         if (branch_dist < min_branch_dist) then
            min_branch_dist = branch_dist
            best_branch_map%subset_size = 0
            call submap_merge(best_branch_map, branch_map)
         end if

         ! Clear the vertex directories of the links of this branch
         cache_arrays%itemdir1_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
         cache_arrays%itemdir2_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
      end do

      ! Keep the best choice for this child
      call submap_merge(best_map, best_branch_map)
      accumulated_dist = accumulated_dist + min_branch_dist
   end do
end subroutine

subroutine assign_atoms_local(coords1, coords2, cache_arrays, mapping1, total_dist)
! Optimal assignment at fixed orientation (see recur_assign_atoms_local);
! total_dist is its sum of squared distances
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), dimension(:), allocatable, intent(out) :: mapping1
   real(rk), intent(out) :: total_dist
   ! Local variables
   type(partmap_t) :: best_map
   integer(ik) :: n_atoms

   n_atoms = cache_arrays%n_atoms1

   call submap_init(best_map, n_atoms)

   ! Pairs already fixed by the self-consistent partition
   call collect_leaf_assignments(cache_arrays, 1, best_map)
   n_combinations = 0
   total_dist = sqdistsum(best_map%subset(:best_map%subset_size), best_map%mapping, coords1, coords2)

   ! Search from the root node
   call recur_assign_atoms_local(coords1, coords2, cache_arrays, 1, best_map, total_dist)

   call get_mapping(best_map, mapping1)
end subroutine

function estimate_unassigned_lower_bound(coords1, coords2, cache_arrays, branch_idx) result(lower_bound)
! Lower bound of the distance still to be added below node branch_idx: each
! atom of molecule 1 in the split part of a child contributes its squared
! distance to the closest atom of molecule 2 in the same part
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(in) :: cache_arrays
   integer(ik), intent(in) :: branch_idx
   real(rk) :: dist, min_dist, lower_bound
   integer(ik) :: child_idx, split_part_idx, item1_idx, item2_idx
   integer(ik) :: items1_offset, items2_offset, n_items1, n_items2
   integer(ik) :: i, j, k
   
   lower_bound = 0.0_rk
   
   do i = 1, cache_arrays%assigntree(branch_idx)%n_children
      child_idx = cache_arrays%assigntree(branch_idx)%child_indices(i)
      split_part_idx = cache_arrays%assigntree(child_idx)%split_part_idx
      
      items1_offset = cache_arrays%partree(split_part_idx)%items1_offset
      n_items1 = cache_arrays%partree(split_part_idx)%n_items1
      items2_offset = cache_arrays%partree(split_part_idx)%items2_offset
      n_items2 = cache_arrays%partree(split_part_idx)%n_items2
      
      do j = 1, n_items1
         item1_idx = cache_arrays%atomidcs1(items1_offset + j)
         min_dist = huge(min_dist)
         
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
! Branch-and-bound version of recur_assign_atoms_local. total_budget is an
! upper bound of the optimal total distance; a subtree is abandoned as soon
! as the distance accumulated so far plus a lower bound of the rest
! (estimate_unassigned_lower_bound) reaches it. success is false when no
! assignment within the budget exists below this node.
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: cache_arrays
   integer(ik), intent(in) :: branch_idx
   real(rk), intent(in) :: total_budget
   type(partmap_t), intent(inout) :: best_map
   real(rk), intent(inout) :: accumulated_dist
   logical(lk), intent(inout) :: success

   integer(ik) :: child_branch_idx, first_link_idx, split_part_idx, i, n_items2, j
   integer(ik) :: branch_link_offset, branch_num_links
   type(partmap_t) :: best_branch_map, branch_map
   real(rk) :: branch_dist, min_branch_dist
   logical(lk) :: branch_success, child_success
   real(rk) :: remaining_budget, lower_bound_estimate
   integer(ik) :: n_atoms

   n_atoms = cache_arrays%n_atoms1

   if (cache_arrays%assigntree(branch_idx)%n_children == 0) then
      n_combinations = n_combinations + 1
      success = .TRUE.
      return
   end if

   remaining_budget = total_budget - accumulated_dist
   
   ! Prune if even the lower bound exceeds the budget
   lower_bound_estimate = estimate_unassigned_lower_bound(coords1, coords2, cache_arrays, branch_idx)
   if (lower_bound_estimate >= remaining_budget) then
      success = .FALSE.
      return
   end if

   ! Allocated once per call and emptied for each child
   call submap_init(best_branch_map, n_atoms)
   call submap_init(branch_map, n_atoms)

   do i = 1, cache_arrays%assigntree(branch_idx)%n_children
      ! Budget exhausted by the previous children
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
      min_branch_dist = huge(min_branch_dist)
      branch_success = .FALSE.

      best_branch_map%subset_size = 0
      branch_map%subset_size = 0

      ! Pair the first atom of molecule 1 with each atom of molecule 2
      do j = 1, n_items2
         branch_map%subset_size = 0
         branch_dist = 0

         ! Individualize the pair and propagate it through the branch
         call assign_branch_atoms(cache_arrays, split_part_idx, child_branch_idx, &
               first_link_idx, 1, j, branch_map)

         ! Distance of the pairs fixed at this level
         branch_dist = branch_dist + &
               sqdistsum(branch_map%subset(:branch_map%subset_size), branch_map%mapping, coords1, coords2)

         ! Explore the subtree only while within budget
         if (branch_dist < remaining_budget) then
            child_success = .FALSE.
            call recur_assign_atoms_local_pruned(coords1, coords2, cache_arrays, &
                  child_branch_idx, total_budget, branch_map, branch_dist, child_success)
            if (child_success) then
               if (branch_dist < min_branch_dist) then
                  min_branch_dist = branch_dist
                  best_branch_map%subset_size = 0
                  call submap_merge(best_branch_map, branch_map)
                  branch_success = .TRUE.
               end if
            end if
         end if

         ! Clear the vertex directories of the links of this branch
         cache_arrays%itemdir1_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
         cache_arrays%itemdir2_entries(:, branch_link_offset + 1 : branch_link_offset + branch_num_links) = 0
      end do

      ! Any child without a solution within budget fails the whole node
      if (.not. branch_success) then
         success = .FALSE.
         return
      end if

      ! Keep the best choice for this child and charge it to the budget
      call submap_merge(best_map, best_branch_map)
      accumulated_dist = accumulated_dist + min_branch_dist
      remaining_budget = remaining_budget - min_branch_dist
   end do

   success = (remaining_budget >= 0)
end subroutine

subroutine assign_atoms_local_pruned(coords1, coords2, cache_arrays, mapping1, mapdist)
! Optimal assignment at fixed orientation with branch-and-bound pruning. On
! entry mapdist is an upper bound of the optimal distance (from
! assign_atoms_greedy or a previous assignment); SQDIST_TOL is added so the
! assignment that attains it is not pruned. On exit mapdist is the optimal
! sum of squared distances.
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

   call submap_init(best_map, n_atoms)

   ! Pairs already fixed by the self-consistent partition
   call collect_leaf_assignments(cache_arrays, 1, best_map)
   mapdist = sqdistsum(best_map%subset(:best_map%subset_size), best_map%mapping, coords1, coords2)

   success = .FALSE.
   n_combinations = 0

   ! Search from the root node
   call recur_assign_atoms_local_pruned(coords1, coords2, cache_arrays, 1, total_budget, &
         best_map, mapdist, success)

   if (.not. success) error stop 'Assignment failed'

   call get_mapping(best_map, mapping1)
end subroutine

end module
