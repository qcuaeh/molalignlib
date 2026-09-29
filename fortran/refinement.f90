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

module refinement
! Hierarchical Neighborhood of Atoms (HNA) partitioning of two molecules and
! construction of the assignment tree (Pseudocodes S1 and S2 of
! J. Chem. Theory Comput., doi:10.1021/acs.jctc.6c00545).
!
! The partition is kept as a part tree: refining a leaf part turns it into
! a branch whose children group its atoms by signature. Refinement levels
! are recorded as the links of a chain; each link lists the parts of that
! level and holds the vertex directory (itemdir) of the atoms placed at
! that level. Assignment tree nodes are chains too: their links list the
! parts refined at each level after the individualization of their split
! part, so the search can replay the refinement.
use parameters
use common_types
use random
use chemdata
use adjacency
use molecule
use linked_list_types
use indexed_list_types
use flags
use error_codes
implicit none
private
public refine_hna_part
public refine_hna_partition
public compute_sc_hna_chain
public build_assignment_tree

contains

function hna_signature(adjc, itemdir) result(signature)
! Signature of an atom (SIG in the paper): one edge_code per neighbor with a
! part in itemdir, combining the part index and the bond type, compared as
! a multiset. Neighbors without a part in itemdir are left out, as in
! assignment_conformer::update_hna_part.
   type(adjc_t), intent(in) :: adjc
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir
   integer(ik), dimension(:), allocatable :: signature
   ! Local variables
   integer(ik) :: codes(MAX_COORDNUM)
   integer(ik) :: n_codes, j
   type(partition_node_t), pointer :: neighbor_part

   n_codes = 0
   do j = 1, adjc%cn
      neighbor_part => itemdir(adjc%list(j))%ptr
      if (associated(neighbor_part)) then
         n_codes = n_codes + 1
         codes(n_codes) = edge_code(neighbor_part%global_idx, adjc%bondtype(j))
      end if
   end do

   signature = codes(:n_codes)
end function

subroutine refine_hna_part(adjcs1, adjcs2, itemdir1, itemdir2, part, link)
! Subdivide the leaf part `part` by signature (Pseudocode S1): one child per
! distinct signature, added to link, and every atom placed in its child and
! in the vertex directory of link. A part whose atoms share one signature
! gets a single child.
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir1, itemdir2
   type(partition_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(item_node_t), pointer :: item
   type(partition_node_t), pointer :: child_part
   integer(ik), dimension(:), allocatable :: signature

   ! Atoms of molecule 1
   item => part%first_item1
   do while (associated(item))
      signature = hna_signature(adjcs1(item%idx), itemdir1)
      child_part => find_child_part(part, signature)
      if (.not. associated(child_part)) then
         child_part => new_child_part(part)
         allocate (child_part%signature, source=signature)
         call link_part(link, child_part)
      end if
      call add_new_item1(child_part, item%idx)
      link%itemdir1(item%idx)%ptr => child_part
      item => item%next_item
   end do

   ! Atoms of molecule 2
   item => part%first_item2
   do while (associated(item))
      signature = hna_signature(adjcs2(item%idx), itemdir2)
      child_part => find_child_part(part, signature)
      if (.not. associated(child_part)) then
         child_part => new_child_part(part)
         allocate (child_part%signature, source=signature)
         call link_part(link, child_part)
      end if
      call add_new_item2(child_part, item%idx)
      link%itemdir2(item%idx)%ptr => child_part
      item => item%next_item
   end do
end subroutine

subroutine refine_hna_partition(adjcs1, adjcs2, hna_chain, n_splits)
! One refinement step of the whole partition: every part of the last link of
! hna_chain is subdivided into a new link. n_splits is the number of parts
! gained, zero once the partition is self-consistent.
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chaintree_node_t), pointer, intent(inout) :: hna_chain
   integer(ik), intent(out) :: n_splits
   ! Local variables
   type(chain_node_t), pointer :: link, new_link
   type(partref_node_t), pointer :: partref

   n_splits = 0

   link => hna_chain%last_link
   new_link => new_chain_link(hna_chain)

   partref => link%first_partref
   do while (associated(partref))
      call refine_hna_part(adjcs1, adjcs2, link%itemdir1, link%itemdir2, partref%part, new_link)
      n_splits = n_splits + partref%part%n_children - 1
      partref => partref%nextref
   end do
end subroutine

subroutine compute_sc_hna_chain(adjcs1, adjcs2, atomtypes, hna_chain, error_code)
! Self-consistent HNA partition of both molecules: starting from the atom
! types, refine until no part splits. The result is the last link of
! hna_chain. Fails with MOLALIGN_ERROR_NOT_CONFORMERS when a part holds
! different numbers of atoms of each molecule (non-isomorphic graphs).
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(partition_t), intent(in) :: atomtypes
   type(chaintree_node_t), pointer, intent(out) :: hna_chain
   integer(ik), intent(out) :: error_code
   ! Local variables
   integer(ik) :: n_splits

   error_code = MOLALIGN_SUCCESS

   hna_chain => chain_from_partition( atomtypes)

   do
      call refine_hna_partition(adjcs1, adjcs2, hna_chain, n_splits)
      if (n_splits == 0) exit
   end do

   if (is_partition_uneven(hna_chain%last_link)) then
      error_code = MOLALIGN_ERROR_NOT_CONFORMERS
      return
   end if
end subroutine

subroutine update_hna_part(adjcs1, adjcs2, itemdir1, itemdir2, part, link)
! Linked-list counterpart of assignment_conformer::update_hna_part:
! redistribute the atoms of part among its existing children by signature,
! overwriting their item nodes instead of allocating new ones
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir1, itemdir2
   type(partition_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(item_node_t), pointer :: item
   integer(ik), dimension(:), allocatable :: signature
   type(partition_node_t), pointer :: child_part

   ! last_item is used as the write cursor of each child
   child_part => part%first_child_part
   do while (associated(child_part))
      child_part%last_item1 => null()
      child_part%last_item2 => null()
      child_part => child_part%next_sibling_part
   end do

   ! Atoms of molecule 1
   item => part%first_item1
   do while (associated(item))
      signature = hna_signature(adjcs1(item%idx), itemdir1)
      child_part => find_child_part(part, signature)
      if (DEBUG_TESTS) then
         if (.not. associated(child_part)) then
            call print_part_signature(signature)
            error stop 'Child part not found'
         end if
      end if
      if (.not. associated(child_part%last_item1)) then
         child_part%last_item1 => child_part%first_item1
      else
         child_part%last_item1 => child_part%last_item1%next_item
      end if
      child_part%last_item1%idx = item%idx
      link%itemdir1(item%idx)%ptr => child_part
      item => item%next_item
   end do

   ! Atoms of molecule 2
   item => part%first_item2
   do while (associated(item))
      signature = hna_signature(adjcs2(item%idx), itemdir2)
      child_part => find_child_part(part, signature)
      if (DEBUG_TESTS) then
         if (.not. associated(child_part)) then
            call print_part_signature(signature)
            error stop 'Child part not found'
         end if
      end if
      if (.not. associated(child_part%last_item2)) then
         child_part%last_item2 => child_part%first_item2
      else
         child_part%last_item2 => child_part%last_item2%next_item
      end if
      child_part%last_item2%idx = item%idx
      link%itemdir2(item%idx)%ptr => child_part
      item => item%next_item
   end do
end subroutine

subroutine update_hna_partition(adjcs1, adjcs2, link)
! Redistribute the atoms of all parts of link into the next link
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(partref_node_t), pointer :: partref

   partref => link%first_partref
   do while (associated(partref))
      call update_hna_part(adjcs1, adjcs2, link%itemdir1, link%itemdir2, partref%part, &
            link%next_link)
      partref => partref%nextref
   end do
end subroutine

subroutine assign_branch_atoms(adjcs1, adjcs2, branch)
! Individualize a random pair of the split part of branch and replay the
! refinement recorded in its links
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chaintree_node_t), pointer, intent(inout) :: branch
   ! Local variables
   type(chain_node_t), pointer :: link
   integer(ik) :: link_idx, rand_idx1, rand_idx2

   rand_idx1 = random_uniform_integer(1, branch%split_part%n_items1)
   rand_idx2 = random_uniform_integer(1, branch%split_part%n_items2)
   call split_part_indexed(branch%split_part, branch%first_link, rand_idx1, rand_idx2)

   link_idx = 1
   link => branch%first_link
   do while (associated(link))
      call update_hna_partition(adjcs1, adjcs2, link)
      link_idx = link_idx + 1
      link => link%next_link
   end do
end subroutine

subroutine split_part_indexed(part, link, index1, index2)
! Linked-list counterpart of assignment_conformer::assign_pair_to_children:
! the atoms at positions index1 and index2 of part go to its first child
! and the rest to its second child, overwriting their item nodes
   type(partition_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   integer(ik), intent(in) :: index1, index2
   ! Local variables
   type(partition_node_t), pointer :: child_part1, child_part2
   type(item_node_t), pointer :: item1, item2
   integer(ik) :: current_index

   child_part1 => part%first_child_part
   child_part2 => part%first_child_part%next_sibling_part

   ! Chosen atom of molecule 1 to the first child
   item1 => part%first_item1
   current_index = 1
   do while (current_index < index1)
      item1 => item1%next_item
      current_index = current_index + 1
   end do
   child_part1%first_item1%idx = item1%idx
   link%itemdir1(item1%idx)%ptr => child_part1

   ! Chosen atom of molecule 2 to the first child
   item2 => part%first_item2
   current_index = 1
   do while (current_index < index2)
      item2 => item2%next_item
      current_index = current_index + 1
   end do
   child_part1%first_item2%idx = item2%idx
   link%itemdir2(item2%idx)%ptr => child_part1

   ! Remaining atoms of molecule 1 to the second child
   item1 => part%first_item1
   current_index = 1
   child_part2%last_item1 => null()

   do while (associated(item1))
      if (current_index /= index1) then
         if (.not. associated(child_part2%last_item1)) then
            child_part2%last_item1 => child_part2%first_item1
         else
            child_part2%last_item1 => child_part2%last_item1%next_item
         end if
         child_part2%last_item1%idx = item1%idx
         link%itemdir1(item1%idx)%ptr => child_part2
      end if

      item1 => item1%next_item
      current_index = current_index + 1
   end do

   ! Remaining atoms of molecule 2 to the second child
   item2 => part%first_item2
   current_index = 1
   child_part2%last_item2 => null()

   do while (associated(item2))
      if (current_index /= index2) then
         if (.not. associated(child_part2%last_item2)) then
            child_part2%last_item2 => child_part2%first_item2
         else
            child_part2%last_item2 => child_part2%last_item2%next_item
         end if
         child_part2%last_item2%idx = item2%idx
         link%itemdir2(item2%idx)%ptr => child_part2
      end if

      item2 => item2%next_item
      current_index = current_index + 1
   end do
end subroutine

recursive subroutine distribute_items(adjcs1, adjcs2, branch)
! Assign a random pair at every node of the assignment tree, depth first.
! Not used by the searches, which work on the array representation.
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chaintree_node_t), pointer, intent(inout) :: branch
   type(chaintree_node_t), pointer :: child_branch
   type(partition_node_t), pointer :: child_part

   child_branch => branch%first_child_chain
   do while (associated(child_branch))
      child_part => child_branch%split_part%first_child_part
      call assign_branch_atoms(adjcs1, adjcs2, child_branch)
      call distribute_items(adjcs1, adjcs2, child_branch)
      child_branch => child_branch%next_sibling_chain
   end do
end subroutine

function would_part_split(adjcs1, adjcs2, itemdir1, itemdir2, part) result(would_split)
! Whether the atoms of part have more than one signature under itemdir
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir1, itemdir2
   type(partition_node_t), pointer, intent(inout) :: part
   ! Local variables
   logical(lk) :: would_split
   type(item_node_t), pointer :: item1, item2
   integer(ik), dimension(:), allocatable :: signature, reference

   would_split = .FALSE.
   item1 => part%first_item1
   item2 => part%first_item2

   ! Reference signature from the first atom
   if (associated(item1)) then
      reference = hna_signature(adjcs1(item1%idx), itemdir1)
      item1 => item1%next_item
   else if (associated(item2)) then
      reference = hna_signature(adjcs2(item2%idx), itemdir2)
      item2 => item2%next_item
   else
      return  ! No items to process
   end if

   ! Compare the remaining atoms of both molecules
   do while (associated(item1))
      signature = hna_signature(adjcs1(item1%idx), itemdir1)
      if (.not. (signature .equiv. reference)) then
         would_split = .TRUE.
         return
      end if
      item1 => item1%next_item
   end do

   do while (associated(item2))
      signature = hna_signature(adjcs2(item2%idx), itemdir2)
      if (.not. (signature .equiv. reference)) then
         would_split = .TRUE.
         return
      end if
      item2 => item2%next_item
   end do
end function

subroutine refine_branched_hna_partition(adjcs1, adjcs2, hna_chain, branch, branch_parts, n_splits)
! Refinement step during assignment tree construction (Pseudocode S2, lines
! 10-12). Unlike refine_hna_partition, only parts that actually split are
! subdivided; the others are carried over to the new level link without
! entering its vertex directory. Each part that splits is recorded in the
! current link of the node `branch`, so the search can replay the step, and
! its children are added to branch_parts as candidates for later
! individualization. n_splits is the number of parts that split.
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chaintree_node_t), pointer, intent(inout) :: hna_chain
   type(chaintree_node_t), pointer, intent(inout) :: branch
   type(chain_node_t), pointer, intent(inout) :: branch_parts
   integer(ik), intent(out) :: n_splits
   ! Local variables
   type(chain_node_t), pointer :: level_link, next_level_link, branch_link
   type(partref_node_t), pointer :: partref
   type(partition_node_t), pointer :: child_part
   logical(lk), dimension(:), allocatable :: will_split
   logical(lk) :: any_splits
   integer(ik) :: i

   n_splits = 0
   any_splits = .FALSE.

   level_link => hna_chain%last_link
   allocate(will_split(level_link%n_parts))

   ! Find the parts that split
   partref => level_link%first_partref
   do i = 1, level_link%n_parts
      will_split(i) = would_part_split(adjcs1, adjcs2, level_link%itemdir1, level_link%itemdir2, &
            partref%part)
      if (will_split(i)) any_splits = .TRUE.
      partref => partref%nextref
   end do

   ! New links are only created when some part splits
   if (any_splits) then
      next_level_link => new_chain_link(hna_chain)

      partref => level_link%first_partref
      do i = 1, level_link%n_parts
         if (will_split(i)) then
            call refine_hna_part(adjcs1, adjcs2, level_link%itemdir1, level_link%itemdir2, &
                  partref%part, next_level_link)
            call link_part(branch%last_link, partref%part)
            child_part => partref%part%first_child_part
            do while (associated(child_part))
               call add_branch_part(branch_parts, child_part)
               child_part => child_part%next_sibling_part
            end do
            n_splits = n_splits + 1
         else
            call link_part(next_level_link, partref%part)
         end if

         partref => partref%nextref
      end do
      branch_link => new_chain_link(branch)
   end if

   deallocate(will_split)
end subroutine

subroutine split_part_first(part, link)
! Individualize the first pair of part (Pseudocode S2, lines 8-9): its first
! atoms go to a new first child and the rest to a second child, both added
! to link and to its vertex directory
   type(partition_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(partition_node_t), pointer :: child_part
   type(item_node_t), pointer :: item1, item2

   ! Individualized pair
   child_part => new_child_part(part)
   call link_part(link, child_part)
   call add_new_item1(child_part, part%first_item1%idx)
   call add_new_item2(child_part, part%first_item2%idx)
   link%itemdir1(part%first_item1%idx)%ptr => child_part
   link%itemdir2(part%first_item2%idx)%ptr => child_part

   ! Remaining atoms
   child_part => new_child_part(part)
   call link_part(link, child_part)
   item1 => part%first_item1%next_item
   item2 => part%first_item2%next_item
   do while (associated(item1))
      call add_new_item1(child_part, item1%idx)
      call add_new_item2(child_part, item2%idx)
      link%itemdir1(item1%idx)%ptr => child_part
      link%itemdir2(item2%idx)%ptr => child_part
      item1 => item1%next_item
      item2 => item2%next_item
   end do
end subroutine

recursive subroutine split_dependent_parts(adjcs1, adjcs2, hna_chain, branch, branch_parts, &
      branching_part, part_to_split, error_code)
! Create a child node of branch for part_to_split, individualize its first
! pair and refine to self-consistency (Pseudocode S2, lines 6-14). While a
! part with several pairs descending from branching_part remains, it is
! split in turn under a new nested node, so the degeneracy of branching_part
! is resolved along a single path of nodes. On return branch points to the
! deepest node of that path, and branch_parts holds the children of all
! parts that split on the way.
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chaintree_node_t), pointer, intent(inout) :: hna_chain
   type(chaintree_node_t), pointer, intent(inout) :: branch
   type(chain_node_t), pointer, intent(inout) :: branch_parts
   type(partition_node_t), pointer, intent(in) :: branching_part
   type(partition_node_t), pointer, intent(inout) :: part_to_split
   integer(ik), intent(out) :: error_code
   ! Local variables
   type(partref_node_t), pointer :: partref
   type(partition_node_t), pointer :: next_part_to_split
   type(chain_node_t), pointer :: level_link, next_level_link
   type(chain_node_t), pointer :: first_branch_link
   type(partition_node_t), pointer :: child_part
   integer(ik) :: n_splits

   error_code = MOLALIGN_SUCCESS

   level_link => hna_chain%last_link
   next_level_link => new_chain_link(hna_chain)

   ! New assignment tree node for this split part
   branch => new_child_chain(branch, part_to_split)
   first_branch_link => new_chain_link(branch)

   ! Carry the other parts over to the new level
   partref => level_link%first_partref
   do while (associated(partref))
      if (.not. associated(partref%part, part_to_split)) then
         call link_part(next_level_link, partref%part)
      end if
      partref => partref%nextref
   end do

   call split_part_first(part_to_split, next_level_link)

   child_part => part_to_split%first_child_part
   do while (associated(child_part))
      call add_branch_part(branch_parts, child_part)
      child_part => child_part%next_sibling_part
   end do

   ! Refine to self-consistency
   do
      call refine_branched_hna_partition(adjcs1, adjcs2, hna_chain, branch, branch_parts, n_splits)
      if (n_splits == 0) exit
   end do

   if (is_partition_uneven(hna_chain%last_link)) then
      error_code = MOLALIGN_ERROR_NOT_CONFORMERS
      return
   end if

   ! Next degenerate part descending from branching_part
   next_part_to_split => null()
   partref => hna_chain%last_link%first_partref
   do while (associated(partref) .and. .not. associated(next_part_to_split))
      if (partref%part%n_items1 >= 2) then
         if (isdescendant(partref%part, branching_part)) then
            next_part_to_split => partref%part
         end if
      end if
      partref => partref%nextref
   end do

   if (associated(next_part_to_split)) then
      call split_dependent_parts(adjcs1, adjcs2, hna_chain, branch, branch_parts, branching_part, &
            next_part_to_split, error_code)
   end if
end subroutine

recursive subroutine split_independent_parts(adjcs1, adjcs2, hna_chain, branch, branch_parts, error_code)
! Build the subtrees of the independent parts in branch_parts (outer loop
! and recursion of Pseudocode S2). Each part that is still a leaf (not
! subdivided by earlier individualizations) starts a path of nodes below
! branch; the parts split along that path are then processed recursively
! below its deepest node. branch_parts is sorted by part size, so the
! smallest parts are individualized first.
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chaintree_node_t), pointer, intent(inout) :: hna_chain, branch
   type(chain_node_t), pointer, intent(in) :: branch_parts
   integer(ik), intent(out) :: error_code
   ! Local variables
   type(chain_node_t), pointer :: new_branch_parts
   type(chaintree_node_t), pointer :: new_branch
   type(partref_node_t), pointer :: partref

   error_code = MOLALIGN_SUCCESS

   partref => branch_parts%first_partref
   do while (associated(partref))
      if (partref%part%n_children == 0) then
         new_branch => branch
         new_branch_parts => new_bare_link()
         call split_dependent_parts(adjcs1, adjcs2, hna_chain, new_branch, new_branch_parts, &
               partref%part, partref%part, error_code)
         if (error_code /= MOLALIGN_SUCCESS) return
         call split_independent_parts(adjcs1, adjcs2, hna_chain, new_branch, new_branch_parts, error_code)
         if (error_code /= MOLALIGN_SUCCESS) return
      end if
      partref => partref%nextref
   end do
end subroutine

subroutine build_assignment_tree( adjcs1, adjcs2, hna_link, cache_arrays, error_code)
! Build the assignment tree from the self-consistent HNA partition hna_link
! (Pseudocode S2) and store it, together with the partition tree and the
! adjacency lists, in the fixed-size arrays used by the searches. The
! partition is rebuilt under a fresh root with one child per part of
! hna_link; all parts with several pairs are the initial candidates for
! individualization.
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chain_node_t), pointer, intent(in) :: hna_link
   type(array_trees_t), intent(out) :: cache_arrays
   integer(ik), intent(out) :: error_code
   ! Local variables
   type(partition_node_t), pointer :: partition_tree
   type(chaintree_node_t), pointer :: assignment_tree
   type(chaintree_node_t), pointer :: hna_chain
   type(chain_node_t), pointer :: branch_parts
   type(partition_node_t), pointer :: child_part
   type(chain_node_t), pointer :: first_link
   type(partref_node_t), pointer :: partref

   error_code = MOLALIGN_SUCCESS

   partition_tree => new_root_part()
   branch_parts => new_bare_link()
   assignment_tree => new_root_chain( size(hna_link%itemdir1), size(hna_link%itemdir2))
   hna_chain => new_root_chain( size(hna_link%itemdir1), size(hna_link%itemdir2))
   first_link => new_chain_link( hna_chain)

   partref => hna_link%first_partref
   do while (associated( partref))
      child_part => new_child_part( partition_tree)
      call link_part( first_link, child_part)
      call copy_part_items( partref%part, child_part)
      call add_branch_part( branch_parts, child_part)
      partref => partref%nextref
   end do

   call split_independent_parts( adjcs1, adjcs2, hna_chain, assignment_tree, branch_parts, error_code)
   if (error_code /= MOLALIGN_SUCCESS) return

   ! Array representation for the searches
   call cache_partition_tree( partition_tree, cache_arrays)
   call cache_assignment_tree( assignment_tree, cache_arrays)
   call cache_adjacency_lists( adjcs1, adjcs2, cache_arrays)
end subroutine

end module
