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
use parameters
use types_basic
use random
use chemistry
use adjacency
use molecule
use types_linked
use types_indexed
use options
implicit none
private
public refine_hna_part
public refine_hna_partition
public compute_scna_partition
public build_assignment_tree

contains

subroutine refine_hna_part(adjcs1, adjcs2, itemdir1, itemdir2, part, link)
! Create children for different signatures - caller decides what to do with them
! Note: part is always a leaf part with no existing children
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir1, itemdir2
   type(partition_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(item_node_t), pointer :: item
   type(partition_node_t), pointer :: child_part
   type(part_nodeptr_t), dimension(:), allocatable :: signature

   ! Process first molecule items - create children for each unique signature
   item => part%first_item1
   do while (associated(item))
      signature = itemdir1(adjcs1(item%idx)%list(:adjcs1(item%idx)%cn))
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

   ! Process second molecule items - create children for each unique signature
   item => part%first_item2
   do while (associated(item))
      signature = itemdir2(adjcs2(item%idx)%list(:adjcs2(item%idx)%cn))
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

subroutine refine_hna_partition(adjcs1, adjcs2, hna_chain, num_splits)
! Compute next level HNAs - always keeps all children (original behavior)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chaintree_node_t), pointer, intent(inout) :: hna_chain
   integer, intent(out) :: num_splits
   ! Local variables
   type(chain_node_t), pointer :: link, new_link
   type(partref_node_t), pointer :: partref

   num_splits = 0

   ! Save the last link before creating a new one
   link => hna_chain%last_link
   new_link => new_chain_link(hna_chain)

   ! Process all parts in the current partition
   partref => link%first_partref
   do while (associated(partref))
      ! Create children based on signatures
      call refine_hna_part(adjcs1, adjcs2, link%itemdir1, link%itemdir2, partref%part, new_link)
      ! Count splits (children beyond the original part)
      num_splits = num_splits + partref%part%num_children - 1
      partref => partref%nextref
   end do
end subroutine

subroutine compute_scna_partition(adjcs1, adjcs2, atomtypes, hna_chain)
! Iteratively compute HNAs until convergence
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(partition_t), intent(in) :: atomtypes
   type(chaintree_node_t), pointer, intent(out) :: hna_chain
   ! Local variables
   integer :: num_splits

!   hna_chain => collect_atomtypes_linked( adjcs1, adjcs2)
   hna_chain => chain_from_partition( atomtypes)

   do
      ! Call refine_hna_partition and get the number of splits
      call refine_hna_partition(adjcs1, adjcs2, hna_chain, num_splits)
      ! Exit loop if no splits occurred
      if (num_splits == 0) exit
   end do

   ! Verify that molecules are conformers
   if (is_partition_uneven(hna_chain%last_link)) then
      write(stderr, '(A)') 'STOP Molecules are not conformers!'
!      call print_partition_chain(hna_chain)
      stop
   end if
end subroutine

subroutine update_hna_part(adjcs1, adjcs2, itemdir1, itemdir2, part, link)
! Update item values of existing item nodes instead of adding new item nodes
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir1, itemdir2
   type(partition_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(item_node_t), pointer :: item
   type(part_nodeptr_t), dimension(:), allocatable :: signature
   type(partition_node_t), pointer :: child_part

   ! Reset last item pointers for all children
   child_part => part%first_child_part
   do while (associated(child_part))
      child_part%last_item1 => null()
      child_part%last_item2 => null()
      child_part => child_part%next_sibling_part
   end do

   ! Process first molecule items - update existing item nodes
   item => part%first_item1
   do while (associated(item))
      signature = itemdir1(adjcs1(item%idx)%list(:adjcs1(item%idx)%cn))
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

   ! Process second molecule items - update existing item nodes
   item => part%first_item2
   do while (associated(item))
      signature = itemdir2(adjcs2(item%idx)%list(:adjcs2(item%idx)%cn))
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
! Compute next level HNAs - always keeps all children (original behavior)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(partref_node_t), pointer :: partref

   ! Process all parts in the current partition
   partref => link%first_partref
   do while (associated(partref))
!      write (stderr,'(A,1X,A)') 'Part', address(partref%part)
      ! Distribute items based on signatures
      call update_hna_part(adjcs1, adjcs2, link%itemdir1, link%itemdir2, partref%part, &
            link%next_link)
      partref => partref%nextref
   end do
end subroutine

subroutine assign_branch_atoms(adjcs1, adjcs2, branch)
! Iteratively compute HNAs until convergence
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chaintree_node_t), pointer, intent(inout) :: branch
   ! Local variables
   type(chain_node_t), pointer :: link
   integer :: link_idx, rand_idx1, rand_idx2

   rand_idx1 = random_uniform_integer(1, branch%split_part%num_items1)
   rand_idx2 = random_uniform_integer(1, branch%split_part%num_items2)
   call split_part_indexed(branch%split_part, branch%first_link, rand_idx1, rand_idx2)
!   call split_part_indexed(branch%split_part, branch%first_link, 1, 1)

   link_idx = 1
   link => branch%first_link
   do while (associated(link))
      ! Recompute next level HNAs
!      write (stderr,'(A,1X,I0)') 'Link', link_idx
!      call print_link_itemdir(link)
      call update_hna_partition(adjcs1, adjcs2, link)
      link_idx = link_idx + 1
      link => link%next_link
   end do
end subroutine

subroutine split_part_indexed(part, link, index1, index2)
   type(partition_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   integer, intent(in) :: index1, index2
   ! Local variables
   type(partition_node_t), pointer :: child_part1, child_part2
   type(item_node_t), pointer :: item1, item2
   integer :: current_index

   ! Get pointers to first and second child parts
   child_part1 => part%first_child_part
   child_part2 => part%first_child_part%next_sibling_part

   ! Find and assign the item at index1 from items1 to first child
   item1 => part%first_item1
   current_index = 1
   do while (current_index < index1)
      item1 => item1%next_item
      current_index = current_index + 1
   end do
   child_part1%first_item1%idx = item1%idx
   link%itemdir1(item1%idx)%ptr => child_part1

   ! Find and assign the item at index2 from items2 to first child
   item2 => part%first_item2
   current_index = 1
   do while (current_index < index2)
      item2 => item2%next_item
      current_index = current_index + 1
   end do
   child_part1%first_item2%idx = item2%idx
   link%itemdir2(item2%idx)%ptr => child_part1

   ! Now assign all other items from items1 to second child
   item1 => part%first_item1
   current_index = 1
   child_part2%last_item1 => null()

   do while (associated(item1))
      if (current_index /= index1) then
         ! Move to next position in second child
         if (.not. associated(child_part2%last_item1)) then
            child_part2%last_item1 => child_part2%first_item1
         else
            child_part2%last_item1 => child_part2%last_item1%next_item
         end if
         ! This is not the selected item, assign to second child
         child_part2%last_item1%idx = item1%idx
         link%itemdir1(item1%idx)%ptr => child_part2
      end if

      ! Move to next item in parent
      item1 => item1%next_item
      current_index = current_index + 1
   end do

   ! Now assign all other items from items2 to second child
   item2 => part%first_item2
   current_index = 1
   child_part2%last_item2 => null()

   do while (associated(item2))
      if (current_index /= index2) then
         ! Move to next position in second child
         if (.not. associated(child_part2%last_item2)) then
            child_part2%last_item2 => child_part2%first_item2
         else
            child_part2%last_item2 => child_part2%last_item2%next_item
         end if
         ! This is not the selected item, assign to second child
         child_part2%last_item2%idx = item2%idx
         link%itemdir2(item2%idx)%ptr => child_part2
      end if

      ! Move to next item in parent
      item2 => item2%next_item
      current_index = current_index + 1
   end do
end subroutine

recursive subroutine distribute_items(adjcs1, adjcs2, branch)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chaintree_node_t), pointer, intent(inout) :: branch
   type(chaintree_node_t), pointer :: child_branch
   type(partition_node_t), pointer :: child_part

!   call random_init(.TRUE., .TRUE.)

   ! Process all children of this branch
   child_branch => branch%first_child_chain
   do while (associated(child_branch))
      ! Process this child branch's split_part
!      write (stderr,*)
!      write (stderr,'(A,1X,A)') 'Split Part', address(child_branch%split_part)
      ! Add target part children to links
      child_part => child_branch%split_part%first_child_part
      call assign_branch_atoms(adjcs1, adjcs2, child_branch)

      ! Recursively process this child's descendants (depth-first)
      call distribute_items(adjcs1, adjcs2, child_branch)

      ! Move to next sibling
      child_branch => child_branch%next_sibling_chain
   end do
end subroutine

function would_part_split(adjcs1, adjcs2, itemdir1, itemdir2, part) result(would_split)
! Check if a part would split by comparing signatures
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir1, itemdir2
   type(partition_node_t), pointer, intent(inout) :: part
   ! Local variables
   logical :: would_split
   type(item_node_t), pointer :: item1, item2
   type(part_nodeptr_t), dimension(:), allocatable :: signature, reference

   would_split = .FALSE.
   item1 => part%first_item1
   item2 => part%first_item2

   ! Set first signature from first available item
   if (associated(item1)) then
      reference = itemdir1(adjcs1(item1%idx)%list(:adjcs1(item1%idx)%cn))
      item1 => item1%next_item
   else if (associated(item2)) then
      reference = itemdir2(adjcs2(item2%idx)%list(:adjcs2(item2%idx)%cn))
      item2 => item2%next_item
   else
      return  ! No items to process
   end if

   ! Check remaining items in first molecule
   do while (associated(item1))
      signature = itemdir1(adjcs1(item1%idx)%list(:adjcs1(item1%idx)%cn))
      if (.not. (signature .equiv. reference)) then
         would_split = .TRUE.
         return
      end if
      item1 => item1%next_item
   end do

   ! Check remaining items in second molecule
   do while (associated(item2))
      signature = itemdir2(adjcs2(item2%idx)%list(:adjcs2(item2%idx)%cn))
      if (.not. (signature .equiv. reference)) then
         would_split = .TRUE.
         return
      end if
      item2 => item2%next_item
   end do
end function

subroutine refine_branched_hna_partition(adjcs1, adjcs2, hna_chain, branch, branch_parts, num_splits)
! Compute next level HNAs - only keeps children if real split occurred
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chaintree_node_t), pointer, intent(inout) :: hna_chain
   type(chaintree_node_t), pointer, intent(inout) :: branch
   type(chain_node_t), pointer, intent(inout) :: branch_parts
   integer, intent(out) :: num_splits
   ! Local variables
   type(chain_node_t), pointer :: level_link, next_level_link, branch_link
   type(partref_node_t), pointer :: partref
   type(partition_node_t), pointer :: child_part
   logical, dimension(:), allocatable :: will_split
   logical :: any_splits
   integer :: i

   num_splits = 0
   any_splits = .FALSE.

   ! Get the current level link
   level_link => hna_chain%last_link

   ! Allocate array to cache split results
   allocate(will_split(level_link%num_parts))

   ! Single pass: check which parts would split and cache results
   partref => level_link%first_partref
   do i = 1, level_link%num_parts
      will_split(i) = would_part_split(adjcs1, adjcs2, level_link%itemdir1, level_link%itemdir2, &
            partref%part)
      if (will_split(i)) any_splits = .TRUE.
      partref => partref%nextref
   end do

   ! Only create new links if splits will occur
   if (any_splits) then
      ! Create new level link for hna_chain
      next_level_link => new_chain_link(hna_chain)

      ! Process all parts using cached split results
      partref => level_link%first_partref
      do i = 1, level_link%num_parts
         if (will_split(i)) then
            ! Create children based on signatures
            call refine_hna_part(adjcs1, adjcs2, level_link%itemdir1, level_link%itemdir2, &
                  partref%part, next_level_link)
            ! Link part to branch link
            call link_part(branch%last_link, partref%part)
            ! Add part children to branch part list
            child_part => partref%part%first_child_part
            do while (associated(child_part))
               call add_branch_part(branch_parts, child_part)
               child_part => child_part%next_sibling_part
            end do
            num_splits = num_splits + 1
         else
            ! Link original part
            call link_part(next_level_link, partref%part)
         end if

         partref => partref%nextref
      end do
      branch_link => new_chain_link(branch)
   end if

   ! Clean up
   deallocate(will_split)
end subroutine

subroutine split_part_first(part, link)
   type(partition_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(partition_node_t), pointer :: child_part
   type(item_node_t), pointer :: item1, item2

   ! Create first child and add first item from each molecule
   child_part => new_child_part(part)
   call link_part(link, child_part)
   call add_new_item1(child_part, part%first_item1%idx)
   call add_new_item2(child_part, part%first_item2%idx)
   link%itemdir1(part%first_item1%idx)%ptr => child_part
   link%itemdir2(part%first_item2%idx)%ptr => child_part

   ! Create second child and add remaining items
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

! Modified split_dependent_parts incorporating split_single_part functionality
recursive subroutine split_dependent_parts(adjcs1, adjcs2, hna_chain, branch, branch_parts, &
      branching_part, part_to_split)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chaintree_node_t), pointer, intent(inout) :: hna_chain
   type(chaintree_node_t), pointer, intent(inout) :: branch
   type(chain_node_t), pointer, intent(inout) :: branch_parts
   type(partition_node_t), pointer, intent(in) :: branching_part
   type(partition_node_t), pointer, intent(inout) :: part_to_split
   ! Local variables
   type(partref_node_t), pointer :: partref
   type(partition_node_t), pointer :: next_part_to_split
   type(chain_node_t), pointer :: level_link, next_level_link
   type(chain_node_t), pointer :: first_branch_link
   type(partition_node_t), pointer :: child_part
   integer :: num_splits

   ! Incorporate split_single_part logic
   ! Save the last link before creating a new one
   level_link => hna_chain%last_link
   next_level_link => new_chain_link(hna_chain)

   ! Create a new child branch for this splitting part
   branch => new_child_chain(branch, part_to_split)
   first_branch_link => new_chain_link(branch)

   ! Add non splitting parts to new link
   partref => level_link%first_partref
   do while (associated(partref))
      if (.not. associated(partref%part, part_to_split)) then
         call link_part(next_level_link, partref%part)
      end if
      partref => partref%nextref
   end do

   ! Split splitting part
   call split_part_first(part_to_split, next_level_link)

   ! Add part children to branch part list
   child_part => part_to_split%first_child_part
   do while (associated(child_part))
      call add_branch_part(branch_parts, child_part)
      child_part => child_part%next_sibling_part
   end do

   ! Compute self-consistent HNAs
   do
      ! Refine HNA partition and get the number of splits
      call refine_branched_hna_partition(adjcs1, adjcs2, hna_chain, branch, branch_parts, num_splits)
      ! Exit loop if no splits occurred
      if (num_splits == 0) exit
   end do

   ! Verify that molecules are conformers
   if (is_partition_uneven(hna_chain%last_link)) then
      write(stderr, '(A)') 'STOP Molecules are not conformers!'
!      call print_partition_chain(hna_chain)
      stop
   end if

   ! Find a degenerate descendant part to split
   next_part_to_split => null()
   partref => hna_chain%last_link%first_partref
   do while (associated(partref) .and. .not. associated(next_part_to_split))
      if (partref%part%num_items1 >= 2) then
         if (isdescendant(partref%part, branching_part)) then
            next_part_to_split => partref%part
         end if
      end if
      partref => partref%nextref
   end do

   ! Perform split if target found
   if (associated(next_part_to_split)) then
      ! Call itself again to split the next degenerate descendant part
      call split_dependent_parts(adjcs1, adjcs2, hna_chain, branch, branch_parts, branching_part, &
            next_part_to_split)
   end if
end subroutine

! Updated split_independent_parts to use the merged function signature
recursive subroutine split_independent_parts(adjcs1, adjcs2, hna_chain, branch, branch_parts)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chaintree_node_t), pointer, intent(inout) :: hna_chain, branch
   type(chain_node_t), pointer, intent(in) :: branch_parts
   ! Local variables
   type(chain_node_t), pointer :: new_branch_parts
   type(chaintree_node_t), pointer :: new_branch
   type(partref_node_t), pointer :: partref

   ! Process each part in branch_parts
   partref => branch_parts%first_partref
   do while (associated(partref))
      if (partref%part%num_children == 0) then
         ! Start leaf chain from current branch chain
         new_branch => branch
         ! Create a new part registry for this branch part
         new_branch_parts => new_bare_link()
         ! Split the target part and continue splitting descendants until convergence
         call split_dependent_parts(adjcs1, adjcs2, hna_chain, new_branch, new_branch_parts, &
               partref%part, partref%part)
         ! Recursively process the resulting branch parts
         call split_independent_parts(adjcs1, adjcs2, hna_chain, new_branch, new_branch_parts)
      end if
      partref => partref%nextref
   end do
end subroutine

subroutine build_assignment_tree( adjcs1, adjcs2, hna_link, cache_arrays)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chain_node_t), pointer, intent(in) :: hna_link
   type(array_trees_t), intent(out) :: cache_arrays
   ! Local variables
   type(partition_node_t), pointer :: partition_tree
   type(chaintree_node_t), pointer :: assignment_tree
   type(chaintree_node_t), pointer :: hna_chain
   type(chain_node_t), pointer :: branch_parts
   type(partition_node_t), pointer :: child_part
   type(chain_node_t), pointer :: first_link
   type(partref_node_t), pointer :: partref

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

   call split_independent_parts( adjcs1, adjcs2, hna_chain, assignment_tree, branch_parts)
!   call distribute_items( adjcs1, adjcs2, assignment_tree)

   ! Cache required data for DFS in fixed size arrays
   call cache_partition_tree( partition_tree, cache_arrays)
   call cache_assignment_tree( assignment_tree, cache_arrays)
   call cache_adjacency_lists( adjcs1, adjcs2, cache_arrays)
end subroutine

end module
