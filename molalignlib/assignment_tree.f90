module assignment_tree
use parameters
use random
use molecule
use lcrs_trees
use lcrs_arrays
use hna
implicit none
private
public build_assignment_tree

contains

subroutine update_hna_part(adjcs1, adjcs2, itemdir1, itemdir2, part, link)
! Update item values of existing item nodes instead of adding new item nodes
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir1, itemdir2
   type(partree_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(item_node_t), pointer :: item
   type(part_nodeptr_t), dimension(:), allocatable :: signature
   type(partree_node_t), pointer :: child_part

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
      signature = itemdir1(adjcs1(item%idx)%adjlist)
      child_part => find_child_part(part, signature)
      if (.not. associated(child_part)) then
         call print_part_signature(signature)
         error stop 'part not found'
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
      signature = itemdir2(adjcs2(item%idx)%adjlist)
      child_part => find_child_part(part, signature)
      if (.not. associated(child_part)) then
         call print_part_signature(signature)
         error stop 'part not found'
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
! Compute next level HNA types - always keeps all children (original behavior)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(partref_node_t), pointer :: partref

   ! Process all parts in the current partition
   partref => link%first_partref
   do while (associated(partref))
!      write (stderr,'(A,1X,A)') 'Part', address(partref%part)
      ! Distribute items based on signatures
      call update_hna_part(adjcs1, adjcs2, link%itemdir1, link%itemdir2, partref%part, link%next_link)
      partref => partref%nextref
   end do
end subroutine

subroutine assign_branch_atoms(adjcs1, adjcs2, branch)
! Iteratively compute HNA types until convergence
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(assigntree_node_t), pointer, intent(inout) :: branch
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
   type(partree_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   integer, intent(in) :: index1, index2
   ! Local variables
   type(partree_node_t), pointer :: child_part1, child_part2
   type(item_node_t), pointer :: item1, item2
   integer :: current_index

   ! Validate indices
   if (index1 < 1 .or. index1 > part%num_items1) error stop "index1 out of range"
   if (index2 < 1 .or. index2 > part%num_items2) error stop "index2 out of range"

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
   type(assigntree_node_t), pointer, intent(inout) :: branch
   type(assigntree_node_t), pointer :: child_branch
   type(partree_node_t), pointer :: child_part

!   call random_init(.true., .true.)

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
   type(partree_node_t), pointer, intent(inout) :: part
   ! Local variables
   logical :: would_split
   type(item_node_t), pointer :: item1, item2
   type(part_nodeptr_t), dimension(:), allocatable :: signature, first_signature

   would_split = .false.
   item1 => part%first_item1
   item2 => part%first_item2

   ! Set first signature from first available item
   if (associated(item1)) then
      first_signature = itemdir1(adjcs1(item1%idx)%adjlist)
      item1 => item1%next_item
   else if (associated(item2)) then
      first_signature = itemdir2(adjcs2(item2%idx)%adjlist)
      item2 => item2%next_item
   else
      return  ! No items to process
   end if

   ! Check remaining items in first molecule
   do while (associated(item1))
      signature = itemdir1(adjcs1(item1%idx)%adjlist)
      if (.not. (signature .equiv. first_signature)) then
         would_split = .true.
         return
      end if
      item1 => item1%next_item
   end do

   ! Check remaining items in second molecule
   do while (associated(item2))
      signature = itemdir2(adjcs2(item2%idx)%adjlist)
      if (.not. (signature .equiv. first_signature)) then
         would_split = .true.
         return
      end if
      item2 => item2%next_item
   end do
end function

subroutine recompute_consistent_hna_partition(adjcs1, adjcs2, hnachain, branch, branch_parts, num_splits)
! Compute next level HNA types - only keeps children if real split occurred
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(assigntree_node_t), pointer, intent(inout) :: hnachain
   type(assigntree_node_t), pointer, intent(inout) :: branch
   type(chain_node_t), pointer, intent(inout) :: branch_parts
   integer, intent(out) :: num_splits
   ! Local variables
   type(chain_node_t), pointer :: level_link, next_level_link, branch_link
   type(partref_node_t), pointer :: partref
   type(partree_node_t), pointer :: child_part
   logical, dimension(:), allocatable :: will_split
   logical :: any_splits
   integer :: i

   num_splits = 0
   any_splits = .false.

   ! Get the current level link
   level_link => hnachain%last_link

   ! Allocate array to cache split results
   allocate(will_split(level_link%num_parts))

   ! Single pass: check which parts would split and cache results
   partref => level_link%first_partref
   do i = 1, level_link%num_parts
      will_split(i) = would_part_split(adjcs1, adjcs2, level_link%itemdir1, level_link%itemdir2, partref%part)
      if (will_split(i)) any_splits = .true.
      partref => partref%nextref
   end do

   ! Only create new links if splits will occur
   if (any_splits) then
      ! Create new level link for hnachain
      next_level_link => new_chain_link(hnachain)

      ! Process all parts using cached split results
      partref => level_link%first_partref
      do i = 1, level_link%num_parts
         if (will_split(i)) then
            ! Create children based on signatures
            call refine_hna_part(adjcs1, adjcs2, level_link%itemdir1, level_link%itemdir2, partref%part, next_level_link)
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
   type(partree_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(partree_node_t), pointer :: child_part
   type(item_node_t), pointer :: item1, item2

   ! Verify that both molecules are conformers
   if (part%num_items1 /= part%num_items2) then
      write(stderr, '(A)') "ERROR: Molecules are not conformers!"
      stop
   end if

   if (DO_DEBUG_TESTS) then
      ! Verify that the item chains are properly linked
      item1 => part%first_item1%next_item
      item2 => part%first_item2%next_item
      do while (associated(item1) .and. associated(item2))
         item1 => item1%next_item
         item2 => item2%next_item
      end do
      if (associated(item1) .neqv. associated(item2)) then
         error stop "FATAL: Item chains have different lengths."
      end if
   end if

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
recursive subroutine split_dependent_parts(adjcs1, adjcs2, hnachain, branch, branch_parts, branching_part, part_to_split)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(assigntree_node_t), pointer, intent(inout) :: hnachain
   type(assigntree_node_t), pointer, intent(inout) :: branch
   type(chain_node_t), pointer, intent(inout) :: branch_parts
   type(partree_node_t), pointer, intent(in) :: branching_part
   type(partree_node_t), pointer, intent(inout) :: part_to_split
   ! Local variables
   type(partref_node_t), pointer :: partref
   type(partree_node_t), pointer :: next_part_to_split
   type(chain_node_t), pointer :: level_link, next_level_link
   type(chain_node_t), pointer :: first_branch_link
   type(partree_node_t), pointer :: child_part
   integer :: num_splits

   ! Incorporate split_single_part logic
   ! Save the last link before creating a new one
   level_link => hnachain%last_link
   next_level_link => new_chain_link(hnachain)

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
      ! Call recompute_consistent_hna_partition and get the number of splits
      call recompute_consistent_hna_partition(adjcs1, adjcs2, hnachain, branch, branch_parts, num_splits)

      ! Exit loop if no splits occurred
      if (num_splits == 0) exit
   end do

   ! Find a degenerate descendant part to split
   next_part_to_split => null()
   partref => hnachain%last_link%first_partref
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
      call split_dependent_parts(adjcs1, adjcs2, hnachain, branch, branch_parts, branching_part, next_part_to_split)
   end if
end subroutine

! Updated split_independent_parts to use the merged function signature
recursive subroutine split_independent_parts(adjcs1, adjcs2, hnachain, branch, branch_parts)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(assigntree_node_t), pointer, intent(inout) :: hnachain, branch
   type(chain_node_t), pointer, intent(in) :: branch_parts
   ! Local variables
   type(chain_node_t), pointer :: new_branch_parts
   type(assigntree_node_t), pointer :: new_branch
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
         call split_dependent_parts(adjcs1, adjcs2, hnachain, new_branch, new_branch_parts, partref%part, partref%part)
         ! Recursively process the resulting branch parts
         call split_independent_parts(adjcs1, adjcs2, hnachain, new_branch, new_branch_parts)
      end if
      partref => partref%nextref
   end do
end subroutine

subroutine build_assignment_tree( adjcs1, adjcs2, hnalink, assign_arrays)
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   type(chain_node_t), pointer, intent(in) :: hnalink
   type(array_trees_t), intent(out) :: assign_arrays
   ! Local variables
   type(partree_node_t), pointer :: part_tree
   type(assigntree_node_t), pointer :: assign_tree
   type(assigntree_node_t), pointer :: hnachain
   type(chain_node_t), pointer :: branch_parts
   type(partree_node_t), pointer :: child_part
   type(chain_node_t), pointer :: first_link
   type(partref_node_t), pointer :: partref

   part_tree => new_root_part()
   branch_parts => new_bare_link()
   assign_tree => new_root_chain( size(hnalink%itemdir1), size(hnalink%itemdir2))
   hnachain => new_root_chain( size(hnalink%itemdir1), size(hnalink%itemdir2))
   first_link => new_chain_link( hnachain)

   partref => hnalink%first_partref
   do while (associated( partref))
      child_part => new_child_part( part_tree)
      call link_part( first_link, child_part)
      call copy_part_items( partref%part, child_part)
      call add_branch_part( branch_parts, child_part)
      partref => partref%nextref
   end do

   call split_independent_parts( adjcs1, adjcs2, hnachain, assign_tree, branch_parts)
!   call distribute_items( adjcs1, adjcs2, assign_tree)
   call convert_trees_to_arrays( adjcs1, adjcs2, part_tree, assign_tree, assign_arrays)
!   call validate_conversion(part_tree, assign_tree, assign_arrays)
end subroutine

end module
