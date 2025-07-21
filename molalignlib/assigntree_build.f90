module assigntree_build
use parameters
use molecule
use lcrs_tree
use lcrs_frame
use atom_mnas
implicit none
private
public build_assignment_tree

contains

function would_part_split(atoms1, atoms2, itemdir1, itemdir2, part) result(would_split)
! Check if a part would split by comparing signatures
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
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
      first_signature = itemdir1(atoms1(item1%value)%adjlist)
      item1 => item1%next_item
   else if (associated(item2)) then
      first_signature = itemdir2(atoms2(item2%value)%adjlist)
      item2 => item2%next_item
   else
      return  ! No items to process
   end if

   ! Check remaining items in first molecule
   do while (associated(item1))
      signature = itemdir1(atoms1(item1%value)%adjlist)
      if (.not. (signature .equiv. first_signature)) then
         would_split = .true.
         return
      end if
      item1 => item1%next_item
   end do

   ! Check remaining items in second molecule
   do while (associated(item2))
      signature = itemdir2(atoms2(item2%value)%adjlist)
      if (.not. (signature .equiv. first_signature)) then
         would_split = .true.
         return
      end if
      item2 => item2%next_item
   end do
end function

subroutine compute_mna_partition(atoms1, atoms2, mnachain, branch, branch_parts, num_splits)
! Compute next level MNA types - only keeps children if real split occurred
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(assigntree_node_t), pointer, intent(inout) :: mnachain
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
   level_link => mnachain%last_link

   ! Allocate array to cache split results
   allocate(will_split(level_link%num_parts))

   ! Single pass: check which parts would split and cache results
   partref => level_link%first_partref
   do i = 1, level_link%num_parts
      will_split(i) = would_part_split(atoms1, atoms2, level_link%itemdir1, level_link%itemdir2, partref%part)
      if (will_split(i)) any_splits = .true.
      partref => partref%nextref
   end do

   ! Only create new links if splits will occur
   if (any_splits) then
      ! Create new level link for mnachain
      next_level_link => new_chain_link(mnachain)

      ! Process all parts using cached split results
      partref => level_link%first_partref
      do i = 1, level_link%num_parts
         if (will_split(i)) then
            ! Create children based on signatures
            call split_part_mna(atoms1, atoms2, level_link%itemdir1, level_link%itemdir2, partref%part, next_level_link)
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

   ! Create first child and add first item from each molecule
   child_part => new_child_part(part)
   call link_part(link, child_part)
   call add_new_item1(child_part, part%first_item1%value)
   call add_new_item2(child_part, part%first_item2%value)
   link%itemdir1(part%first_item1%value)%ptr => child_part
   link%itemdir2(part%first_item2%value)%ptr => child_part

   ! Create second child and add remaining items
   child_part => new_child_part(part)
   call link_part(link, child_part)
   item1 => part%first_item1%next_item
   item2 => part%first_item2%next_item
   do while (associated(item1))
      call add_new_item1(child_part, item1%value)
      call add_new_item2(child_part, item2%value)
      link%itemdir1(item1%value)%ptr => child_part
      link%itemdir2(item2%value)%ptr => child_part
      item1 => item1%next_item
      item2 => item2%next_item
   end do
end subroutine

! Modified split_dependent_parts incorporating split_single_part functionality
recursive subroutine split_dependent_parts(atoms1, atoms2, mnachain, branch, branch_parts, branching_part, part_to_split)
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(assigntree_node_t), pointer, intent(inout) :: mnachain
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
   level_link => mnachain%last_link
   next_level_link => new_chain_link(mnachain)

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

   ! Compute self-consistent MNAs
   do
      ! Call compute_mna_partition and get the number of splits
      call compute_mna_partition(atoms1, atoms2, mnachain, branch, branch_parts, num_splits)

      ! Exit loop if no splits occurred
      if (num_splits == 0) exit
   end do

   ! Find a degenerate descendant part to split
   next_part_to_split => null()
   partref => mnachain%last_link%first_partref
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
      call split_dependent_parts(atoms1, atoms2, mnachain, branch, branch_parts, branching_part, next_part_to_split)
   end if
end subroutine

! Updated split_independent_parts to use the merged function signature
recursive subroutine split_independent_parts(atoms1, atoms2, mnachain, branch, branch_parts)
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(assigntree_node_t), pointer, intent(inout) :: mnachain, branch
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
         call split_dependent_parts(atoms1, atoms2, mnachain, new_branch, new_branch_parts, partref%part, partref%part)
         ! Recursively process the resulting branch parts
         call split_independent_parts(atoms1, atoms2, mnachain, new_branch, new_branch_parts)
      end if
      partref => partref%nextref
   end do
end subroutine

subroutine build_assignment_tree( atoms1, atoms2, mnalink, assign_frame)
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(chain_node_t), pointer, intent(in) :: mnalink
   type(array_trees_t), intent(out) :: assign_frame
   ! Local variables
   type(partree_node_t), pointer :: part_tree
   type(assigntree_node_t), pointer :: assign_tree
   type(assigntree_node_t), pointer :: mnachain
   type(chain_node_t), pointer :: branch_parts
   type(partree_node_t), pointer :: child_part
   type(chain_node_t), pointer :: first_link
   type(partref_node_t), pointer :: partref

   part_tree => new_root_part()
   branch_parts => new_bare_link()
   assign_tree => new_root_chain( size(mnalink%itemdir1), size(mnalink%itemdir2))
   mnachain => new_root_chain( size(mnalink%itemdir1), size(mnalink%itemdir2))
   first_link => new_chain_link( mnachain)

   partref => mnalink%first_partref
   do while (associated( partref))
      child_part => new_child_part( part_tree)
      call link_part( first_link, child_part)
      call copy_part_items( partref%part, child_part)
      call add_branch_part( branch_parts, child_part)
      partref => partref%nextref
   end do

   call split_independent_parts( atoms1, atoms2, mnachain, assign_tree, branch_parts)
!block
!   use assigntree_distribute_linked
!   call distribute_items( atoms1, atoms2, assign_tree)
!end block
   call convert_trees_to_arrays( atoms1, atoms2, part_tree, assign_tree, assign_frame)
!   call validate_conversion(part_tree, assign_tree, assign_frame)
end subroutine

end module
