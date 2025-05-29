module mna_compute
use parameters
use molecule
use lcrs_tree
use eltype_compute
implicit none

contains

function would_part_split(mol1, mol2, leaf_part) result(would_split)
! Check if a part would split by comparing signatures
   type(mol_type), intent(in) :: mol1, mol2
   type(part_node_t), pointer, intent(in) :: leaf_part
   ! Local variables
   logical :: would_split
   type(item_node_t), pointer :: item1, item2
   type(part_nodeptr_t), dimension(:), allocatable :: signature, first_signature

   would_split = .false.
   item1 => leaf_part%first_item1
   item2 => leaf_part%first_item2

   ! Set first signature from first available item
   if (associated(item1)) then
      first_signature = leaf_part%parent_link%itemdir1(mol1%atoms(item1%value)%adjlist)
      item1 => item1%next_item
   else if (associated(item2)) then
      first_signature = leaf_part%parent_link%itemdir2(mol2%atoms(item2%value)%adjlist)
      item2 => item2%next_item
   else
      return  ! No items to process
   end if

   ! Check remaining items in first molecule
   do while (associated(item1))
      signature = leaf_part%parent_link%itemdir1(mol1%atoms(item1%value)%adjlist)
      if (.not. (signature .equiv. first_signature)) then
         would_split = .true.
         return
      end if
      item1 => item1%next_item
   end do

   ! Check remaining items in second molecule
   do while (associated(item2))
      signature = leaf_part%parent_link%itemdir2(mol2%atoms(item2%value)%adjlist)
      if (.not. (signature .equiv. first_signature)) then
         would_split = .true.
         return
      end if
      item2 => item2%next_item
   end do
end function

subroutine split_part_mnas(mol1, mol2, part, link, branch_link)
! Split a part by creating children with different signatures
   type(mol_type), intent(in) :: mol1, mol2
   type(part_node_t), pointer, intent(inout) :: part
   type(link_node_t), pointer, intent(inout) :: link, branch_link
   ! Local variables
   type(item_node_t), pointer :: item
   type(part_nodeptr_t), dimension(:), allocatable :: signature
   type(part_node_t), pointer :: child_part

   ! Process first molecule items
   item => part%first_item1
   do while (associated(item))
      signature = part%parent_link%itemdir1(mol1%atoms(item%value)%adjlist)
      child_part => find_child_part_node(part, signature)
      if (.not. associated(child_part)) then
         child_part => add_new_part(link, part)
         child_part%signature = signature
         call add_part(branch_link, part)
      end if
      call add_new_item1(child_part, item%value)
      item => item%next_item
   end do

   ! Process second molecule items
   item => part%first_item2
   do while (associated(item))
      signature = part%parent_link%itemdir2(mol2%atoms(item%value)%adjlist)
      child_part => find_child_part_node(part, signature)
      if (.not. associated(child_part)) then
         child_part => add_new_part(link, part)
         child_part%signature = signature
         call add_part(branch_link, part)
      end if
      call add_new_item2(child_part, item%value)
      item => item%next_item
   end do
end subroutine

subroutine compute_nextlevel_mnas(mol1, mol2, mnachain, branch, num_splits)
! Compute next level MNA types by processing the last partition
   type(mol_type), intent(in) :: mol1, mol2
   type(chain_root_t), pointer, intent(inout) :: mnachain
   type(branch_node_t), pointer, intent(inout) :: branch
   integer, intent(out) :: num_splits ! New argument for split count
   ! Local variables
   type(link_node_t), pointer :: last_link, new_link, branch_link
   type(partref_node_t), pointer :: partref

   num_splits = 0 ! Initialize split count

   ! Get the last link in the chain
   last_link => mnachain%last_link

   ! Create a new link for the next level
   new_link => add_new_link(mnachain)
   branch_link => add_new_link_branch(branch)

   ! Process all parts in the current partition
   partref => last_link%first_partref
   do while (associated(partref))
      ! All parts in the last link are leaf parts - split if necessary, otherwise reuse
      if (would_part_split(mol1, mol2, partref%part_node)) then
         call split_part_mnas(mol1, mol2, partref%part_node, new_link, branch_link)
         call update_branch_parts(branch, partref%part_node)
         num_splits = num_splits + 1 ! Increment split count
      else
         call add_part(new_link, partref%part_node)
      end if
      partref => partref%next_partref
   end do
end subroutine

subroutine compute_consistent_mnas(mol1, mol2, mnachain, branch)
! Iteratively compute MNA types until convergence
   type(mol_type), intent(in) :: mol1, mol2
   type(chain_root_t), pointer, intent(inout) :: mnachain
   type(branch_node_t), pointer, intent(inout) :: branch
   ! Local variables
   integer :: num_splits

   do
      ! Call compute_nextlevel_mnas and get the number of splits
      call compute_nextlevel_mnas(mol1, mol2, mnachain, branch, num_splits)

      ! Exit loop if no splits occurred in the last iteration
      if (num_splits == 0) exit
   end do
end subroutine

subroutine split_part_item(link, branch_link, part)
   type(link_node_t), pointer, intent(inout) :: link, branch_link
   type(part_node_t), pointer, intent(inout) :: part
   ! Local variables
   type(part_node_t), pointer :: child_part
   type(item_node_t), pointer :: item1, item2

   child_part => add_new_part(link, part)
   call add_part(branch_link, part)
   call add_new_item1(child_part, part%first_item1%value)
   call add_new_item2(child_part, part%first_item2%value)
   child_part => add_new_part(link, part)
   call add_part(branch_link, part)
   item1 => part%first_item1%next_item
   item2 => part%first_item2%next_item
   do while (associated(item1))
      call add_new_item1(child_part, item1%value)
      call add_new_item2(child_part, item2%value)
      item1 => item1%next_item
      item2 => item2%next_item
   end do
end subroutine

function should_split(down_part, top_part)
   type(part_node_t), pointer, intent(in) :: down_part, top_part
   ! Local variables
   logical :: should_split
   type(part_node_t), pointer :: up_part

   if (down_part%num_items1 >= 2) then
      up_part => down_part
      do while (up_part%depth >= top_part%depth)
         if (associated(up_part, top_part)) then
            should_split = .true.
            return
         end if
         up_part => up_part%parent_part
      end do
   end if

   should_split = .false.
end function

recursive subroutine split_branch_items(mol1, mol2, mnachain, branch)
   type(mol_type), intent(in) :: mol1, mol2
   type(chain_root_t), pointer, intent(inout) :: mnachain
   type(branch_node_t), pointer, intent(inout) :: branch
   ! Local variables
   type(link_node_t), pointer :: new_link, branch_link
   type(partref_node_t), pointer :: partref

   ! Get parts from current last link
   partref => mnachain%last_link%first_partref

   ! Create new link (this updates mnachain%last_link)
   new_link => add_new_link(mnachain)
   branch_link => add_new_link_branch(branch)

   ! Process all parts, looking for one to split
   do while (associated(partref))
      if (should_split(partref%part_node, branch%part_node)) then
         ! Found a part to split
         call split_part_item(new_link, partref%part_node, branch_link)
         call update_branch_parts(branch, partref%part_node)

         ! Add all remaining parts
         partref => partref%next_partref
         do while (associated(partref))
            call add_part(new_link, partref%part_node)
            partref => partref%next_partref
         end do

         ! Compute consistent MNAs and recurse
         call compute_consistent_mnas(mol1, mol2, mnachain, branch)
         call split_branch_items(mol1, mol2, mnachain, branch)
         return
      end if

      ! No split needed, just add this part
      call add_part(new_link, partref%part_node)
      partref => partref%next_partref
   end do

   ! Base case: no splits found, recursion ends
end subroutine

recursive subroutine build_permutation_tree( mol1, mol2, mnachain, branch)
   type(mol_type), intent(in) :: mol1, mol2
   type(chain_root_t), pointer, intent(out) :: mnachain
   type(branch_node_t), pointer, intent(out) :: branch
   ! Local variables
   type(partref_node_t), pointer :: partref
   type(branch_node_t), pointer :: new_branch

   ! Traverse all parts referenced by this branch
   partref => branch%first_partref
   do while (associated(partref))
      ! Only process leaf parts (parts with no children)
      if (partref%part_node%num_children == 0) then
         ! Create a new child branch for this leaf part
         new_branch => add_new_branch(branch, partref%part_node)

         ! Split items in this part if needed, creating new partition levels
         call split_branch_items(mol1, mol2, mnachain, new_branch)

         ! Recursively process the new branch to find more leaf parts
         call build_permutation_tree( mol1, mol2, mnachain, new_branch)
      end if
      partref => partref%next_partref
   end do

   ! Base case: when no leaf parts are found, recursion ends
end subroutine

subroutine assign_conform_atoms( mol1, mol2, mnachain, root_branch)
   type(mol_type), intent(in) :: mol1, mol2
   type(chain_root_t), pointer, intent(inout) :: mnachain
   type(branch_node_t), pointer, intent(inout) :: root_branch
   ! Local variables

!   call print_atoms( mol1)
!   call print_atoms( mol2)

   call compute_consistent_mnas( mol1, mol2, mnachain, root_branch)
   call build_permutation_tree( mol1, mol2, mnachain, root_branch)

!   call print_chain( mnachain)
   call print_part_tree( mnachain%first_link%first_partref%part_node)
   call print_branch_tree( root_branch)
   call print_branch_contents( root_branch)
end subroutine

end module
