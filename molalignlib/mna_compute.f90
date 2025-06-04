module mna_compute
use parameters
use molecule
use lcrs_tree
use eltype_compute
implicit none

contains

subroutine split_part_mna(atoms1, atoms2, itemdir1, itemdir2, part)
! Create children for different signatures - caller decides what to do with them
! Note: part is always a leaf part with no existing children
   type(atom_type), dimension(:), intent(in) :: atoms1, atoms2
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir1, itemdir2
   type(part_node_t), pointer, intent(inout) :: part
   ! Local variables
   type(item_node_t), pointer :: item
   type(part_nodeptr_t), dimension(:), allocatable :: signature
   type(part_node_t), pointer :: child_part

   ! Process first molecule items - create children for each unique signature
   item => part%first_item1
   do while (associated(item))
      signature = itemdir1(atoms1(item%value)%adjlist)
      child_part => find_child_part(part, signature)
      if (.not. associated(child_part)) then
         child_part => new_child_part(part)
         child_part%signature = signature
      end if
      call add_new_item1(child_part, item%value)
      item => item%next_item
   end do

   ! Process second molecule items - create children for each unique signature
   item => part%first_item2
   do while (associated(item))
      signature = itemdir2(atoms2(item%value)%adjlist)
      child_part => find_child_part(part, signature)
      if (.not. associated(child_part)) then
         child_part => new_child_part(part)
         child_part%signature = signature
      end if
      call add_new_item2(child_part, item%value)
      item => item%next_item
   end do
end subroutine

subroutine compute_nextlevel_mnas(mol1, mol2, mnachain, num_splits)
! Compute next level MNA types - always keeps all children (original behavior)
   type(mol_type), intent(in) :: mol1, mol2
   type(branch_node_t), pointer, intent(inout) :: mnachain
   integer, intent(out) :: num_splits
   ! Local variables
   type(link_node_t), pointer :: last_link, new_link
   type(partref_node_t), pointer :: partref

   num_splits = 0

   ! Save the last link before creating a new one
   last_link => mnachain%last_link
   new_link => new_chain_link(mnachain)

   ! Process all parts in the current partition
   partref => last_link%first_partref
   do while (associated(partref))
      ! Create children based on signatures
      call split_part_mna(mol1%atoms, mol2%atoms, last_link%itemdir1, last_link%itemdir2, partref%part)

      ! Always link all children to new_link (preserves original behavior)
      call update_link_parts(new_link, partref%part)

      ! Count splits (children beyond the original part)
      num_splits = num_splits + partref%part%num_children - 1

      partref => partref%nextref
   end do
end subroutine

subroutine compute_nextlevel_mnas_branch(mol1, mol2, mnachain, mnabranch, branch_parts, num_splits)
! Compute next level MNA types - only keeps children if real split occurred
   type(mol_type), intent(in) :: mol1, mol2
   type(branch_node_t), pointer, intent(inout) :: mnachain
   type(branch_node_t), pointer, intent(inout) :: mnabranch
   type(link_node_t), pointer, intent(inout) :: branch_parts
   integer, intent(out) :: num_splits
   ! Local variables
   type(link_node_t), pointer :: last_link, new_link, leaf_link
   type(partref_node_t), pointer :: partref

   num_splits = 0

   ! Save the last link before creating a new one
   last_link => mnachain%last_link
   new_link => new_chain_link(mnachain)
   leaf_link => new_chain_link(mnabranch)

   ! Process all parts in the current partition (all are leaf parts)
   partref => last_link%first_partref
   do while (associated(partref))
      ! Create children based on signatures
      call split_part_mna(mol1%atoms, mol2%atoms, last_link%itemdir1, last_link%itemdir2, partref%part)

      ! Check if real split occurred (more than one child)
      if (partref%part%num_children > 1) then
         ! Real split - link all children to new_link
         call update_link_parts(new_link, partref%part)
         call add_link_part(leaf_link, partref%part)
!         call update_link_parts(leaf_link, partref%part)
         call update_branch_parts(branch_parts, partref%part)
         num_splits = num_splits + 1
      else
         ! No real split (num_children == 1) - remove the single child and reuse original part
         call remove_onlychild_part(partref%part)
         call add_link_part(new_link, partref%part)
      end if

      partref => partref%nextref
   end do

   ! If no splits occurred, remove the newly created links
   if (num_splits == 0) then
      call remove_last_link(mnachain)
      call remove_last_link(mnabranch)
   end if
end subroutine

subroutine compute_consistent_mnas(mol1, mol2, mnachain)
! Iteratively compute MNA types until convergence
   type(mol_type), intent(in) :: mol1, mol2
   type(branch_node_t), pointer, intent(inout) :: mnachain
   ! Local variables
   integer :: num_splits

   do
      ! Call compute_nextlevel_mnas and get the number of splits
      call compute_nextlevel_mnas(mol1, mol2, mnachain, num_splits)

      ! Exit loop if no splits occurred in the last iteration
      if (num_splits == 0) exit
   end do
end subroutine

subroutine compute_consistent_mnas_branch(mol1, mol2, mnachain, mnabranch, branch_parts)
! Iteratively compute MNA types until convergence
   type(mol_type), intent(in) :: mol1, mol2
   type(branch_node_t), pointer, intent(inout) :: mnachain
   type(branch_node_t), pointer, intent(inout) :: mnabranch
   type(link_node_t), pointer, intent(inout) :: branch_parts
   ! Local variables
   integer :: num_splits

   do
      ! Call compute_nextlevel_mnas and get the number of splits
      call compute_nextlevel_mnas_branch(mol1, mol2, mnachain, mnabranch, branch_parts, num_splits)

      ! Exit loop if no splits occurred in the last iteration
      if (num_splits == 0) exit
   end do
end subroutine

subroutine split_part_item(part)
   type(part_node_t), pointer, intent(inout) :: part
   ! Local variables
   type(part_node_t), pointer :: child_part
   type(item_node_t), pointer :: item1, item2

   ! Create first child and add first item from each molecule
   child_part => new_child_part(part)
   call add_new_item1(child_part, part%first_item1%value)
   call add_new_item2(child_part, part%first_item2%value)

   ! Create second child and add remaining items
   child_part => new_child_part(part)
   item1 => part%first_item1%next_item
   item2 => part%first_item2%next_item
   do while (associated(item1))
      call add_new_item1(child_part, item1%value)
      call add_new_item2(child_part, item2%value)
      item1 => item1%next_item
      item2 => item2%next_item
   end do
end subroutine

! Modified split_single_part as a function that returns the created child branch
function split_single_part(mnachain, branch, part, branch_parts) result(child_branch)
   type(branch_node_t), pointer, intent(inout) :: mnachain
   type(branch_node_t), pointer, intent(inout) :: branch
   type(part_node_t), pointer, intent(inout) :: part
   type(link_node_t), pointer, intent(inout) :: branch_parts
   type(branch_node_t), pointer :: child_branch
   ! Local variables
   type(link_node_t), pointer :: last_link, new_link, leaf_link
   type(partref_node_t), pointer :: partref

   ! Save the last link before creating a new one
   last_link => mnachain%last_link
   new_link => new_chain_link(mnachain)

   ! Add non-target parts to new link
   partref => last_link%first_partref
   do while (associated(partref))
      if (.not. associated(partref%part, part)) then
         call add_link_part(new_link, partref%part)
      end if
      partref => partref%nextref
   end do

   ! Split target part and add its children to the new link
   call split_part_item(part)
   call update_link_parts(new_link, part)

   ! Create a new child branch for this leaf part
   child_branch => new_child_branch(branch)
   leaf_link => new_chain_link(child_branch)
   
   ! Update branch parts
   call add_link_part(leaf_link, part)
   call update_branch_parts(branch_parts, part)
end function

! Modified split_dependent_parts as a function that returns the final branch
recursive function split_dependent_parts(mol1, mol2, mnachain, branch, fork_part, branch_parts) result(branch_tip)
   type(mol_type), intent(in) :: mol1, mol2
   type(branch_node_t), pointer, intent(inout) :: mnachain
   type(branch_node_t), pointer, intent(inout) :: branch
   type(part_node_t), pointer, intent(in) :: fork_part
   type(link_node_t), pointer, intent(inout) :: branch_parts
   type(branch_node_t), pointer :: branch_tip
   ! Local variables
   type(partref_node_t), pointer :: partref
   type(part_node_t), pointer :: split_part

   ! Start with the input branch
   branch_tip => branch

   ! Compute consistent MNAs
   call compute_consistent_mnas_branch(mol1, mol2, mnachain, branch_tip, branch_parts)

   ! Find a degenerate descendant part to split
   split_part => null()
   partref => mnachain%last_link%first_partref
   do while (associated(partref) .and. .not. associated(split_part))
      if (partref%part%num_items1 >= 2) then
         if (isdescendant(partref%part, fork_part)) then
            split_part => partref%part
         end if
      end if
      partref => partref%nextref
   end do

   ! Perform split if target found
   if (associated(split_part)) then
      ! Split the target part and get the new child branch
      branch_tip => split_single_part(mnachain, branch_tip, split_part, branch_parts)
      ! Call itself again to split the next degenerate descendant part
      branch_tip => split_dependent_parts(mol1, mol2, mnachain, branch_tip, fork_part, branch_parts)
   end if
end function

! Updated split_independent_parts to use the new function signatures
recursive subroutine split_independent_parts(mol1, mol2, mnachain, branch, branch_parts)
   type(mol_type), intent(in) :: mol1, mol2
   type(branch_node_t), pointer, intent(inout) :: mnachain, branch
   type(link_node_t), pointer, intent(in) :: branch_parts
   ! Local variables
   type(link_node_t), pointer :: new_branch_parts
   type(branch_node_t), pointer :: branch_tip
   type(partref_node_t), pointer :: partref

   ! Process each part in branch_parts
   partref => branch_parts%first_partref
   do while (associated(partref))
      if (partref%part%num_children == 0) then
         ! Create a new part registry for this branch part
         new_branch_parts => new_root_link()
         ! Split the target part and get the new child branch
         branch_tip => split_single_part(mnachain, branch, partref%part, new_branch_parts)
         ! Split items for this specific part until convergence
         branch_tip => split_dependent_parts(mol1, mol2, mnachain, branch_tip, partref%part, new_branch_parts)
         ! Recursively process the resulting branch parts
         call split_independent_parts(mol1, mol2, mnachain, branch_tip, new_branch_parts)
      end if
      partref => partref%nextref
   end do
end subroutine

subroutine assign_conform_atoms( mol1, mol2, root_part, mnachain, root_branch)
   type(mol_type), intent(in) :: mol1, mol2
   type(part_node_t), pointer, intent(inout) :: root_part
   type(branch_node_t), pointer, intent(inout) :: mnachain
   type(branch_node_t), pointer, intent(out) :: root_branch
   type(link_node_t), pointer :: leaf_link, branch_parts
   ! Local variables

!   call print_atoms( mol1)
!   call print_atoms( mol2)

   root_branch => new_root_branch( mnachain%tot_items1, mnachain%tot_items2)
   leaf_link => new_chain_link( root_branch)
   call add_link_part(leaf_link, root_part)
   branch_parts => new_root_link()

   call compute_consistent_mnas_branch( mol1, mol2, mnachain, root_branch, branch_parts)
   call split_independent_parts( mol1, mol2, mnachain, root_branch, branch_parts)

!   call print_chain( mnachain)
   call print_part_tree( root_part)
   call print_branch_tree( root_branch)
end subroutine

end module
