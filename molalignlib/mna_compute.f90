module mna_compute
use parameters
use molecule
use lcrs_tree
implicit none

contains

subroutine compute_nextlevel_mnas(mol1, mol2, mnachain)
! Compute next level MNA types
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_node_t), pointer, intent(in) :: mnachain
   ! Local variables
   type(link_node_t), pointer :: last_link, new_link
   type(part_node_t), pointer :: part, child_part
   type(item_node_t), pointer :: item
   type(part_nodeptr_t), dimension(:), allocatable :: neighbors

   ! Create a new branch as a child of the current branch
   last_link => mnachain%last_link
   new_link => add_new_link(mnachain)

   part => last_link%first_part
   do while (associated(part))
      ! First molecule
      item => part%first_item1
      do while (associated(item))
         neighbors = last_link%itemdir1(mol1%atoms(item%value)%adjlist)
         child_part => find_child_part_node(part, neighbors)
         if (.not. associated(child_part)) then
            child_part => add_new_part(new_link, part, neighbors)
         end if
         call add_new_item1(child_part, item%value)
         item => item%next_item
      end do
      ! Second molecule
      item => part%first_item2
      do while (associated(item))
         neighbors = last_link%itemdir2(mol2%atoms(item%value)%adjlist)
         child_part => find_child_part_node(part, neighbors)
         if (.not. associated(child_part)) then
            child_part => add_new_part(new_link, part, neighbors)
         end if
         call add_new_item2(child_part, item%value)
         item => item%next_item
      end do
      part => part%next_part
   end do
end subroutine

subroutine compute_consistent_mnas(mol1, mol2, mnachain)
! Iteratively compute MNA types
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_node_t), pointer, intent(inout) :: mnachain
   integer :: prev_num_parts

   do
      ! Save the current number of parts before computation
      prev_num_parts = mnachain%last_link%num_parts

      ! Compute next level and update current branch
      call compute_nextlevel_mnas(mol1, mol2, mnachain)

      ! Exit the loop if no change
      if (mnachain%last_link%num_parts == prev_num_parts) exit
   end do
end subroutine

subroutine copy_partition_items(last_link, new_link)
   type(link_node_t), pointer, intent(in) :: last_link
   type(link_node_t), pointer, intent(inout) :: new_link
   ! Local variables
   type(part_node_t), pointer :: part
   type(part_node_t), pointer :: child_part
   type(item_node_t), pointer :: item1, item2

   part => last_link%first_part
   do while (associated(part))
      child_part => add_new_part(new_link, part)
      item1 => part%first_item1
      item2 => part%first_item2
      do while (associated(item1))
         call add_new_item1(child_part, item1%value)
         call add_new_item2(child_part, item2%value)
         item1 => item1%next_item
         item2 => item2%next_item
      end do
      part => part%next_part
   end do
end subroutine

function isancestor(part, branch_part)
   type(part_node_t), pointer, intent(in) :: part, branch_part
   ! Local variables
   logical :: isancestor
   type(part_node_t), pointer :: ancestor

   ancestor => part%parent
   do while (associated(ancestor))
      if (associated(ancestor, branch_part)) then
         isancestor = .true.
         return
      end if
      ancestor => ancestor%parent
   end do

   isancestor = .false.
end function

subroutine assign_partition_item(new_link, branch_part, assigned)
   type(link_node_t), pointer, intent(inout) :: new_link
   type(part_node_t), pointer, intent(in) :: branch_part
   logical, intent(out) :: assigned
   ! Local variables
   type(part_node_t), pointer :: part, child_part

   part => new_link%first_part
   do while (associated(part))
      ! Check if this leaf needs processing
      if (part%num_items1 >= 2) then
         if (isancestor(part, branch_part)) then
            ! Create a new part and move the first item pair to it
            child_part => add_new_part(new_link, part%parent)
            call move_first_item1(part, child_part)
            call move_first_item2(part, child_part)
            assigned = .true.
            return
         end if
      end if
      part => part%next_part
   end do
   assigned = .false.
end subroutine

subroutine assign_chain_item(new_branch, branch_part, assigned)
   type(tree_node_t), pointer, intent(inout) :: new_branch
   type(part_node_t), pointer, intent(in) :: branch_part
   logical, intent(out) :: assigned
   ! Local variables
   type(link_node_t), pointer :: last_link, new_link

   last_link => new_branch%last_link
   new_link => add_new_link(new_branch)

   call copy_partition_items(last_link, new_link)
   call assign_partition_item(new_link, branch_part, assigned)
end subroutine

subroutine split_mnas(mol1, mol2, mnatree)
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_node_t), pointer, intent(inout) :: mnatree
   ! Local variables
   type(tree_node_t), pointer :: new_branch
   type(link_node_t), pointer :: new_link
   type(part_node_t), pointer :: part
   logical :: assigned

   ! Sort the linked list of parts in place
   call sort_parts_by_size(mnatree%last_link)

   part => mnatree%last_link%first_part
   do while (associated(part))
      if (part%num_items1 >= 2) then
         if (part%num_leaves == 1) then
            ! Create a new branch for the current part
            new_branch => add_new_branch(mnatree)
            new_link => add_new_link(new_branch)
            call copy_partition_items(mnatree%last_link, new_link)
            call assign_partition_item(new_link, part, assigned)
            ! Process current part
            do while (assigned)
               call compute_consistent_mnas(mol1, mol2, new_branch)
               call assign_chain_item(new_branch, part, assigned)
            end do
         end if
      end if
      part => part%next_part
   end do
end subroutine

end module
