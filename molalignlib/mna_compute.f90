module mna_compute
use parameters
use molecule
use lcrs_tree
implicit none

contains

function would_part_split(mol1, mol2, leaf_part) result(would_split)
! Check if a part would split by comparing signatures
   type(mol_type), intent(in) :: mol1, mol2
   type(part_node_t), pointer, intent(in) :: leaf_part
   logical :: would_split
   type(item_node_t), pointer :: item1, item2
   type(part_nodeptr_t), dimension(:), allocatable :: signature, first_signature

   would_split = .false.
   item1 => leaf_part%first_item1
   item2 => leaf_part%first_item2

   ! Set first signature from first available item
   if (associated(item1)) then
      first_signature = leaf_part%partition_root%itemdir1(mol1%atoms(item1%value)%adjlist)
      item1 => item1%next_item
   else if (associated(item2)) then
      first_signature = leaf_part%partition_root%itemdir2(mol2%atoms(item2%value)%adjlist)
      item2 => item2%next_item
   else
      return  ! No items to process
   end if

   ! Check remaining items in first molecule
   do while (associated(item1))
      signature = leaf_part%partition_root%itemdir1(mol1%atoms(item1%value)%adjlist)
      if (.not. (signature .equiv. first_signature)) then
         would_split = .true.
         return
      end if
      item1 => item1%next_item
   end do

   ! Check remaining items in second molecule
   do while (associated(item2))
      signature = leaf_part%partition_root%itemdir2(mol2%atoms(item2%value)%adjlist)
      if (.not. (signature .equiv. first_signature)) then
         would_split = .true.
         return
      end if
      item2 => item2%next_item
   end do
end function

subroutine split_part(mol1, mol2, part, link)
! Split a part by creating children with different signatures
   type(mol_type), intent(in) :: mol1, mol2
   type(part_node_t), pointer, intent(inout) :: part
   type(link_node_t), pointer, intent(inout) :: link
   type(item_node_t), pointer :: item
   type(part_nodeptr_t), dimension(:), allocatable :: signature
   type(part_node_t), pointer :: child_part

   ! Process first molecule items
   item => part%first_item1
   do while (associated(item))
      signature = part%partition_root%itemdir1(mol1%atoms(item%value)%adjlist)
      child_part => find_child_part_node(part, signature)
      if (.not. associated(child_part)) then
         child_part => add_new_part(link, part)
         child_part%signature = signature
      end if
      call add_new_item1(child_part, item%value)
      item => item%next_item
   end do

   ! Process second molecule items
   item => part%first_item2
   do while (associated(item))
      signature = part%partition_root%itemdir2(mol2%atoms(item%value)%adjlist)
      child_part => find_child_part_node(part, signature)
      if (.not. associated(child_part)) then
         child_part => add_new_part(link, part)
         child_part%signature = signature
      end if
      call add_new_item2(child_part, item%value)
      item => item%next_item
   end do
end subroutine

subroutine compute_nextlevel_mnas(mol1, mol2, chain_root)
! Compute next level MNA types by processing the last partition
   type(mol_type), intent(in) :: mol1, mol2
   type(chain_root_t), pointer, intent(inout) :: chain_root
   type(link_node_t), pointer :: last_link, new_link
   type(partref_node_t), pointer :: partref

   ! Get the last link in the chain
   last_link => chain_root%last_link

   ! Create a new link for the next level
   new_link => add_new_link(chain_root)

   ! Process all parts in the current partition
   partref => last_link%first_partref
   do while (associated(partref))
      ! All parts in the last link are leaf parts - split if necessary, otherwise reuse
      if (would_part_split(mol1, mol2, partref%part_ptr)) then
         call split_part(mol1, mol2, partref%part_ptr, new_link)
      else
         call add_part(new_link, partref%part_ptr)
      end if
      partref => partref%next_partref
   end do
end subroutine

subroutine compute_consistent_mnas(mol1, mol2, chain_root)
! Iteratively compute MNA types until convergence
   type(mol_type), intent(in) :: mol1, mol2
   type(chain_root_t), pointer, intent(inout) :: chain_root
   integer :: prev_num_parts

   do
      prev_num_parts = chain_root%last_link%num_parts
      call compute_nextlevel_mnas(mol1, mol2, chain_root)
      if (chain_root%last_link%num_parts == prev_num_parts) exit
   end do
end subroutine

subroutine split_first_item(part, new_link)
   type(part_node_t), pointer, intent(inout) :: part
   type(link_node_t), pointer, intent(inout) :: new_link
   ! Local variables
   type(part_node_t), pointer :: child_part
   type(item_node_t), pointer :: item1, item2

   child_part => add_new_part(new_link, part)
   item1 => part%first_item1
   item2 => part%first_item2
   call add_new_item1(child_part, item1%value)
   call add_new_item2(child_part, item2%value)
   item1 => item1%next_item
   item2 => item2%next_item
   child_part => add_new_part(new_link, part)
   do while (associated(item1))
      call add_new_item1(child_part, item1%value)
      call add_new_item2(child_part, item2%value)
      item1 => item1%next_item
      item2 => item2%next_item
   end do
end subroutine

subroutine copy_items(part, new_link)
   type(part_node_t), pointer, intent(inout) :: part
   type(link_node_t), pointer, intent(inout) :: new_link
   ! Local variables
   type(part_node_t), pointer :: child_part
   type(item_node_t), pointer :: item1, item2

   child_part => add_new_part(new_link, part)
   item1 => part%first_item1
   item2 => part%first_item2
   do while (associated(item1))
      call add_new_item1(child_part, item1%value)
      call add_new_item2(child_part, item2%value)
      item1 => item1%next_item
      item2 => item2%next_item
   end do
end subroutine

subroutine split_mnas(mol1, mol2, mnachain)
   type(mol_type), intent(in) :: mol1, mol2
   type(chain_root_t), pointer, intent(inout) :: mnachain
   ! Local variables
   type(link_node_t), pointer :: new_link
   type(partref_node_t), pointer :: partref

   partref => mnachain%last_link%first_partref
   new_link => add_new_link(mnachain)
   call split_first_item(partref%part_ptr, new_link)
   partref => partref%next_partref
   do while (associated(partref))
      call copy_items(partref%part_ptr, new_link)
      partref => partref%next_partref
   end do
   call compute_consistent_mnas(mol1, mol2, mnachain)
end subroutine

end module
