module mna_compute
use parameters
use molecule
use lcrs_tree
implicit none

contains

subroutine split_part_mna(atoms1, atoms2, itemdir1, itemdir2, part, link)
! Create children for different signatures - caller decides what to do with them
! Note: part is always a leaf part with no existing children
   type(atom_type), dimension(:), intent(in) :: atoms1, atoms2
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir1, itemdir2
   type(part_node_t), pointer, intent(inout) :: part
   type(link_node_t), pointer, intent(inout) :: link
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
         call link_part(link, child_part)
      end if
      call add_new_item1(child_part, item%value)
      link%itemdir1(item%value)%ptr => child_part
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
         call link_part(link, child_part)
      end if
      call add_new_item2(child_part, item%value)
      link%itemdir2(item%value)%ptr => child_part
      item => item%next_item
   end do
end subroutine

subroutine compute_nextlevel_mnas(mol1, mol2, mnachain, num_splits)
! Compute next level MNA types - always keeps all children (original behavior)
   type(mol_type), intent(in) :: mol1, mol2
   type(chain_node_t), pointer, intent(inout) :: mnachain
   integer, intent(out) :: num_splits
   ! Local variables
   type(link_node_t), pointer :: link, new_link
   type(partref_node_t), pointer :: partref

   num_splits = 0

   ! Save the last link before creating a new one
   link => mnachain%last_link
   new_link => new_chain_link(mnachain)

   ! Process all parts in the current partition
   partref => link%first_partref
   do while (associated(partref))
      ! Create children based on signatures
      call split_part_mna(mol1%atoms, mol2%atoms, link%itemdir1, link%itemdir2, partref%part, new_link)

      ! Count splits (children beyond the original part)
      num_splits = num_splits + partref%part%num_children - 1

      partref => partref%nextref
   end do
end subroutine

subroutine compute_consistent_mnas(mol1, mol2, mnachain)
! Iteratively compute MNA types until convergence
   type(mol_type), intent(in) :: mol1, mol2
   type(chain_node_t), pointer, intent(inout) :: mnachain
   ! Local variables
   integer :: num_splits

   do
      ! Call compute_nextlevel_mnas and get the number of splits
      call compute_nextlevel_mnas(mol1, mol2, mnachain, num_splits)

      ! Exit loop if no splits occurred
      if (num_splits == 0) exit
   end do
end subroutine

end module
