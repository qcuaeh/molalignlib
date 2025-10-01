module hna
use parameters
use derived_types
use molecule
use lcrs_trees
implicit none
private
public refine_hna_part
public compute_hna_partition

contains

subroutine refine_hna_part(atoms1, atoms2, itemdir1, itemdir2, part, link)
! Create children for different signatures - caller decides what to do with them
! Note: part is always a leaf part with no existing children
   type(adjc_t), dimension(:), intent(in) :: atoms1, atoms2
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir1, itemdir2
   type(partree_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(item_node_t), pointer :: item
   type(partree_node_t), pointer :: child_part
   type(part_nodeptr_t), target :: signature_alloc(MAX_COORD)
   type(part_nodeptr_t), pointer :: signature(:)

   ! Process first molecule items - create children for each unique signature
   item => part%first_item1
   do while (associated(item))
      signature => signature_alloc(1:size(atoms1(item%idx)%adjlist))
      signature = itemdir1(atoms1(item%idx)%adjlist)
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
      signature => signature_alloc(1:size(atoms2(item%idx)%adjlist))
      signature = itemdir2(atoms2(item%idx)%adjlist)
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

subroutine refine_hna_partition(atoms1, atoms2, hnachain, num_splits)
! Compute next level HNA types - always keeps all children (original behavior)
   type(adjc_t), dimension(:), intent(in) :: atoms1, atoms2
   type(assigntree_node_t), pointer, intent(inout) :: hnachain
   integer, intent(out) :: num_splits
   ! Local variables
   type(chain_node_t), pointer :: link, new_link
   type(partref_node_t), pointer :: partref

   num_splits = 0

   ! Save the last link before creating a new one
   link => hnachain%last_link
   new_link => new_chain_link(hnachain)

   ! Process all parts in the current partition
   partref => link%first_partref
   do while (associated(partref))
      ! Create children based on signatures
      call refine_hna_part(atoms1, atoms2, link%itemdir1, link%itemdir2, partref%part, new_link)

      ! Count splits (children beyond the original part)
      num_splits = num_splits + partref%part%num_children - 1

      partref => partref%nextref
   end do
end subroutine

subroutine compute_hna_partition(atoms1, atoms2, atomtypes, hnachain)
! Iteratively compute HNA types until convergence
   type(adjc_t), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_t), intent(in) :: atomtypes
   type(assigntree_node_t), pointer, intent(out) :: hnachain
   ! Local variables
   integer :: num_splits

!   hnachain => collect_atomtypes_linked( atoms1, atoms2)
   hnachain => chain_from_partition( atomtypes)

   do
      ! Call refine_hna_partition and get the number of splits
      call refine_hna_partition(atoms1, atoms2, hnachain, num_splits)

      ! Exit loop if no splits occurred
      if (num_splits == 0) exit
   end do
end subroutine

end module
