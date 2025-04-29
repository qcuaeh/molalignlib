module metapartitioning
use parameters
use molecule
use lcrs_tree
implicit none

contains

subroutine collect_mnas(mol1, root)
   type(mol_type), intent(in) :: mol1
   type(tree_node), intent(inout) :: root
!   type(tree_node), pointer, intent(out) :: metatypes
   ! Local variables
    type(tree_node), pointer :: child
!   type(item_partition) :: subtypes
!   type(tree_node), pointer :: node, jnode

!   metatypes => new_tree()

!   do h = 1, mnatypes%num_parts
!      if (mnatypes%parts(h)%part_size > 1) then
!         subtypes = mnatypes
!         call split_mnatype(h, subtypes)
!         call compute_consistent_mnas(mol, subtypes)
!         if (.not. (subtypes .in. partitiondict)) then
!            partlist(partitiondict%new_index(subtypes))%ptr => &
!               metatypes%new_leaf(subtypes)
!         end if
!         call partlist(partitiondict%get_index(subtypes))%ptr%add_item(h)
!      end if
!   end do

   child => root%first_child
   do while (associated(child))
      if (child%num_items1 > 1) then
         call split_mnatype(child)
         call compute_consistent_mnatypes1(mol1, root)
         call print_subtree1(root, 1)
         stop
      end if
      child => child%next_sibling
   end do

!   node => metatypes%first_child
!   loop1: do while (associated(node))
!      jnode => node%next_sibling
!      do while (associated(jnode))
!         if (node%subtypes < jnode%subtypes) then
!            call delete_branch(jnode)
!         else if (jnode%subtypes < node%subtypes) then
!            call delete_branch(node)
!            cycle loop1
!         else
!            jnode => jnode%next_sibling
!         end if
!      end do
!      node => node%next_sibling
!   end do loop1
end subroutine

subroutine split_mnatype(leaf)
   type(tree_node), intent(inout) :: leaf
   type(tree_node), pointer :: child

   do while (associated(leaf%first_item1))
      child => add_new_child(leaf)
      call move_next_item1(leaf, child)
      call move_next_item2(leaf, child)
   end do
end subroutine

subroutine compute_consistent_mnatypes1(mol1, mnatypetree)
! Iteratively recompute MNA types

   type(mol_type), intent(in) :: mol1
   type(tree_node), intent(inout) :: mnatypetree
   ! Local variables
   type(tree_node_ptr), dimension(:), allocatable :: itemdir1

   do

      itemdir1 = mnatypetree%itemdir1

      ! Compute MNA upper level types
      call compute_nextlevel_mnatypes1(mol1, itemdir1, mnatypetree)
!      call print_tree(mnatypetree)

      ! Exit loop if types did not change
      if (all(mnatypetree%itemdir1 == itemdir1)) exit

   end do
end subroutine

recursive subroutine compute_nextlevel_mnatypes1(mol1, itemdir1, node)
! Recompute next level MNA types

   type(mol_type), intent(in) :: mol1
   type(tree_node_ptr), dimension(:), intent(in) :: itemdir1
   type(tree_node), intent(inout) :: node
   ! Local variables
   type(tree_node), pointer :: child

   if (associated(node%first_child)) then
      ! Internal node - process children
      child => node%first_child
      do
         call compute_nextlevel_mnatypes1(mol1, itemdir1, child)
         if (.not. associated(child%next_sibling)) return
         child => child%next_sibling
      end do
   else
      ! Leaf node
      if (node%num_items1 + node%num_items2 > 1) then
         ! Multiple items - update MNA type
         call compute_mnatype1(mol1, itemdir1, node)
      end if
   end if
end subroutine

subroutine compute_mnatype1(mol1, itemdir1, node)
   type(mol_type), intent(in) :: mol1
   type(tree_node_ptr), dimension(:), intent(in) :: itemdir1
   type(tree_node), intent(inout) :: node
   ! Local variables
   type(tree_node), pointer :: child
   type(tree_node_ptr), dimension(:), allocatable :: typehood

   allocate (node%typehoodtable%items(node%num_items1))
   node%typehoodtable%num_items = 0

   ! First molecule
   do while (associated(node%first_item1))
      typehood = itemdir1(mol1%atoms(node%first_item1%index)%adjlist)
      child => find_heir(node, typehood)
      if (.not. associated(child)) then
         child => add_new_child(node)
         call add_heir(node, child, typehood)
      end if
      call move_next_item1(node, child)
   end do

   ! Revert changes if only child
   if (node%num_childs == 1) then
      call move_node_items(node%first_child, node)
      deallocate (node%first_child)
      node%num_childs = 0
   end if
end subroutine

end module
