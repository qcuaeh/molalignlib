module metapartitioning
use parameters
use molecule
use lcrs_tree
use partitioning
implicit none

contains

subroutine collect_mnatypes(mol1, root)
   type(mol_type), intent(in) :: mol1
   type(tree_node), intent(inout) :: root
!   type(tree_node), pointer, intent(out) :: metatypes
   ! Local variables
    type(tree_node), pointer :: child
!   type(partition_type) :: subtypes
!   type(tree_node), pointer :: inode, jnode

!   metatypes => new_tree()

!   do h = 1, mnatypes%num_parts
!      if (mnatypes%parts(h)%part_size > 1) then
!         subtypes = mnatypes
!         call split_mnatype(h, subtypes)
!         call compute_consistent_mnatypes(mol, subtypes)
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
         call recompute_consistent_mnatypes1(mol1, root)
         call print_subtree1(root, 1)
         stop
      end if
      child => child%next_sibling
   end do

!   inode => metatypes%first_child
!   loop1: do while (associated(inode))
!      jnode => inode%next_sibling
!      do while (associated(jnode))
!         if (inode%subtypes < jnode%subtypes) then
!            call delete_branch(jnode)
!         else if (jnode%subtypes < inode%subtypes) then
!            call delete_branch(inode)
!            cycle loop1
!         else
!            jnode => jnode%next_sibling
!         end if
!      end do
!      inode => inode%next_sibling
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

end module
