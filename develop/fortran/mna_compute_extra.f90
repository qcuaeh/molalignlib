subroutine assign_next_item(mnatree, assigned)
   type(branch_node), pointer, intent(inout) :: mnatree
   logical, intent(inout) :: assigned
   type(branch_node), pointer :: current_branch, new_branch
   type(leaf_node), pointer :: parent_leaf, child_leaf
   type(item_node), pointer :: item1, item2

   ! Get the deepest branch
   current_branch => find_deepest_branch(mnatree)

   ! Create a new branch as child
   new_branch => add_new_branch(current_branch)

   parent_leaf => current_branch%first_branch_leaf
   do while (associated(parent_leaf))
      ! Check if this leaf needs processing
      item1 => parent_leaf%first_item1
      item2 => parent_leaf%first_item2
      if (.not. assigned .and. parent_leaf%num_items1 > 1) then
         assigned = .true.
         ! Create new leaf and move one item from each list to it
         child_leaf => add_new_leaf(new_branch, parent_leaf)
         call add_new_item1(child_leaf, item1%value)
         call add_new_item2(child_leaf, item2%value)
         item1 => item1%next_item
         item2 => item2%next_item
      end if
      child_leaf => add_new_leaf(new_branch, parent_leaf)
      do while (associated(item1))
         call add_new_item1(child_leaf, item1%value)
         call add_new_item2(child_leaf, item2%value)
         item1 => item1%next_item
         item2 => item2%next_item
      end do
      parent_leaf => parent_leaf%next_branch_leaf
   end do
end subroutine

subroutine separate_first_item(parent_leaf, new_branch)
   type(leaf_node), intent(inout) :: parent_leaf
   type(branch_node), intent(inout) :: new_branch
   type(leaf_node), pointer :: child_leaf
   type(item_node), pointer :: item1, item2

   child_leaf => add_new_leaf(new_branch, parent_leaf)
   call add_new_item1(child_leaf, parent_leaf%first_item1%value)
   call add_new_item2(child_leaf, parent_leaf%first_item2%value)
   child_leaf => add_new_leaf(new_branch, parent_leaf)
   item1 => parent_leaf%first_item1%next_item
   item2 => parent_leaf%first_item2%next_item
   do while (associated(item1))
      call add_new_item1(child_leaf, item1%value)
      call add_new_item2(child_leaf, item2%value)
      item1 => item1%next_item
      item2 => item2%next_item
   end do
end subroutine

subroutine copy_leaf_items(parent_leaf, new_branch)
   type(leaf_node), intent(inout) :: parent_leaf
   type(branch_node), intent(inout) :: new_branch
   type(leaf_node), pointer :: child_leaf
   type(item_node), pointer :: item1, item2

   child_leaf => add_new_leaf(new_branch, parent_leaf)
   item1 => parent_leaf%first_item1
   item2 => parent_leaf%first_item2
   do while (associated(item1))
      call add_new_item1(child_leaf, item1%value)
      call add_new_item2(child_leaf, item2%value)
      item1 => item1%next_item
      item2 => item2%next_item
   end do
end subroutine

subroutine split_mnas(mol1, mol2, mnatree)
   type(mol_type), intent(in) :: mol1, mol2
   type(branch_node), pointer, intent(inout) :: mnatree
   ! Local variables
   type(branch_node), pointer :: deepest_branch, new_branch
   type(leaf_node), pointer :: ileaf, jleaf
   logical :: assigned

   ! Find the deepest branch once at the start
   deepest_branch => find_deepest_branch(mnatree)

   ileaf => deepest_branch%first_branch_leaf
   do while (associated(ileaf))
      if (ileaf%num_items1 > 1) then
         ! Create a new branch for the current leaf
         new_branch => add_new_branch(deepest_branch)
         ! Process current leaf
         call separate_first_item(ileaf, new_branch)
         ! Process other leaves
         jleaf => deepest_branch%first_branch_leaf
         do while (associated(jleaf))
            if (.not. associated(ileaf, jleaf)) then
               call copy_leaf_items(jleaf, new_branch)
            end if
            jleaf => jleaf%next_branch_leaf
         end do
         ! Compute MNA types for the current branch
         do
            assigned = .false.
            call assign_next_item(new_branch, assigned)
            if (.not. assigned) exit
            call compute_consistent_mnas(mol1, mol2, new_branch)
         end do
      end if
      ileaf => ileaf%next_branch_leaf
   end do
end subroutine


