module mna_recompute
use parameters
use molecule
use lcrs_tree
implicit none

contains

subroutine resplit_part_mna(atoms1, atoms2, itemdir1, itemdir2, part)
! Update item values of existing item nodes instead of adding new item nodes
   type(atom_type), dimension(:), intent(in) :: atoms1, atoms2
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir1, itemdir2
   type(part_node_t), pointer, intent(inout) :: part
   ! Local variables
   type(item_node_t), pointer :: item
   type(part_nodeptr_t), dimension(:), allocatable :: signature
   type(part_node_t), pointer :: child_part

   ! Reset last item pointers for all children
   child_part => part%first_child_part
   do while (associated(child_part))
      child_part%last_item1 => null()
      child_part%last_item2 => null()
      child_part => child_part%next_sibling_part
   end do

   ! Process first molecule items - update existing item nodes
   item => part%first_item1
   do while (associated(item))
      signature = itemdir1(atoms1(item%value)%adjlist)
      child_part => find_child_part(part, signature)
      if (.not. associated(child_part)) then
         call print_part_signature(signature)
         error stop 'part not found'
      end if
      if (.not. associated(child_part%last_item1)) then
         child_part%last_item1 => child_part%first_item1
      else
         child_part%last_item1 => child_part%last_item1%next_item
      end if
      child_part%last_item1%value = item%value
      item => item%next_item
   end do

   ! Process second molecule items - update existing item nodes
   item => part%first_item2
   do while (associated(item))
      signature = itemdir2(atoms2(item%value)%adjlist)
      child_part => find_child_part(part, signature)
      if (.not. associated(child_part)) then
         call print_part_signature(signature)
         error stop 'part not found'
      end if
      if (.not. associated(child_part%last_item2)) then
         child_part%last_item2 => child_part%first_item2
      else
         child_part%last_item2 => child_part%last_item2%next_item
      end if
      child_part%last_item2%value = item%value
      item => item%next_item
   end do
end subroutine

subroutine recompute_nextlevel_mnas(mol1, mol2, link)
! Compute next level MNA types - always keeps all children (original behavior)
   type(mol_type), intent(in) :: mol1, mol2
   type(link_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(partref_node_t), pointer :: partref
   type(part_node_t), pointer :: child_part

   ! Process all parts in the current partition
   partref => link%first_partref
   do while (associated(partref))
!      write (stderr,'(A,1X,A)') 'Part', address(partref%part)
      ! Distribute items based on signatures
      call resplit_part_mna(mol1%atoms, mol2%atoms, link%itemdir1, link%itemdir2, partref%part)
      child_part => partref%part%first_child_part
      do while (associated(child_part))
         call update_itemdir(link%next_link, child_part)
         child_part => child_part%next_sibling_part
      end do
      partref => partref%nextref
   end do
end subroutine

subroutine recompute_consistent_mnas(mol1, mol2, branch)
! Iteratively compute MNA types until convergence
   type(mol_type), intent(in) :: mol1, mol2
   type(split_node_t), pointer, intent(inout) :: branch
   ! Local variables
   type(link_node_t), pointer :: link
   integer link_idx

   link_idx = 1
   link => branch%first_link
   do while (associated(link))
      ! Recompute next level MNAs
!      write (stderr,'(A,1X,I0)') 'Link', link_idx
!      call print_link_itemdir(link)
      call recompute_nextlevel_mnas(mol1, mol2, link)
      link_idx = link_idx + 1
      link => link%next_link
   end do
end subroutine

subroutine resplit_part_first(part)
   type(part_node_t), pointer, intent(inout) :: part
   ! Local variables
   type(part_node_t), pointer :: child_part
   type(item_node_t), pointer :: item1, item2

   ! Add first item to first child
   child_part => part%first_child_part
   child_part%first_item1%value = part%first_item1%value
   child_part%first_item2%value = part%first_item2%value

   ! Add second item to second child
   child_part => part%first_child_part%next_sibling_part
   child_part%first_item1%value = part%first_item1%next_item%value
   child_part%first_item2%value = part%first_item2%next_item%value

   ! Add remaining items to second child
   item1 => part%first_item1%next_item%next_item
   item2 => part%first_item2%next_item%next_item
   child_part%last_item1 => child_part%first_item1
   child_part%last_item2 => child_part%first_item2
   do while (associated(item1))
      child_part%last_item1 => child_part%last_item1%next_item
      child_part%last_item2 => child_part%last_item2%next_item
      child_part%last_item1%value = item1%value
      child_part%last_item2%value = item2%value
      item1 => item1%next_item
      item2 => item2%next_item
   end do
end subroutine

subroutine resplit_part_random(part)
   type(part_node_t), pointer, intent(inout) :: part
   ! Local variables
   type(part_node_t), pointer :: child_part1, child_part2
   type(item_node_t), pointer :: item1, item2
   integer :: num_items1, num_items2, random_index1, random_index2, current_index
   real :: random_real

   ! Count the number of items in each list
   num_items1 = part%num_items1
   num_items2 = part%num_items2
   if (num_items1 < 1 .or. num_items2 < 1) error stop "Cannot split part with less than 1 item in either list"

   ! Generate random index for items1 (1 to num_items1)
   call random_number(random_real)
   random_index1 = int(random_real * num_items1) + 1

   ! Generate random index for items2 (1 to num_items2)
   call random_number(random_real)
   random_index2 = int(random_real * num_items2) + 1

   ! Get pointers to first and second child parts
   child_part1 => part%first_child_part
   child_part2 => part%first_child_part%next_sibling_part

   ! Find and assign the random item from items1 to first child
   item1 => part%first_item1
   current_index = 1
   do while (current_index < random_index1)
      item1 => item1%next_item
      current_index = current_index + 1
   end do
   child_part1%first_item1%value = item1%value

   ! Find and assign the random item from items2 to first child
   item2 => part%first_item2
   current_index = 1
   do while (current_index < random_index2)
      item2 => item2%next_item
      current_index = current_index + 1
   end do
   child_part1%first_item2%value = item2%value

   ! Now assign all other items from items1 to second child
   item1 => part%first_item1
   current_index = 1
   child_part2%last_item1 => child_part2%first_item1

   do while (associated(item1))
      if (current_index /= random_index1) then
         ! This is not the random item, assign to second child
         child_part2%last_item1%value = item1%value
         
         ! Move to next position in second child (if not the last item)
         if (associated(child_part2%last_item1%next_item)) then
            child_part2%last_item1 => child_part2%last_item1%next_item
         end if
      end if
      
      ! Move to next item in parent
      item1 => item1%next_item
      current_index = current_index + 1
   end do

   ! Now assign all other items from items2 to second child
   item2 => part%first_item2
   current_index = 1
   child_part2%last_item2 => child_part2%first_item2

   do while (associated(item2))
      if (current_index /= random_index2) then
         ! This is not the random item, assign to second child
         child_part2%last_item2%value = item2%value
         
         ! Move to next position in second child (if not the last item)
         if (associated(child_part2%last_item2%next_item)) then
            child_part2%last_item2 => child_part2%last_item2%next_item
         end if
      end if
      
      ! Move to next item in parent
      item2 => item2%next_item
      current_index = current_index + 1
   end do
end subroutine

recursive subroutine redistribute_items(mol1, mol2, branch)
   type(mol_type), intent(in) :: mol1, mol2
   type(split_node_t), pointer, intent(inout) :: branch
   type(split_node_t), pointer :: child_branch
   type(part_node_t), pointer :: child_part

   ! Process all children of this branch
   child_branch => branch%first_child_branch
   do while (associated(child_branch))
      ! Process this child branch's split_part
!      write (stderr,*)
!      write (stderr,'(A,1X,A)') 'Split Part', address(child_branch%split_part)
!      call resplit_part_first(child_branch%split_part)
      call resplit_part_random(child_branch%split_part)
      ! Add target part children to links
      child_part => child_branch%split_part%first_child_part
      do while (associated(child_part))
         call update_itemdir(child_branch%first_link, child_part)
         child_part => child_part%next_sibling_part
      end do
      call recompute_consistent_mnas(mol1, mol2, child_branch)
      
      ! Recursively process this child's descendants (depth-first)
      call redistribute_items(mol1, mol2, child_branch)
      
      ! Move to next sibling
      child_branch => child_branch%next_sibling_branch
   end do
end subroutine

end module
