module assigntree_distribute_linked
use parameters
use random
use molecule
use lcrs_tree
implicit none
private
public distribute_items

contains

subroutine resplit_part_mna(vertices1, vertices2, itemdir1, itemdir2, part, link)
! Update item values of existing item nodes instead of adding new item nodes
   type(vertex_t), dimension(:), intent(in) :: vertices1, vertices2
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir1, itemdir2
   type(partree_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(item_node_t), pointer :: item
   type(part_nodeptr_t), dimension(:), allocatable :: signature
   type(partree_node_t), pointer :: child_part

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
      signature = itemdir1(vertices1(item%vertidx)%adjlist)
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
      child_part%last_item1%vertidx = item%vertidx
      link%itemdir1(item%vertidx)%ptr => child_part
      item => item%next_item
   end do

   ! Process second molecule items - update existing item nodes
   item => part%first_item2
   do while (associated(item))
      signature = itemdir2(vertices2(item%vertidx)%adjlist)
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
      child_part%last_item2%vertidx = item%vertidx
      link%itemdir2(item%vertidx)%ptr => child_part
      item => item%next_item
   end do
end subroutine

subroutine recompute_mna_partition(vertices1, vertices2, link)
! Compute next level MNA types - always keeps all children (original behavior)
   type(vertex_t), dimension(:), intent(in) :: vertices1, vertices2
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(partref_node_t), pointer :: partref

   ! Process all parts in the current partition
   partref => link%first_partref
   do while (associated(partref))
!      write (stderr,'(A,1X,A)') 'Part', address(partref%part)
      ! Distribute items based on signatures
      call resplit_part_mna(vertices1, vertices2, link%itemdir1, link%itemdir2, partref%part, link%next_link)
      partref => partref%nextref
   end do
end subroutine

subroutine distribute_part_items(vertices1, vertices2, branch)
! Iteratively compute MNA types until convergence
   type(vertex_t), dimension(:), intent(in) :: vertices1, vertices2
   type(assigntree_node_t), pointer, intent(inout) :: branch
   ! Local variables
   type(chain_node_t), pointer :: link
   integer :: link_idx, rand_idx1, rand_idx2

   rand_idx1 = random_uniform_integer(1, branch%split_part%num_items1)
   rand_idx2 = random_uniform_integer(1, branch%split_part%num_items2)
   call resplit_part_indexed(branch%split_part, branch%first_link, rand_idx1, rand_idx2)
!   call resplit_part_indexed(branch%split_part, branch%first_link, 1, 1)

   link_idx = 1
   link => branch%first_link
   do while (associated(link))
      ! Recompute next level MNAs
!      write (stderr,'(A,1X,I0)') 'Link', link_idx
!      call print_link_itemdir(link)
      call recompute_mna_partition(vertices1, vertices2, link)
      link_idx = link_idx + 1
      link => link%next_link
   end do
end subroutine

subroutine resplit_part_indexed(part, link, index1, index2)
   type(partree_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   integer, intent(in) :: index1, index2
   ! Local variables
   type(partree_node_t), pointer :: child_part1, child_part2
   type(item_node_t), pointer :: item1, item2
   integer :: current_index

   ! Validate indices
   if (index1 < 1 .or. index1 > part%num_items1) error stop "index1 out of range"
   if (index2 < 1 .or. index2 > part%num_items2) error stop "index2 out of range"

   ! Get pointers to first and second child parts
   child_part1 => part%first_child_part
   child_part2 => part%first_child_part%next_sibling_part

   ! Find and assign the item at index1 from items1 to first child
   item1 => part%first_item1
   current_index = 1
   do while (current_index < index1)
      item1 => item1%next_item
      current_index = current_index + 1
   end do
   child_part1%first_item1%vertidx = item1%vertidx
   link%itemdir1(item1%vertidx)%ptr => child_part1

   ! Find and assign the item at index2 from items2 to first child
   item2 => part%first_item2
   current_index = 1
   do while (current_index < index2)
      item2 => item2%next_item
      current_index = current_index + 1
   end do
   child_part1%first_item2%vertidx = item2%vertidx
   link%itemdir2(item2%vertidx)%ptr => child_part1

   ! Now assign all other items from items1 to second child
   item1 => part%first_item1
   current_index = 1
   child_part2%last_item1 => null()

   do while (associated(item1))
      if (current_index /= index1) then
         ! Move to next position in second child
         if (.not. associated(child_part2%last_item1)) then
            child_part2%last_item1 => child_part2%first_item1
         else
            child_part2%last_item1 => child_part2%last_item1%next_item
         end if
         ! This is not the selected item, assign to second child
         child_part2%last_item1%vertidx = item1%vertidx
         link%itemdir1(item1%vertidx)%ptr => child_part2
      end if

      ! Move to next item in parent
      item1 => item1%next_item
      current_index = current_index + 1
   end do

   ! Now assign all other items from items2 to second child
   item2 => part%first_item2
   current_index = 1
   child_part2%last_item2 => null()

   do while (associated(item2))
      if (current_index /= index2) then
         ! Move to next position in second child
         if (.not. associated(child_part2%last_item2)) then
            child_part2%last_item2 => child_part2%first_item2
         else
            child_part2%last_item2 => child_part2%last_item2%next_item
         end if
         ! This is not the selected item, assign to second child
         child_part2%last_item2%vertidx = item2%vertidx
         link%itemdir2(item2%vertidx)%ptr => child_part2
      end if

      ! Move to next item in parent
      item2 => item2%next_item
      current_index = current_index + 1
   end do
end subroutine

recursive subroutine distribute_items(vertices1, vertices2, branch)
   type(vertex_t), dimension(:), intent(in) :: vertices1, vertices2
   type(assigntree_node_t), pointer, intent(inout) :: branch
   type(assigntree_node_t), pointer :: child_branch
   type(partree_node_t), pointer :: child_part

!   call random_init(.true., .true.)

   ! Process all children of this branch
   child_branch => branch%first_child_chain
   do while (associated(child_branch))
      ! Process this child branch's split_part
!      write (stderr,*)
!      write (stderr,'(A,1X,A)') 'Split Part', address(child_branch%split_part)
      ! Add target part children to links
      child_part => child_branch%split_part%first_child_part
      call distribute_part_items(vertices1, vertices2, child_branch)

      ! Recursively process this child's descendants (depth-first)
      call distribute_items(vertices1, vertices2, child_branch)

      ! Move to next sibling
      child_branch => child_branch%next_sibling_chain
   end do
end subroutine

end module
