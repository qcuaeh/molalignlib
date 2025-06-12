module mna_recompute_arrays
use parameters
use molecule
use array_trees
implicit none

contains

subroutine resplit_part_mna_array(atoms1, atoms2, array_trees, part_idx, read_link_idx, write_link_idx)
! Array-based version of resplit_part_mna
   type(atom_type), dimension(:), intent(in) :: atoms1, atoms2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx, read_link_idx, write_link_idx
   ! Local variables
   integer :: item_idx, child_idx, target_part_idx
   integer, dimension(:), allocatable :: signature

   ! Reset last item pointers for all children
   child_idx = array_trees%parts(part_idx)%first_child_idx
   do while (child_idx > 0)
      array_trees%parts(child_idx)%last_item1_idx = 0
      array_trees%parts(child_idx)%last_item2_idx = 0
      child_idx = array_trees%parts(child_idx)%next_sibling_idx
   end do

   ! Process first molecule items
   item_idx = array_trees%parts(part_idx)%first_item1_idx
   do while (item_idx > 0)
      signature = get_signature_from_itemdir1_array(array_trees, read_link_idx, atoms1, array_trees%items1(item_idx)%value)
      target_part_idx = find_child_part_array(array_trees, part_idx, signature)
      if (target_part_idx == 0) error stop 'part not found'

      call update_item1_in_child_array(array_trees, target_part_idx, array_trees%items1(item_idx)%value)
      call update_itemdir1_entry_array(array_trees, write_link_idx, array_trees%items1(item_idx)%value, target_part_idx)

      item_idx = array_trees%items1(item_idx)%next_item_idx
   end do

   ! Process second molecule items
   item_idx = array_trees%parts(part_idx)%first_item2_idx
   do while (item_idx > 0)
      signature = get_signature_from_itemdir2_array(array_trees, read_link_idx, atoms2, array_trees%items2(item_idx)%value)
      target_part_idx = find_child_part_array(array_trees, part_idx, signature)
      if (target_part_idx == 0) error stop 'part not found'

      call update_item2_in_child_array(array_trees, target_part_idx, array_trees%items2(item_idx)%value)
      call update_itemdir2_entry_array(array_trees, write_link_idx, array_trees%items2(item_idx)%value, target_part_idx)

      item_idx = array_trees%items2(item_idx)%next_item_idx
   end do
end subroutine

function get_signature_from_itemdir1_array(array_trees, link_idx, atoms, item_value) result(signature)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: link_idx, item_value
   type(atom_type), dimension(:), intent(in) :: atoms
   integer, dimension(:), allocatable :: signature
   integer :: i, adjlist_size, itemdir1_start

   adjlist_size = size(atoms(item_value)%adjlist)
   allocate(signature(adjlist_size))
   itemdir1_start = array_trees%links(link_idx)%itemdir1_start_idx

   ! Direct array access pattern instead of complex pointer chasing
   do i = 1, adjlist_size
      signature(i) = array_trees%itemdir_entries(itemdir1_start + atoms(item_value)%adjlist(i) - 1)%part_idx
   end do
end function

function get_signature_from_itemdir2_array(array_trees, link_idx, atoms, item_value) result(signature)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: link_idx, item_value
   type(atom_type), dimension(:), intent(in) :: atoms
   integer, dimension(:), allocatable :: signature
   integer :: i, adjlist_size, itemdir2_start

   adjlist_size = size(atoms(item_value)%adjlist)
   allocate(signature(adjlist_size))
   itemdir2_start = array_trees%links(link_idx)%itemdir2_start_idx

   do i = 1, adjlist_size
      signature(i) = array_trees%itemdir_entries(itemdir2_start + atoms(item_value)%adjlist(i) - 1)%part_idx
   end do
end function

function find_child_part_array(array_trees, parent_idx, signature) result(child_idx)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: parent_idx
   integer, dimension(:), intent(in) :: signature
   integer :: child_idx

   child_idx = array_trees%parts(parent_idx)%first_child_idx
   do while (child_idx > 0)
      if (signature_equivalence_array(array_trees, child_idx, signature)) then
         return
      end if
      child_idx = array_trees%parts(child_idx)%next_sibling_idx
   end do

   child_idx = 0
end function

function signature_equivalence_array(array_trees, part_idx, target_signature) result(equiv)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   integer, dimension(:), intent(in) :: target_signature
   logical :: equiv
   integer :: matches1, matches2, i, j, sig_len

   sig_len = array_trees%parts(part_idx)%signature_length
   if (sig_len /= size(target_signature)) then
      equiv = .false.
      return
   end if

   ! Direct signature comparison - much simpler with fixed-size arrays!
   do i = 1, sig_len
      matches1 = 0
      matches2 = 0
      do j = 1, sig_len
         if (array_trees%parts(part_idx)%signature(i) == target_signature(j)) matches1 = matches1 + 1
         if (array_trees%parts(part_idx)%signature(i) == array_trees%parts(part_idx)%signature(j)) matches2 = matches2 + 1
      end do
      if (matches1 /= matches2) then
         equiv = .false.
         return
      end if
   end do

   equiv = .true.
end function

subroutine update_item1_in_child_array(array_trees, child_idx, item_value)
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: child_idx, item_value
   integer :: current_item_idx

   if (array_trees%parts(child_idx)%last_item1_idx == 0) then
      current_item_idx = array_trees%parts(child_idx)%first_item1_idx
   else
      current_item_idx = array_trees%items1(array_trees%parts(child_idx)%last_item1_idx)%next_item_idx
   end if

   array_trees%items1(current_item_idx)%value = item_value
   array_trees%parts(child_idx)%last_item1_idx = current_item_idx
end subroutine

subroutine update_item2_in_child_array(array_trees, child_idx, item_value)
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: child_idx, item_value
   integer :: current_item_idx

   if (array_trees%parts(child_idx)%last_item2_idx == 0) then
      current_item_idx = array_trees%parts(child_idx)%first_item2_idx
   else
      current_item_idx = array_trees%items2(array_trees%parts(child_idx)%last_item2_idx)%next_item_idx
   end if

   array_trees%items2(current_item_idx)%value = item_value
   array_trees%parts(child_idx)%last_item2_idx = current_item_idx
end subroutine

subroutine update_itemdir1_entry_array(array_trees, link_idx, item_value, part_idx)
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: link_idx, item_value, part_idx

   array_trees%itemdir_entries(array_trees%links(link_idx)%itemdir1_start_idx + item_value - 1)%part_idx = part_idx
end subroutine

subroutine update_itemdir2_entry_array(array_trees, link_idx, item_value, part_idx)
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: link_idx, item_value, part_idx

   array_trees%itemdir_entries(array_trees%links(link_idx)%itemdir2_start_idx + item_value - 1)%part_idx = part_idx
end subroutine

subroutine recompute_nextlevel_mnas_array(mol1, mol2, array_trees, link_idx)
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: link_idx
   ! Local variables
   integer :: partref_idx, next_link_idx

   next_link_idx = array_trees%links(link_idx)%next_link_idx

   ! Process all partrefs using direct array iteration instead of linked list traversal
   partref_idx = array_trees%links(link_idx)%first_partref_idx
   do while (partref_idx > 0)
      call resplit_part_mna_array(mol1%atoms, mol2%atoms, array_trees, &
                                  array_trees%partrefs(partref_idx)%part_idx, &
                                  link_idx, next_link_idx)
      partref_idx = array_trees%partrefs(partref_idx)%next_ref_idx
   end do
end subroutine

subroutine recompute_consistent_mnas_array(mol1, mol2, array_trees, branch_idx)
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: branch_idx
   integer :: link_idx

   ! Iterate through all links in branch using clean array access
   link_idx = array_trees%chains(branch_idx)%first_link_idx
   do while (link_idx > 0)
      call recompute_nextlevel_mnas_array(mol1, mol2, array_trees, link_idx)
      link_idx = array_trees%links(link_idx)%next_link_idx
   end do
end subroutine

subroutine resplit_part_first_array(array_trees, part_idx, write_link_idx)
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx, write_link_idx
   integer :: child_part1, child_part2, item1_idx, item2_idx

   ! Get child parts using direct array access
   child_part1 = array_trees%parts(part_idx)%first_child_idx
   child_part2 = array_trees%parts(child_part1)%next_sibling_idx

   ! Get first items
   item1_idx = array_trees%parts(part_idx)%first_item1_idx
   item2_idx = array_trees%parts(part_idx)%first_item2_idx

   ! Add first item to first child
   array_trees%items1(array_trees%parts(child_part1)%first_item1_idx)%value = array_trees%items1(item1_idx)%value
   array_trees%items2(array_trees%parts(child_part1)%first_item2_idx)%value = array_trees%items2(item2_idx)%value

   call update_itemdir1_entry_array(array_trees, write_link_idx, array_trees%items1(item1_idx)%value, child_part1)
   call update_itemdir2_entry_array(array_trees, write_link_idx, array_trees%items2(item2_idx)%value, child_part1)

   ! Add remaining items to second child
   array_trees%parts(child_part2)%last_item1_idx = 0
   array_trees%parts(child_part2)%last_item2_idx = 0

   item1_idx = array_trees%items1(item1_idx)%next_item_idx
   item2_idx = array_trees%items2(item2_idx)%next_item_idx

   do while (item1_idx > 0)
      call update_item1_in_child_array(array_trees, child_part2, array_trees%items1(item1_idx)%value)
      call update_item2_in_child_array(array_trees, child_part2, array_trees%items2(item2_idx)%value)
      call update_itemdir1_entry_array(array_trees, write_link_idx, array_trees%items1(item1_idx)%value, child_part2)
      call update_itemdir2_entry_array(array_trees, write_link_idx, array_trees%items2(item2_idx)%value, child_part2)

      item1_idx = array_trees%items1(item1_idx)%next_item_idx
      item2_idx = array_trees%items2(item2_idx)%next_item_idx
   end do
end subroutine

subroutine resplit_part_random_array(array_trees, part_idx, write_link_idx)
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx, write_link_idx
   integer :: child_part1, child_part2, item1_idx, item2_idx
   integer :: num_items1, num_items2, random_index1, random_index2, current_index
   integer :: chosen_item1_value, chosen_item2_value
   real :: random_real

   num_items1 = array_trees%parts(part_idx)%num_items1
   num_items2 = array_trees%parts(part_idx)%num_items2
   if (num_items1 < 1 .or. num_items2 < 1) error stop "Cannot split part with less than 1 item in either list"

   ! Generate random indices
   call random_number(random_real)
   random_index1 = int(random_real * num_items1) + 1
   call random_number(random_real)
   random_index2 = int(random_real * num_items2) + 1

   child_part1 = array_trees%parts(part_idx)%first_child_idx
   child_part2 = array_trees%parts(child_part1)%next_sibling_idx

   ! Find random items using optimized traversal
   item1_idx = array_trees%parts(part_idx)%first_item1_idx
   current_index = 1
   do while (current_index < random_index1)
      item1_idx = array_trees%items1(item1_idx)%next_item_idx
      current_index = current_index + 1
   end do
   chosen_item1_value = array_trees%items1(item1_idx)%value

   item2_idx = array_trees%parts(part_idx)%first_item2_idx
   current_index = 1
   do while (current_index < random_index2)
      item2_idx = array_trees%items2(item2_idx)%next_item_idx
      current_index = current_index + 1
   end do
   chosen_item2_value = array_trees%items2(item2_idx)%value

   write (stderr,*) part_idx, chosen_item1_value, chosen_item2_value

   ! Assign to first child
   array_trees%items1(array_trees%parts(child_part1)%first_item1_idx)%value = chosen_item1_value
   array_trees%items2(array_trees%parts(child_part1)%first_item2_idx)%value = chosen_item2_value
   call update_itemdir1_entry_array(array_trees, write_link_idx, chosen_item1_value, child_part1)
   call update_itemdir2_entry_array(array_trees, write_link_idx, chosen_item2_value, child_part1)

   ! Assign remaining items to second child using clean iteration pattern
   item1_idx = array_trees%parts(part_idx)%first_item1_idx
   current_index = 1
   array_trees%parts(child_part2)%last_item1_idx = 0

   do while (item1_idx > 0)
      if (current_index /= random_index1) then
         call update_item1_in_child_array(array_trees, child_part2, array_trees%items1(item1_idx)%value)
         call update_itemdir1_entry_array(array_trees, write_link_idx, array_trees%items1(item1_idx)%value, child_part2)
      end if
      item1_idx = array_trees%items1(item1_idx)%next_item_idx
      current_index = current_index + 1
   end do

   item2_idx = array_trees%parts(part_idx)%first_item2_idx
   current_index = 1
   array_trees%parts(child_part2)%last_item2_idx = 0

   do while (item2_idx > 0)
      if (current_index /= random_index2) then
         call update_item2_in_child_array(array_trees, child_part2, array_trees%items2(item2_idx)%value)
         call update_itemdir2_entry_array(array_trees, write_link_idx, array_trees%items2(item2_idx)%value, child_part2)
      end if
      item2_idx = array_trees%items2(item2_idx)%next_item_idx
      current_index = current_index + 1
   end do
end subroutine

recursive subroutine redistribute_items_array(mol1, mol2, array_trees, branch_idx)
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: branch_idx
   integer :: child_branch_idx, first_link_idx

   ! Process all child branches using clean iteration pattern
   child_branch_idx = array_trees%chains(branch_idx)%first_child_idx
   do while (child_branch_idx > 0)
      first_link_idx = array_trees%chains(child_branch_idx)%first_link_idx
      call resplit_part_random_array(array_trees, array_trees%chains(child_branch_idx)%split_part_idx, first_link_idx)
      call recompute_consistent_mnas_array(mol1, mol2, array_trees, child_branch_idx)
      call redistribute_items_array(mol1, mol2, array_trees, child_branch_idx)
      child_branch_idx = array_trees%chains(child_branch_idx)%next_sibling_idx
   end do
end subroutine

end module
