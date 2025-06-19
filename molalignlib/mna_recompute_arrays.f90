module mna_recompute_arrays
use parameters
use molecule
use array_trees
implicit none

! Module-level signature workspace to eliminate allocations
integer :: signature_array(MAX_COORD)
integer :: signature_length

! Module-level redistribution counters (replaces items1_fill_count/items2_fill_count in parts)
! These track the current number of items placed in each part during redistribution
integer, parameter :: MAX_PARTS = 10000
integer :: part_items1_fill(MAX_PARTS)
integer :: part_items2_fill(MAX_PARTS)

! Module-level total squared distance tracking variables
real(rk) :: total_squared_distance
integer :: total_assigned_pairs

contains

function signature_equivalence_array(array_trees, part_idx) result(equiv)
   ! OPTIMIZED: Fast path for length-1 signatures (most common case)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   logical :: equiv
   integer :: signature_frequencies, i, j

   if (signature_length /= array_trees%parts(part_idx)%signature_length) then
      equiv = .false.
      return
   end if

   ! FAST PATH: Direct comparison for length-1 signatures (most common)
   if (array_trees%parts(part_idx)%signature_length == 1) then
      equiv = (signature_array(1) == array_trees%parts(part_idx)%signature_values(1))
      return
   end if

   ! GENERIC PATH
   do i = 1, array_trees%parts(part_idx)%signature_unique_count
      signature_frequencies = 0

      ! Count matches in target signature
      do j = 1, signature_length
         if (signature_array(j) == array_trees%parts(part_idx)%signature_values(i)) then
            signature_frequencies = signature_frequencies + 1
         end if
      end do

      if (signature_frequencies /= array_trees%parts(part_idx)%signature_frequencies(i)) then
         equiv = .false.
         return
      end if
   end do

   equiv = .true.
end function

function find_child_part_array(array_trees, parent_idx) result(child_idx)
   ! OPTIMIZED: Fast path for exactly 2 children (most common case)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: parent_idx
   integer :: child_idx, i, num_children

   num_children = array_trees%parts(parent_idx)%num_children

   ! FAST PATH: Exactly 2 children (most common case)
   if (num_children == 2) then
      ! Check first child
      child_idx = array_trees%parts(parent_idx)%child_indices(1)
      if (signature_equivalence_array(array_trees, child_idx)) then
         return
      end if

      ! Check second child
      child_idx = array_trees%parts(parent_idx)%child_indices(2)
      if (signature_equivalence_array(array_trees, child_idx)) then
         return
      end if

      child_idx = 0
      return
   end if

   ! GENERIC PATH
   do i = 1, num_children
      child_idx = array_trees%parts(parent_idx)%child_indices(i)
      if (signature_equivalence_array(array_trees, child_idx)) then
         return
      end if
   end do

   child_idx = 0
end function

subroutine calculate_leaf_squared_distance_contribution(mol1, mol2, array_trees, part_idx)
   ! Calculate squared distance contribution for a leaf part (items1_count == items2_count == 1)
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   integer :: item1_idx, item2_idx
   real(rk) :: squared_distance
   real(rk) :: dx, dy, dz

   ! Verify this is a leaf part with exactly one item from each molecule
   if (array_trees%parts(part_idx)%items1_count /= 1 .or. &
       array_trees%parts(part_idx)%items2_count /= 1) then
      write(stderr, '(A,I0,A,I0,A,I0)') 'Warning: Part ', part_idx, &
         ' is not a proper leaf (items1=', array_trees%parts(part_idx)%items1_count, &
         ', items2=', array_trees%parts(part_idx)%items2_count, ')'
      return
   end if

   ! Get the assigned items
   item1_idx = array_trees%item1_values(array_trees%parts(part_idx)%items1_offset + 1)
   item2_idx = array_trees%item2_values(array_trees%parts(part_idx)%items2_offset + 1)

   ! Calculate squared distance between assigned atoms
   dx = mol1%atoms(item1_idx)%coords(1) - mol2%atoms(item2_idx)%coords(1)
   dy = mol1%atoms(item1_idx)%coords(2) - mol2%atoms(item2_idx)%coords(2)
   dz = mol1%atoms(item1_idx)%coords(3) - mol2%atoms(item2_idx)%coords(3)
   
   squared_distance = dx*dx + dy*dy + dz*dz

   ! Update progressive total squared distance calculation
   total_squared_distance = total_squared_distance + squared_distance
   total_assigned_pairs = total_assigned_pairs + 1

   ! Optional: Print progressive total squared distance for monitoring
!   write(stderr, '(A,I0,A,I0,A,F8.4)') &
!      'Assignment: atom ', item1_idx, ' -> atom ', item2_idx, &
!      ', pair distance=', squared_distance
end subroutine

subroutine check_leaf_parts_for_squared_distance(mol1, mol2, array_trees, part_idx)
   ! Check if child parts are leaves and calculate squared distance contributions
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   integer :: i, child_idx

   ! Check all children of this part
   do i = 1, array_trees%parts(part_idx)%num_children
      child_idx = array_trees%parts(part_idx)%child_indices(i)
      
      ! If this child is a leaf (no children), calculate squared distance contribution
      if (array_trees%parts(child_idx)%num_children == 0) then
         call calculate_leaf_squared_distance_contribution(mol1, mol2, array_trees, child_idx)
      end if
   end do
end subroutine

subroutine resplit_part_mna_array(mol1, mol2, array_trees, part_idx, read_link_idx, write_link_idx)
! Array-based version of resplit_part_mna with inlined signature generation
! UPDATED: Now uses module-level redistribution counters instead of part-level fill counts
! UPDATED: Added RMSD calculation for newly created leaf parts
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx, read_link_idx, write_link_idx
   ! Local variables
   integer :: i, j, target_part_idx, item_value, target_idx, part_ref_idx
   integer :: items1_offset, items1_count, items2_offset, items2_count
   integer :: itemdir1_offset, itemdir2_offset, read_itemdir1_offset, read_itemdir2_offset
   integer :: child_idx

   ! INITIALIZATION: Reset redistribution counters for all child parts
   do i = 1, array_trees%parts(part_idx)%num_children
      child_idx = array_trees%parts(part_idx)%child_indices(i)
      part_items1_fill(child_idx) = 0
      part_items2_fill(child_idx) = 0
   end do

   ! Extract commonly used values for readability
   items1_offset = array_trees%parts(part_idx)%items1_offset
   items1_count = array_trees%parts(part_idx)%items1_count
   items2_offset = array_trees%parts(part_idx)%items2_offset
   items2_count = array_trees%parts(part_idx)%items2_count
   itemdir1_offset = array_trees%links(write_link_idx)%itemdir1_offset
   itemdir2_offset = array_trees%links(write_link_idx)%itemdir2_offset
   read_itemdir1_offset = array_trees%links(read_link_idx)%itemdir1_offset
   read_itemdir2_offset = array_trees%links(read_link_idx)%itemdir2_offset

   ! Process first molecule items with inlined signature generation
   do i = 1, items1_count
      item_value = array_trees%item1_values(items1_offset + i)

      ! INLINED: Generate compact signature from itemdir1
      signature_length = 0
      do j = 1, size(mol1%atoms(item_value)%adjlist)
         part_ref_idx = array_trees%itemdir_entries(read_itemdir1_offset + mol1%atoms(item_value)%adjlist(j))
         if (part_ref_idx /= 0) then
            signature_length = signature_length + 1
            signature_array(signature_length) = part_ref_idx
         end if
      end do

      target_part_idx = find_child_part_array(array_trees, part_idx)
      if (target_part_idx == 0) error stop 'part not found'

      ! Add item to target child using module-level redistribution counter
      part_items1_fill(target_part_idx) = part_items1_fill(target_part_idx) + 1
      target_idx = array_trees%parts(target_part_idx)%items1_offset + part_items1_fill(target_part_idx)
      array_trees%item1_values(target_idx) = item_value

      ! Update itemdir
      array_trees%itemdir_entries(itemdir1_offset + item_value) = target_part_idx
   end do

   ! Process second molecule items with inlined signature generation
   do i = 1, items2_count
      item_value = array_trees%item2_values(items2_offset + i)

      ! INLINED: Generate compact signature from itemdir2
      signature_length = 0
      do j = 1, size(mol2%atoms(item_value)%adjlist)
         part_ref_idx = array_trees%itemdir_entries(read_itemdir2_offset + mol2%atoms(item_value)%adjlist(j))
         if (part_ref_idx /= 0) then
            signature_length = signature_length + 1
            signature_array(signature_length) = part_ref_idx
         end if
      end do

      target_part_idx = find_child_part_array(array_trees, part_idx)
      if (target_part_idx == 0) error stop 'part not found'

      ! Add item to target child using module-level redistribution counter
      part_items2_fill(target_part_idx) = part_items2_fill(target_part_idx) + 1
      target_idx = array_trees%parts(target_part_idx)%items2_offset + part_items2_fill(target_part_idx)
      array_trees%item2_values(target_idx) = item_value

      ! Update itemdir
      array_trees%itemdir_entries(itemdir2_offset + item_value) = target_part_idx
   end do

   ! Check if the split created any leaf parts and calculate squared distance contributions
   call check_leaf_parts_for_squared_distance(mol1, mol2, array_trees, part_idx)
end subroutine

subroutine recompute_nextlevel_mnas_array(mol1, mol2, array_trees, link_idx)
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: link_idx
   integer :: next_link_idx, i, part_idx
   integer :: num_parts, partref_offset

   next_link_idx = link_idx + 1
   num_parts = array_trees%links(link_idx)%num_parts
   partref_offset = array_trees%links(link_idx)%partref_offset

   do i = 1, num_parts
      part_idx = array_trees%partref_entries(partref_offset + i)
      call resplit_part_mna_array(mol1, mol2, array_trees, part_idx, link_idx, next_link_idx)
   end do
end subroutine

subroutine recompute_consistent_mnas_array(mol1, mol2, array_trees, branch_idx)
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: branch_idx
   integer :: i, link_idx
   integer :: num_links, link_offset

   num_links = array_trees%chains(branch_idx)%num_links
   link_offset = array_trees%chains(branch_idx)%link_offset

   do i = 1, num_links
      link_idx = link_offset + i
      call recompute_nextlevel_mnas_array(mol1, mol2, array_trees, link_idx)
   end do
end subroutine

subroutine resplit_part_first_array(mol1, mol2, array_trees, part_idx, write_link_idx)
   ! OPTIMIZED: No redistribution counters needed - direct placement only
   ! UPDATED: Added RMSD calculation for newly created leaf parts
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx, write_link_idx
   integer :: child_part1, child_part2, first_item1, first_item2, i, item_value, target_idx
   integer :: items1_offset, items1_count, items2_offset, items2_count
   integer :: itemdir1_offset, itemdir2_offset

   ! Extract commonly used offsets and values
   items1_offset = array_trees%parts(part_idx)%items1_offset
   items1_count = array_trees%parts(part_idx)%items1_count
   items2_offset = array_trees%parts(part_idx)%items2_offset
   items2_count = array_trees%parts(part_idx)%items2_count
   itemdir1_offset = array_trees%links(write_link_idx)%itemdir1_offset
   itemdir2_offset = array_trees%links(write_link_idx)%itemdir2_offset

   ! Get child parts using direct array access
   child_part1 = array_trees%parts(part_idx)%child_indices(1)
   child_part2 = array_trees%parts(part_idx)%child_indices(2)

   ! Get first items
   first_item1 = array_trees%item1_values(items1_offset + 1)
   first_item2 = array_trees%item2_values(items2_offset + 1)

   ! Add first items to first child (direct placement)
   array_trees%item1_values(array_trees%parts(child_part1)%items1_offset + 1) = first_item1
   array_trees%item2_values(array_trees%parts(child_part1)%items2_offset + 1) = first_item2
   array_trees%itemdir_entries(itemdir1_offset + first_item1) = child_part1
   array_trees%itemdir_entries(itemdir2_offset + first_item2) = child_part1

   ! Add remaining items to second child (direct placement - no redistribution counters needed)
   do i = 2, items1_count
      item_value = array_trees%item1_values(items1_offset + i)
      target_idx = array_trees%parts(child_part2)%items1_offset + (i - 1)
      array_trees%item1_values(target_idx) = item_value
      array_trees%itemdir_entries(itemdir1_offset + item_value) = child_part2
   end do

   do i = 2, items2_count
      item_value = array_trees%item2_values(items2_offset + i)
      target_idx = array_trees%parts(child_part2)%items2_offset + (i - 1)
      array_trees%item2_values(target_idx) = item_value
      array_trees%itemdir_entries(itemdir2_offset + item_value) = child_part2
   end do

   ! Check if the split created any leaf parts and calculate squared distance contributions
   call check_leaf_parts_for_squared_distance(mol1, mol2, array_trees, part_idx)
end subroutine

subroutine resplit_part_random_array(mol1, mol2, array_trees, part_idx, write_link_idx)
   ! OPTIMIZED: No redistribution counters needed - direct placement only
   ! UPDATED: Added RMSD calculation for newly created leaf parts
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx, write_link_idx
   integer :: child_part1, child_part2, i, item_value, target_idx
   integer :: num_items1, num_items2, random_index1, random_index2
   integer :: chosen_item1, chosen_item2
   integer :: items1_offset, items2_offset, itemdir1_offset, itemdir2_offset
   real :: random_real

   ! Extract commonly used values
   items1_offset = array_trees%parts(part_idx)%items1_offset
   items2_offset = array_trees%parts(part_idx)%items2_offset
   num_items1 = array_trees%parts(part_idx)%items1_count
   num_items2 = array_trees%parts(part_idx)%items2_count
   itemdir1_offset = array_trees%links(write_link_idx)%itemdir1_offset
   itemdir2_offset = array_trees%links(write_link_idx)%itemdir2_offset

   if (num_items1 < 1 .or. num_items2 < 1) error stop "Cannot split part with less than 1 item in either list"

   ! Generate random indices
   call random_number(random_real)
   random_index1 = int(random_real * num_items1) + 1
   call random_number(random_real)
   random_index2 = int(random_real * num_items2) + 1

   ! Get child parts using direct array access
   child_part1 = array_trees%parts(part_idx)%child_indices(1)
   child_part2 = array_trees%parts(part_idx)%child_indices(2)

   ! Get randomly chosen items
   chosen_item1 = array_trees%item1_values(items1_offset + random_index1)
   chosen_item2 = array_trees%item2_values(items2_offset + random_index2)

   ! Assign chosen items to first child (direct placement)
   array_trees%item1_values(array_trees%parts(child_part1)%items1_offset + 1) = chosen_item1
   array_trees%item2_values(array_trees%parts(child_part1)%items2_offset + 1) = chosen_item2
   array_trees%itemdir_entries(itemdir1_offset + chosen_item1) = child_part1
   array_trees%itemdir_entries(itemdir2_offset + chosen_item2) = child_part1

   ! Copy remaining items1 to second child (direct placement - no redistribution counters needed)
   target_idx = array_trees%parts(child_part2)%items1_offset
   do i = 1, num_items1
      if (i /= random_index1) then
         item_value = array_trees%item1_values(items1_offset + i)
         target_idx = target_idx + 1
         array_trees%item1_values(target_idx) = item_value
         array_trees%itemdir_entries(itemdir1_offset + item_value) = child_part2
      end if
   end do

   ! Copy remaining items2 to second child (direct placement - no redistribution counters needed)
   target_idx = array_trees%parts(child_part2)%items2_offset
   do i = 1, num_items2
      if (i /= random_index2) then
         item_value = array_trees%item2_values(items2_offset + i)
         target_idx = target_idx + 1
         array_trees%item2_values(target_idx) = item_value
         array_trees%itemdir_entries(itemdir2_offset + item_value) = child_part2
      end if
   end do

   ! Check if the split created any leaf parts and calculate squared distance contributions
   call check_leaf_parts_for_squared_distance(mol1, mol2, array_trees, part_idx)
end subroutine

recursive subroutine redistribute_items_array_recursive(mol1, mol2, array_trees, branch_idx)
   ! UPDATED: No longer performs global initialization - this is now done in wrapper
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: branch_idx
   integer :: child_branch_idx, first_link_idx, split_part_idx, i

   ! OPTIMIZATION: Use direct array access instead of linked list traversal
   ! This provides much better cache locality and eliminates pointer chasing
   do i = 1, array_trees%chains(branch_idx)%num_children
      child_branch_idx = array_trees%chains(branch_idx)%child_indices(i)

      first_link_idx = array_trees%chains(child_branch_idx)%link_offset + 1
      split_part_idx = array_trees%chains(child_branch_idx)%split_part_idx

!      call resplit_part_first_array(mol1, mol2, array_trees, split_part_idx, first_link_idx)
      call resplit_part_random_array(mol1, mol2, array_trees, split_part_idx, first_link_idx)
      call recompute_consistent_mnas_array(mol1, mol2, array_trees, child_branch_idx)
      call redistribute_items_array_recursive(mol1, mol2, array_trees, child_branch_idx)
   end do
end subroutine

subroutine redistribute_items_array(mol1, mol2, array_trees, branch_idx, final_total_squared_distance)
   ! OPTIMIZED: Simplified wrapper - only resets itemdir since redistribution counters are handled locally
   ! UPDATED: Added total squared distance calculation and returns final value
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: branch_idx
   real(rk), intent(out), optional :: final_total_squared_distance

   ! Initialize total squared distance tracking variables
   total_squared_distance = 0.0_rk
   total_assigned_pairs = 0

   ! Clear all itemdir entries
   array_trees%itemdir_entries = 0

   write(stderr, '(A)') "=== Starting item redistribution with squared distance tracking ==="

   ! Perform redistribution
   call redistribute_items_array_recursive(mol1, mol2, array_trees, branch_idx)

   ! Report final total squared distance
   if (total_assigned_pairs > 0) then
      write(stderr, '(A)') repeat("=", 50)
      write(stderr, '(A,I0)') "Total assigned pairs: ", total_assigned_pairs
      write(stderr, '(A,F10.4)') "Total squared distance: ", total_squared_distance
      write(stderr, '(A)') repeat("=", 50)
   else
      write(stderr, '(A)') "Warning: No leaf assignments found!"
      total_squared_distance = 0.0_rk
   end if

   ! Return final total squared distance if requested
   if (present(final_total_squared_distance)) then
      final_total_squared_distance = total_squared_distance
   end if
end subroutine

end module
