module mna_recompute_arrays
use parameters
use molecule
use array_trees
implicit none

! NOTE: This version assumes offset-based fields in array_trees data structures:
! 
! In part_array_t:
!   - items1_offset instead of items1_start_idx (where offset = start_idx - 1)
!   - items2_offset instead of items2_start_idx (where offset = start_idx - 1)
!
! In link_array_t:
!   - itemdir1_offset instead of itemdir1_start_idx (where offset = start_idx - 1)
!   - itemdir2_offset instead of itemdir2_start_idx (where offset = start_idx - 1)
!   - partref_offset instead of partref_start_idx (where offset = start_idx - 1)
!
! In chain_array_t:
!   - link_offset instead of link_start_idx (where offset = start_idx - 1)
!
! This version uses intermediate variables for:
! - Shorter, more declarative code lines
! - Better readability and maintainability
! - Clearer separation of concerns
!
! OPTIMIZATION: Uses compact translation arrays to eliminate zeros from signatures

! Compact mapping for efficient signature handling
type :: compact_map_t
   integer, allocatable :: part_to_compact(:)   ! maps original part_idx to compact index (0 if not used)
   integer, allocatable :: compact_parts(:)     ! array of non-zero part indices in compact order
   integer :: num_compact                       ! number of non-zero part indices
end type

contains

subroutine resplit_part_mna_array(atoms1, atoms2, array_trees, part_idx, read_link_idx, write_link_idx)
! Array-based version of resplit_part_mna
   type(atom_type), dimension(:), intent(in) :: atoms1, atoms2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: part_idx, read_link_idx, write_link_idx
   ! Local variables
   integer :: i, child_idx, target_part_idx, item_value, target_idx
   integer :: items1_offset, items1_count, items2_offset, items2_count
   integer :: itemdir1_offset, itemdir2_offset
   integer, dimension(:), allocatable :: compact_signature
   type(compact_map_t) :: read_link_compact_map

   ! Build compact map for the link we're reading from
   call build_compact_map_for_link(array_trees, read_link_idx, read_link_compact_map)

   ! Extract commonly used values for readability
   items1_offset = array_trees%parts(part_idx)%items1_offset
   items1_count = array_trees%parts(part_idx)%items1_count
   items2_offset = array_trees%parts(part_idx)%items2_offset
   items2_count = array_trees%parts(part_idx)%items2_count
   itemdir1_offset = array_trees%links(write_link_idx)%itemdir1_offset
   itemdir2_offset = array_trees%links(write_link_idx)%itemdir2_offset

   ! Reset fill counters for all children
   child_idx = array_trees%parts(part_idx)%first_child_idx
   do while (child_idx > 0)
      array_trees%parts(child_idx)%items1_fill_count = 0
      array_trees%parts(child_idx)%items2_fill_count = 0
      child_idx = array_trees%parts(child_idx)%next_sibling_idx
   end do

   ! Process first molecule items
   do i = 1, items1_count
      item_value = array_trees%item1_values(items1_offset + i)
      compact_signature = get_compact_signature_from_itemdir1_array(array_trees, read_link_idx, atoms1, item_value, &
            read_link_compact_map)
      target_part_idx = find_child_part_compact_array(array_trees, part_idx, compact_signature, read_link_compact_map)
      if (target_part_idx == 0) error stop 'part not found'

      ! Add item to target child
      array_trees%parts(target_part_idx)%items1_fill_count = array_trees%parts(target_part_idx)%items1_fill_count + 1
      target_idx = array_trees%parts(target_part_idx)%items1_offset + array_trees%parts(target_part_idx)%items1_fill_count
      array_trees%item1_values(target_idx) = item_value
      
      ! Update itemdir
      array_trees%itemdir_entries(itemdir1_offset + item_value) = target_part_idx
   end do

   ! Process second molecule items
   do i = 1, items2_count
      item_value = array_trees%item2_values(items2_offset + i)
      compact_signature = get_compact_signature_from_itemdir2_array(array_trees, read_link_idx, atoms2, item_value, &
            read_link_compact_map)
      target_part_idx = find_child_part_compact_array(array_trees, part_idx, compact_signature, read_link_compact_map)
      if (target_part_idx == 0) error stop 'part not found'

      ! Add item to target child
      array_trees%parts(target_part_idx)%items2_fill_count = array_trees%parts(target_part_idx)%items2_fill_count + 1
      target_idx = array_trees%parts(target_part_idx)%items2_offset + array_trees%parts(target_part_idx)%items2_fill_count
      array_trees%item2_values(target_idx) = item_value
      
      ! Update itemdir
      array_trees%itemdir_entries(itemdir2_offset + item_value) = target_part_idx
   end do

   ! Build compact map for the link we just wrote to (for future reads)
   call build_compact_map_for_link_itemdir(array_trees, write_link_idx)
end subroutine

subroutine build_compact_map_for_link(array_trees, link_idx, compact_map)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: link_idx
   type(compact_map_t), intent(out) :: compact_map
   integer :: i, compact_idx, part_idx, max_part_idx
   integer :: itemdir1_offset, itemdir2_offset, itemdir1_size, itemdir2_size
   logical, dimension(:), allocatable :: part_used

   ! Get itemdir info for this link
   itemdir1_offset = array_trees%links(link_idx)%itemdir1_offset
   itemdir2_offset = array_trees%links(link_idx)%itemdir2_offset
   itemdir1_size = array_trees%itemdir_size1
   itemdir2_size = array_trees%itemdir_size2

   ! Find maximum part index to size our tracking array
   max_part_idx = array_trees%total_parts
   allocate(part_used(max_part_idx))
   part_used = .false.

   ! Mark which parts are used in this link's itemdir
   do i = 1, itemdir1_size
      part_idx = array_trees%itemdir_entries(itemdir1_offset + i)
      if (part_idx > 0 .and. part_idx <= max_part_idx) then
         part_used(part_idx) = .true.
      end if
   end do

   do i = 1, itemdir2_size
      part_idx = array_trees%itemdir_entries(itemdir2_offset + i)
      if (part_idx > 0 .and. part_idx <= max_part_idx) then
         part_used(part_idx) = .true.
      end if
   end do

   ! Count how many parts are actually used
   compact_map%num_compact = count(part_used)

   ! Allocate compact arrays
   allocate(compact_map%compact_parts(compact_map%num_compact))
   allocate(compact_map%part_to_compact(max_part_idx))
   compact_map%part_to_compact = 0

   ! Build the compact mapping
   compact_idx = 0
   do part_idx = 1, max_part_idx
      if (part_used(part_idx)) then
         compact_idx = compact_idx + 1
         compact_map%compact_parts(compact_idx) = part_idx
         compact_map%part_to_compact(part_idx) = compact_idx
      end if
   end do

   ! Print the translation mapping
   write(stderr, '(A,I0,A,I0,A)') 'Link ', link_idx, ' compact map (', compact_map%num_compact, ' entries):'
   do compact_idx = 1, compact_map%num_compact
      part_idx = compact_map%compact_parts(compact_idx)
      write(stderr, '(A,I0,A,I0)') '  part ', part_idx, ' -> compact ', compact_idx
   end do

   deallocate(part_used)
end subroutine

subroutine build_compact_map_for_link_itemdir(array_trees, link_idx)
   ! This is called after itemdir_entries are populated to prepare for future reads
   ! For now, this is a placeholder - in a full implementation, we might store 
   ! these compact maps in array_trees for reuse
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: link_idx
   ! Placeholder - in full implementation, could store compact maps in array_trees
end subroutine

function get_compact_signature_from_itemdir1_array(array_trees, link_idx, atoms, item_value, compact_map) result(compact_signature)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: link_idx, item_value
   type(atom_type), dimension(:), intent(in) :: atoms
   type(compact_map_t), intent(in) :: compact_map
   integer, dimension(:), allocatable :: compact_signature
   integer :: i, adjlist_size, itemdir_offset, part_idx, compact_idx, sig_idx

   adjlist_size = size(atoms(item_value)%adjlist)
   allocate(compact_signature(adjlist_size))
   itemdir_offset = array_trees%links(link_idx)%itemdir1_offset

   sig_idx = 0
   do i = 1, adjlist_size
      part_idx = array_trees%itemdir_entries(itemdir_offset + atoms(item_value)%adjlist(i))
      if (part_idx > 0) then
         compact_idx = compact_map%part_to_compact(part_idx)
         if (compact_idx > 0) then
            sig_idx = sig_idx + 1
            compact_signature(sig_idx) = compact_idx
         end if
      end if
   end do

   ! Resize to actual compact size (removing unused slots)
   if (sig_idx < adjlist_size) then
      compact_signature = compact_signature(1:sig_idx)
   end if
end function

function get_compact_signature_from_itemdir2_array(array_trees, link_idx, atoms, item_value, compact_map) result(compact_signature)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: link_idx, item_value
   type(atom_type), dimension(:), intent(in) :: atoms
   type(compact_map_t), intent(in) :: compact_map
   integer, dimension(:), allocatable :: compact_signature
   integer :: i, adjlist_size, itemdir_offset, part_idx, compact_idx, sig_idx

   adjlist_size = size(atoms(item_value)%adjlist)
   allocate(compact_signature(adjlist_size))
   itemdir_offset = array_trees%links(link_idx)%itemdir2_offset

   sig_idx = 0
   do i = 1, adjlist_size
      part_idx = array_trees%itemdir_entries(itemdir_offset + atoms(item_value)%adjlist(i))
      if (part_idx > 0) then
         compact_idx = compact_map%part_to_compact(part_idx)
         if (compact_idx > 0) then
            sig_idx = sig_idx + 1
            compact_signature(sig_idx) = compact_idx
         end if
      end if
   end do

   ! Resize to actual compact size (removing unused slots)
   if (sig_idx < adjlist_size) then
      compact_signature = compact_signature(1:sig_idx)
   end if
end function

function find_child_part_compact_array(array_trees, parent_idx, compact_signature, compact_map) result(child_idx)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: parent_idx
   integer, dimension(:), intent(in) :: compact_signature
   type(compact_map_t), intent(in) :: compact_map
   integer :: child_idx

   child_idx = array_trees%parts(parent_idx)%first_child_idx
   
   do while (child_idx > 0)
      if (compact_signature_equivalence_array(array_trees, child_idx, compact_signature, compact_map)) then
         return
      end if
      child_idx = array_trees%parts(child_idx)%next_sibling_idx
   end do

   child_idx = 0
end function

function compact_signature_equivalence_array(array_trees, part_idx, target_compact_signature, compact_map) result(equiv)
   ! OPTIMIZED VERSION: Uses frequency-based hash comparison for O(n) performance
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   integer, dimension(:), intent(in) :: target_compact_signature
   type(compact_map_t), intent(in) :: compact_map
   logical :: equiv
   integer :: i, sig_len, target_len
   integer, dimension(:), allocatable :: part_compact_signature
   integer, dimension(:), allocatable :: part_freq, target_freq

   ! Convert part's signature to compact form
   call convert_part_signature_to_compact(array_trees, part_idx, compact_map, part_compact_signature)
   
   sig_len = size(part_compact_signature)
   target_len = size(target_compact_signature)
   
   if (sig_len /= target_len) then
      equiv = .false.
      return
   end if

   ! Quick exit for empty signatures
   if (sig_len == 0) then
      equiv = .true.
      return
   end if

   ! Build frequency arrays (size = max possible compact index)
   allocate(part_freq(compact_map%num_compact))
   allocate(target_freq(compact_map%num_compact))
   part_freq = 0
   target_freq = 0

   ! Count frequencies in both signatures
   do i = 1, sig_len
      part_freq(part_compact_signature(i)) = part_freq(part_compact_signature(i)) + 1
   end do

   do i = 1, target_len
      target_freq(target_compact_signature(i)) = target_freq(target_compact_signature(i)) + 1
   end do

   ! Compare frequency arrays - O(num_compact) comparison
   equiv = .true.
   do i = 1, compact_map%num_compact
      if (part_freq(i) /= target_freq(i)) then
         equiv = .false.
         exit
      end if
   end do

   deallocate(part_freq, target_freq)
end function

subroutine convert_part_signature_to_compact(array_trees, part_idx, compact_map, part_compact_signature)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   type(compact_map_t), intent(in) :: compact_map
   integer, dimension(:), allocatable, intent(out) :: part_compact_signature
   integer :: i, sig_len, part_sig_val, compact_idx, compact_count

   sig_len = array_trees%parts(part_idx)%signature_length
   
   ! First pass: count non-zero entries that map to compact indices
   compact_count = 0
   do i = 1, sig_len
      part_sig_val = array_trees%parts(part_idx)%signature(i)
      if (part_sig_val > 0 .and. part_sig_val <= size(compact_map%part_to_compact)) then
         compact_idx = compact_map%part_to_compact(part_sig_val)
         if (compact_idx > 0) then
            compact_count = compact_count + 1
         end if
      end if
   end do

   ! Allocate and fill compact signature
   allocate(part_compact_signature(compact_count))
   compact_count = 0
   do i = 1, sig_len
      part_sig_val = array_trees%parts(part_idx)%signature(i)
      if (part_sig_val > 0 .and. part_sig_val <= size(compact_map%part_to_compact)) then
         compact_idx = compact_map%part_to_compact(part_sig_val)
         if (compact_idx > 0) then
            compact_count = compact_count + 1
            part_compact_signature(compact_count) = compact_idx
         end if
      end if
   end do
end subroutine

! Legacy functions for backward compatibility - these are now less efficient
function get_signature_from_itemdir1_array(array_trees, link_idx, atoms, item_value) result(signature)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: link_idx, item_value
   type(atom_type), dimension(:), intent(in) :: atoms
   integer, dimension(:), allocatable :: signature
   integer :: i, adjlist_size, itemdir_offset

   adjlist_size = size(atoms(item_value)%adjlist)
   allocate(signature(adjlist_size))
   itemdir_offset = array_trees%links(link_idx)%itemdir1_offset

   do i = 1, adjlist_size
      signature(i) = array_trees%itemdir_entries(itemdir_offset + atoms(item_value)%adjlist(i))
   end do
end function

function get_signature_from_itemdir2_array(array_trees, link_idx, atoms, item_value) result(signature)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: link_idx, item_value
   type(atom_type), dimension(:), intent(in) :: atoms
   integer, dimension(:), allocatable :: signature
   integer :: i, adjlist_size, itemdir_offset

   adjlist_size = size(atoms(item_value)%adjlist)
   allocate(signature(adjlist_size))
   itemdir_offset = array_trees%links(link_idx)%itemdir2_offset

   do i = 1, adjlist_size
      signature(i) = array_trees%itemdir_entries(itemdir_offset + atoms(item_value)%adjlist(i))
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
   ! OPTIMIZED VERSION: Skip zero values (null pointers)
   type(array_trees_t), intent(in) :: array_trees
   integer, intent(in) :: part_idx
   integer, dimension(:), intent(in) :: target_signature
   logical :: equiv
   integer :: matches1, matches2, i, j, sig_len
   integer :: part_sig_val, target_val

   sig_len = array_trees%parts(part_idx)%signature_length
   
   if (sig_len /= size(target_signature)) then
      equiv = .false.
      return
   end if

   ! OPTIMIZED: Only compare non-zero signature entries
   do i = 1, sig_len
      part_sig_val = array_trees%parts(part_idx)%signature(i)

      ! Skip zero values (null pointers) - this is the key optimization
      if (part_sig_val == 0) cycle

      matches1 = 0
      matches2 = 0

      ! Count matches in target signature (skip zeros)
      do j = 1, sig_len
         target_val = target_signature(j)
         if (target_val /= 0 .and. part_sig_val == target_val) matches1 = matches1 + 1
      end do

      ! Count matches in part signature (skip zeros)
      do j = 1, sig_len
         if (array_trees%parts(part_idx)%signature(j) /= 0 .and. &
             part_sig_val == array_trees%parts(part_idx)%signature(j)) matches2 = matches2 + 1
      end do

      if (matches1 /= matches2) then
         equiv = .false.
         return
      end if
   end do

   equiv = .true.
end function

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
      call resplit_part_mna_array(mol1%atoms, mol2%atoms, array_trees, part_idx, link_idx, next_link_idx)
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

subroutine resplit_part_first_array(array_trees, part_idx, write_link_idx)
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

   ! Get child parts
   child_part1 = array_trees%parts(part_idx)%first_child_idx
   child_part2 = array_trees%parts(child_part1)%next_sibling_idx

   ! Get first items
   first_item1 = array_trees%item1_values(items1_offset + 1)
   first_item2 = array_trees%item2_values(items2_offset + 1)

   ! Add first items to first child
   array_trees%item1_values(array_trees%parts(child_part1)%items1_offset + 1) = first_item1
   array_trees%item2_values(array_trees%parts(child_part1)%items2_offset + 1) = first_item2
   array_trees%itemdir_entries(itemdir1_offset + first_item1) = child_part1
   array_trees%itemdir_entries(itemdir2_offset + first_item2) = child_part1

   ! Add remaining items to second child
   array_trees%parts(child_part2)%items1_fill_count = 0
   array_trees%parts(child_part2)%items2_fill_count = 0

   do i = 2, items1_count
      item_value = array_trees%item1_values(items1_offset + i)
      array_trees%parts(child_part2)%items1_fill_count = array_trees%parts(child_part2)%items1_fill_count + 1
      target_idx = array_trees%parts(child_part2)%items1_offset + array_trees%parts(child_part2)%items1_fill_count
      array_trees%item1_values(target_idx) = item_value
      array_trees%itemdir_entries(itemdir1_offset + item_value) = child_part2
   end do

   do i = 2, items2_count
      item_value = array_trees%item2_values(items2_offset + i)
      array_trees%parts(child_part2)%items2_fill_count = array_trees%parts(child_part2)%items2_fill_count + 1
      target_idx = array_trees%parts(child_part2)%items2_offset + array_trees%parts(child_part2)%items2_fill_count
      array_trees%item2_values(target_idx) = item_value
      array_trees%itemdir_entries(itemdir2_offset + item_value) = child_part2
   end do

   ! Build compact map for future reads from this link
   call build_compact_map_for_link_itemdir(array_trees, write_link_idx)
end subroutine

subroutine resplit_part_random_array(array_trees, part_idx, write_link_idx)
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

   ! Get child parts
   child_part1 = array_trees%parts(part_idx)%first_child_idx
   child_part2 = array_trees%parts(child_part1)%next_sibling_idx

   ! Get randomly chosen items
   chosen_item1 = array_trees%item1_values(items1_offset + random_index1)
   chosen_item2 = array_trees%item2_values(items2_offset + random_index2)

   write (stderr,*) part_idx, chosen_item1, chosen_item2

   ! Assign chosen items to first child
   array_trees%item1_values(array_trees%parts(child_part1)%items1_offset + 1) = chosen_item1
   array_trees%item2_values(array_trees%parts(child_part1)%items2_offset + 1) = chosen_item2
   array_trees%itemdir_entries(itemdir1_offset + chosen_item1) = child_part1
   array_trees%itemdir_entries(itemdir2_offset + chosen_item2) = child_part1

   ! Assign remaining items to second child
   array_trees%parts(child_part2)%items1_fill_count = 0
   do i = 1, num_items1
      if (i /= random_index1) then
         item_value = array_trees%item1_values(items1_offset + i)
         array_trees%parts(child_part2)%items1_fill_count = array_trees%parts(child_part2)%items1_fill_count + 1
         target_idx = array_trees%parts(child_part2)%items1_offset + array_trees%parts(child_part2)%items1_fill_count
         array_trees%item1_values(target_idx) = item_value
         array_trees%itemdir_entries(itemdir1_offset + item_value) = child_part2
      end if
   end do

   array_trees%parts(child_part2)%items2_fill_count = 0
   do i = 1, num_items2
      if (i /= random_index2) then
         item_value = array_trees%item2_values(items2_offset + i)
         array_trees%parts(child_part2)%items2_fill_count = array_trees%parts(child_part2)%items2_fill_count + 1
         target_idx = array_trees%parts(child_part2)%items2_offset + array_trees%parts(child_part2)%items2_fill_count
         array_trees%item2_values(target_idx) = item_value
         array_trees%itemdir_entries(itemdir2_offset + item_value) = child_part2
      end if
   end do

   ! Build compact map for future reads from this link
   call build_compact_map_for_link_itemdir(array_trees, write_link_idx)
end subroutine

recursive subroutine redistribute_items_array(mol1, mol2, array_trees, branch_idx)
   type(mol_type), intent(in) :: mol1, mol2
   type(array_trees_t), intent(inout) :: array_trees
   integer, intent(in) :: branch_idx
   integer :: child_branch_idx, first_link_idx, split_part_idx

   child_branch_idx = array_trees%chains(branch_idx)%first_child_idx
   
   do while (child_branch_idx > 0)
      first_link_idx = array_trees%chains(child_branch_idx)%link_offset + 1
      split_part_idx = array_trees%chains(child_branch_idx)%split_part_idx
      
      call resplit_part_random_array(array_trees, split_part_idx, first_link_idx)
      call recompute_consistent_mnas_array(mol1, mol2, array_trees, child_branch_idx)
      call redistribute_items_array(mol1, mol2, array_trees, child_branch_idx)
      
      child_branch_idx = array_trees%chains(child_branch_idx)%next_sibling_idx
   end do
end subroutine

end module
