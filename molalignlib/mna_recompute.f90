module mna_recompute
use parameters
use molecule
use lcrs_tree
implicit none

interface operator (.equiv.)
   module procedure neighbor_array_equivalence
end interface

contains

subroutine chain_to_chainarray(chain_root, chainarray)
   type(chain_root_t), target, intent(in) :: chain_root
   type(chainarray_t), target, intent(out) :: chainarray
   type(link_node_t), pointer :: curr_link
   type(part_node_t), pointer :: part
   type(item_node_t), pointer :: item
   type(part_node_t), pointer :: child_part
   integer :: link_idx, part_idx, item_idx, neighbor_idx, child_idx

   ! Initialize chainarray structure
   allocate(chainarray%partitions(chain_root%num_links))
   chainarray%num_links = chain_root%num_links
   chainarray%tot_items1 = chain_root%tot_items1
   chainarray%tot_items2 = chain_root%tot_items2

   ! Single direct loop over all links
   curr_link => chain_root%first_link
   do link_idx = 1, chain_root%num_links

      ! Initialize current partition
      chainarray%partitions(link_idx)%num_parts = curr_link%num_parts
      allocate(chainarray%partitions(link_idx)%parts(curr_link%num_parts))
      allocate(chainarray%partitions(link_idx)%itemdir1(chain_root%tot_items1))
      allocate(chainarray%partitions(link_idx)%itemdir2(chain_root%tot_items2))

      ! Process each part in current link using part index
      part => curr_link%first_part
      do while (associated(part))
         part_idx = part%index  ! Leverage the index value from part_node_t

         ! Initialize part structure
         chainarray%partitions(link_idx)%parts(part_idx)%num_items1 = part%num_items1
         chainarray%partitions(link_idx)%parts(part_idx)%num_items2 = part%num_items2
         chainarray%partitions(link_idx)%parts(part_idx)%num_neighbors = part%num_neighbors
         chainarray%partitions(link_idx)%parts(part_idx)%num_children = part%num_children

         ! Allocate arrays
         allocate(chainarray%partitions(link_idx)%parts(part_idx)%items1(part%num_items1))
         allocate(chainarray%partitions(link_idx)%parts(part_idx)%items2(part%num_items2))
         if (part%num_neighbors > 0) then
            allocate(chainarray%partitions(link_idx)%parts(part_idx)%neighbors(part%num_neighbors))
         end if
         if (part%num_children > 0) then
            allocate(chainarray%partitions(link_idx)%parts(part_idx)%children(part%num_children))
         end if

         ! Copy items1 and update itemdir1
         item => part%first_item1
         do item_idx = 1, part%num_items1
            chainarray%partitions(link_idx)%parts(part_idx)%items1(item_idx) = item%value
            chainarray%partitions(link_idx)%itemdir1(item%value) = part_idx
            item => item%next_item
         end do

         ! Copy items2 and update itemdir2
         item => part%first_item2
         do item_idx = 1, part%num_items2
            chainarray%partitions(link_idx)%parts(part_idx)%items2(item_idx) = item%value
            chainarray%partitions(link_idx)%itemdir2(item%value) = part_idx
            item => item%next_item
         end do

         ! Copy neighborhood information
         do neighbor_idx = 1, part%num_neighbors
            chainarray%partitions(link_idx)%parts(part_idx)%neighbors(neighbor_idx) = &
               part%neighbors(neighbor_idx)%ptr%index
         end do

         ! Copy children information (part indices in next link)
         child_part => part%first_child
         do child_idx = 1, part%num_children
            chainarray%partitions(link_idx)%parts(part_idx)%children(child_idx) = child_part%index
            child_part => child_part%next_sibling
         end do

         part => part%next_part
      end do

      curr_link => curr_link%next_link
   end do
end subroutine

! Check if two neighbors arrays are equivalent (similar to .equiv. operator)
! This matches the original implementation by checking for set equivalence
! without considering order
function neighbor_array_equivalence(array1, array2) result(equiv)
   integer, dimension(:), intent(in) :: array1, array2
   logical :: equiv
   integer :: i, j, matches

   if (size(array1) /= size(array2)) then
      equiv = .false.
      return
   end if

   do i = 1, size(array1)
      matches = 0
      do j = 1, size(array1)
         if (array1(i) == array2(j)) matches = matches + 1
         if (array1(i) == array1(j)) matches = matches - 1
      end do
      if (matches /= 0) then
         equiv = .false.
         return
      end if
   end do

   equiv = .true.
end function

function find_child_part_array(chainarray, part_idx, neighbors) result(child_part_idx)
   type(chainarray_t), intent(in) :: chainarray
   integer, intent(in) :: part_idx
   integer, dimension(:), intent(in) :: neighbors
   integer :: child_part_idx
   ! Local variables
   integer :: i

   ! Check each part child of the current part
   do i = 1, chainarray%partitions(chainarray%num_links)%parts(part_idx)%num_children
      child_part_idx = chainarray%partitions(chainarray%num_links)%parts(part_idx)%children(i)
      ! Check for neighbors equivalence
      if (neighbors .equiv. chainarray%partitions(chainarray%num_links+1)%parts(child_part_idx)%neighbors) return
   end do

   ! No matching part child found - this is an error
   error stop 'find_child_part_array: No matching part child found'
end function

subroutine reset_item_counts(chainarray, link_idx)
   type(chainarray_t), intent(inout) :: chainarray
   integer, intent(in) :: link_idx
   integer :: i

   ! Reset item counts for all parts in this link_idx
   do i = 1, chainarray%partitions(link_idx)%num_parts
      chainarray%partitions(link_idx)%parts(i)%num_items1 = 0
      chainarray%partitions(link_idx)%parts(i)%num_items2 = 0
   end do
end subroutine

subroutine recompute_nextlevel_mnas(mol1, mol2, chainarray)
! Compute next level MNA types - replicates compute_nextlevel_mnas logic
   type(mol_type), intent(in) :: mol1, mol2
   type(chainarray_t), target, intent(inout) :: chainarray
   ! Local variables
   integer :: part_idx, item_idx, child_part_idx, i, j, curr_link_idx, next_link_idx
   integer, dimension(:), allocatable :: neighbors
   logical :: child_found
   type(partitionarray_t), pointer :: curr_part, next_part

   curr_link_idx = chainarray%num_links
   next_link_idx = chainarray%num_links + 1
   chainarray%num_links = chainarray%num_links + 1
   curr_part => chainarray%partitions(curr_link_idx)
   next_part => chainarray%partitions(next_link_idx)

   ! Reset item counts in the new link
   call reset_item_counts(chainarray, next_link_idx)

   ! Process each parent part in the current link
   do part_idx = 1, curr_part%num_parts

      ! First molecule - process each item in this parent part
      do i = 1, curr_part%parts(part_idx)%num_items1
         item_idx = curr_part%parts(part_idx)%items1(i)
         neighbors = curr_part%itemdir1(mol1%atoms(item_idx)%adjlist)

         ! Find or create child part with these neighbors
         child_found = .false.
         do j = 1, curr_part%parts(part_idx)%num_children
            child_part_idx = curr_part%parts(part_idx)%children(j)

            if (next_part%parts(child_part_idx)%num_neighbors > 0) then
               if (neighbors .equiv. next_part%parts(child_part_idx)%neighbors) then
                  child_found = .true.
                  exit
               end if
            else
               ! This child doesn't have neighbors set yet, use it
               next_part%parts(child_part_idx)%num_neighbors = size(neighbors)
               if (size(neighbors) > 0) then
                  allocate(next_part%parts(child_part_idx)%neighbors(size(neighbors)))
                  next_part%parts(child_part_idx)%neighbors = neighbors
               end if
               child_found = .true.
               exit
            end if
         end do

         if (.not. child_found) error stop 'No matching child part found'

         ! Add item to the child part
         next_part%parts(child_part_idx)%num_items1 = next_part%parts(child_part_idx)%num_items1 + 1
         next_part%parts(child_part_idx)%items1(next_part%parts(child_part_idx)%num_items1) = item_idx
         next_part%itemdir1(item_idx) = child_part_idx
      end do

      ! Second molecule - process each item in this parent part
      do i = 1, curr_part%parts(part_idx)%num_items2
         item_idx = curr_part%parts(part_idx)%items2(i)
         neighbors = curr_part%itemdir2(mol2%atoms(item_idx)%adjlist)

         ! Find child part with these neighbors
         child_found = .false.
         do j = 1, curr_part%parts(part_idx)%num_children
            child_part_idx = curr_part%parts(part_idx)%children(j)

            if (next_part%parts(child_part_idx)%num_neighbors > 0) then
               if (neighbors .equiv. next_part%parts(child_part_idx)%neighbors) then
                  child_found = .true.
                  exit
               end if
            else
               ! This child doesn't have neighbors set yet, use it
               next_part%parts(child_part_idx)%num_neighbors = size(neighbors)
               if (size(neighbors) > 0) then
                  allocate(next_part%parts(child_part_idx)%neighbors(size(neighbors)))
                  next_part%parts(child_part_idx)%neighbors = neighbors
               end if
               child_found = .true.
               exit
            end if
         end do

         if (.not. child_found) error stop 'No matching child part found'

         ! Add item to the child part
         next_part%parts(child_part_idx)%num_items2 = next_part%parts(child_part_idx)%num_items2 + 1
         next_part%parts(child_part_idx)%items2(next_part%parts(child_part_idx)%num_items2) = item_idx
         next_part%itemdir2(item_idx) = child_part_idx
      end do
   end do
end subroutine

subroutine recompute_consistent_mnas(mol1, mol2, chainarray)
! Iteratively compute MNA types - replicates compute_consistent_mnas logic
   type(mol_type), intent(in) :: mol1, mol2
   type(chainarray_t), intent(inout) :: chainarray
   integer :: prev_num_parts

   do
      ! Save the number of parts before computation
      prev_num_parts = chainarray%partitions(chainarray%num_links)%num_parts

      ! Compute next level and update chainarray
      call recompute_nextlevel_mnas(mol1, mol2, chainarray)

      ! Exit the loop if no change (convergence reached)
      if (chainarray%partitions(chainarray%num_links)%num_parts == prev_num_parts) exit
   end do
end subroutine

end module
