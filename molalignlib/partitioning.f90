module partitioning
use parameters
use derived_types
use chemistry
use adjacency
use molecule
use lcrs_trees
use options
implicit none
private
public collect_atomtypes
public refine_mlna_part
public refine_mlna_partition
public compute_scna_partition

type :: atomtype_item_t
   integer :: elnum
   integer :: typeid
   type(partree_node_t), pointer :: part
end type

type :: atomtype_table_t
   integer :: num_items
   type(atomtype_item_t), dimension(:), allocatable :: items
end type

contains

subroutine add_atomtype(atomtypetable, atom, part)
   type(atomtype_table_t), intent(inout) :: atomtypetable
   type(atom_t), intent(in) :: atom
   type(partree_node_t), pointer, intent(in) :: part

   atomtypetable%num_items = atomtypetable%num_items + 1
   atomtypetable%items(atomtypetable%num_items)%elnum = atom%elnum
   atomtypetable%items(atomtypetable%num_items)%typeid = atom%typeid
   atomtypetable%items(atomtypetable%num_items)%part => part
end subroutine

function find_atomtype(atomtypetable, atom) result(part)
   type(atomtype_table_t), intent(in) :: atomtypetable
   type(atom_t), intent(in) :: atom
   type(partree_node_t), pointer :: part
   integer :: i

   do i = 1, atomtypetable%num_items
      if (atomtypetable%items(i)%elnum == atom%elnum) then
         if (.not. label_flag .or. atomtypetable%items(i)%typeid == atom%typeid) then
            part => atomtypetable%items(i)%part
            return
         end if
      end if
   end do

   part => null()
end function

subroutine build_atomtypes_tree(atomset1, atomset2, atoms1, atoms2, chain_root, root_part)
! Partition atoms by atomic number and label using linked list structures
   integer, dimension(:), intent(in) :: atomset1, atomset2
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(assigntree_node_t), pointer, intent(out) :: chain_root
   type(partree_node_t), pointer, intent(out) :: root_part
   ! Local variables
   type(partree_node_t), pointer :: child_part
   type(chain_node_t), pointer :: new_link
   type(atomtype_table_t) :: atomtypetable
   integer :: num_atoms1, num_atoms2, i, atomidx

   num_atoms1 = size(atoms1)
   num_atoms2 = size(atoms2)

   root_part => new_root_part()
   chain_root => new_root_chain(num_atoms1, num_atoms2)
   new_link => new_chain_link(chain_root)

   allocate(atomtypetable%items(num_atoms1 + num_atoms2))
   atomtypetable%num_items = 0

   ! First molecule
   do i = 1, size(atomset1)
      atomidx = atomset1(i)
      child_part => find_atomtype(atomtypetable, atoms1(atomidx))
      if (.not. associated(child_part)) then
         child_part => new_child_part(root_part)
         call link_part(new_link, child_part)
         call add_atomtype(atomtypetable, atoms1(atomidx), child_part)
      end if
      call add_new_item1(child_part, atomidx)
      new_link%itemdir1(atomidx)%ptr => child_part
   end do

   ! Second molecule
   do i = 1, size(atomset2)
      atomidx = atomset2(i)
      child_part => find_atomtype(atomtypetable, atoms2(atomidx))
      if (.not. associated(child_part)) then
         child_part => new_child_part(root_part)
         call link_part(new_link, child_part)
         call add_atomtype(atomtypetable, atoms2(atomidx), child_part)
      end if
      call add_new_item2(child_part, atomidx)
      new_link%itemdir2(atomidx)%ptr => child_part
   end do

   deallocate(atomtypetable%items)
end subroutine

subroutine collect_atomtypes(atomset1, atomset2, atoms1, atoms2, atomtypes)
! Partition atoms by atomic number and label
! Uses linked list structures internally, then converts to partition array
   integer, dimension(:), intent(in) :: atomset1, atomset2
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_t), intent(out) :: atomtypes
   ! Local variables
   type(assigntree_node_t), pointer :: chain_root
   type(partree_node_t), pointer :: root_part
   type(chain_node_t), pointer :: first_link

   ! Create linked list structures
   call build_atomtypes_tree(atomset1, atomset2, atoms1, atoms2, chain_root, root_part)

   ! Get the first (and only) link from the chain
   first_link => chain_root%first_link

   ! Convert link to partition array structure
   call link_to_partition(first_link, atomtypes)

   ! Clean up tree structures
   call delete_chain(chain_root)
   call delete_part_tree(root_part)
end subroutine

subroutine refine_mlna_part(adjcs1, adjcs2, itemdir1, itemdir2, part, link)
! Create children for different signatures - caller decides what to do with them
! Note: part is always a leaf part with no existing children
   type(adjcs_t), intent(in) :: adjcs1, adjcs2
   type(part_nodeptr_t), dimension(:), intent(in) :: itemdir1, itemdir2
   type(partree_node_t), pointer, intent(inout) :: part
   type(chain_node_t), pointer, intent(inout) :: link
   ! Local variables
   type(item_node_t), pointer :: item
   type(partree_node_t), pointer :: child_part
   type(part_nodeptr_t), dimension(:), allocatable :: signature

   ! Process first molecule items - create children for each unique signature
   item => part%first_item1
   do while (associated(item))
      signature = itemdir1(adjcs1%lists(:adjcs1%cns(item%idx), item%idx))
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
      signature = itemdir2(adjcs2%lists(:adjcs2%cns(item%idx), item%idx))
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

subroutine refine_mlna_partition(adjcs1, adjcs2, mlnachain, num_splits)
! Compute next level MLNA types - always keeps all children (original behavior)
   type(adjcs_t), intent(in) :: adjcs1, adjcs2
   type(assigntree_node_t), pointer, intent(inout) :: mlnachain
   integer, intent(out) :: num_splits
   ! Local variables
   type(chain_node_t), pointer :: link, new_link
   type(partref_node_t), pointer :: partref

   num_splits = 0

   ! Save the last link before creating a new one
   link => mlnachain%last_link
   new_link => new_chain_link(mlnachain)

   ! Process all parts in the current partition
   partref => link%first_partref
   do while (associated(partref))
      ! Create children based on signatures
      call refine_mlna_part(adjcs1, adjcs2, link%itemdir1, link%itemdir2, partref%part, new_link)

      ! Count splits (children beyond the original part)
      num_splits = num_splits + partref%part%num_children - 1

      partref => partref%nextref
   end do
end subroutine

subroutine compute_scna_partition(adjcs1, adjcs2, atomtypes, mlnachain)
! Iteratively compute MLNA types until convergence
   type(adjcs_t), intent(in) :: adjcs1, adjcs2
   type(partition_t), intent(in) :: atomtypes
   type(assigntree_node_t), pointer, intent(out) :: mlnachain
   ! Local variables
   integer :: num_splits

!   mlnachain => collect_atomtypes_linked( adjcs1, adjcs2)
   mlnachain => chain_from_partition( atomtypes)

   do
      ! Call refine_mlna_partition and get the number of splits
      call refine_mlna_partition(adjcs1, adjcs2, mlnachain, num_splits)

      ! Exit loop if no splits occurred
      if (num_splits == 0) exit
   end do
end subroutine

end module
