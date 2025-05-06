module eltype_compute
use parameters
use chemdata
use molecule
use lcrs_tree
implicit none

type :: atomtype_item
   integer :: elnum
   integer :: label
   type(part_node_t), pointer :: part
end type

type :: atomtype_table
   integer :: num_items
   type(atomtype_item), dimension(:), allocatable :: items
end type

contains

subroutine add_atomtype(atomtypetable, elnum, label, part)
   type(atomtype_table), intent(inout) :: atomtypetable
   integer, intent(in) :: elnum
   integer, intent(in) :: label
   type(part_node_t), pointer, intent(in) :: part

   atomtypetable%num_items = atomtypetable%num_items + 1
   atomtypetable%items(atomtypetable%num_items)%elnum = elnum
   atomtypetable%items(atomtypetable%num_items)%label = label
   atomtypetable%items(atomtypetable%num_items)%part => part
end subroutine

function find_atomtype(atomtypetable, elnum, label) result(part)
   type(atomtype_table), intent(in) :: atomtypetable
   integer, intent(in) :: elnum
   integer, intent(in) :: label
   type(part_node_t), pointer :: part
   integer :: i

   do i = 1, atomtypetable%num_items
      if (atomtypetable%items(i)%elnum == elnum .and. &
          atomtypetable%items(i)%label == label) then
         part => atomtypetable%items(i)%part
         return
      end if
   end do

   part => null()
end function

subroutine compute_eltypes(mol1, mol2, eltypes)
! Partition atoms by atomic number and label
   type(mol_type), intent(in) :: mol1, mol2
   type(partitionarray_t), intent(out) :: eltypes
   ! Local variables
   type(tree_node_t), pointer :: tree_root
   type(link_node_t), pointer :: new_link
   type(part_node_t), pointer :: part
   type(atomtype_table) :: atomtypetable
   integer :: i, elnum, label, num_atoms1, num_atoms2

   num_atoms1 = size(mol1%atoms)
   num_atoms2 = size(mol2%atoms)

   tree_root => make_new_tree(num_atoms1, num_atoms2)
   new_link => add_new_link(tree_root)
   allocate(atomtypetable%items(num_atoms1 + num_atoms2))
   atomtypetable%num_items = 0

   ! First molecule
   do i = 1, num_atoms1
      elnum = mol1%atoms(i)%elnum
      label = mol1%atoms(i)%label
      part => find_atomtype(atomtypetable, elnum, label)
      if (.not. associated(part)) then
         part => add_new_part(new_link)
         call add_atomtype(atomtypetable, elnum, label, part)
      end if
      call add_new_item1(part, i)
   end do

   ! Second molecule
   do i = 1, num_atoms2
      elnum = mol2%atoms(i)%elnum
      label = mol2%atoms(i)%label
      part => find_atomtype(atomtypetable, elnum, label)
      if (.not. associated(part)) then
         part => add_new_part(new_link)
         call add_atomtype(atomtypetable, elnum, label, part)
      end if
      call add_new_item2(part, i)
   end do

   call partition_to_partitionarray(new_link, eltypes)

   ! Clean up
   call delete_tree(tree_root)
end subroutine

end module
