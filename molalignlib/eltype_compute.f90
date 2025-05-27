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

subroutine add_atomtype(atomtypetable, atom, part)
   type(atomtype_table), intent(inout) :: atomtypetable
   type(atom_type), intent(in) :: atom
   type(part_node_t), pointer, intent(in) :: part

   atomtypetable%num_items = atomtypetable%num_items + 1
   atomtypetable%items(atomtypetable%num_items)%elnum = atom%elnum
   atomtypetable%items(atomtypetable%num_items)%label = atom%label
   atomtypetable%items(atomtypetable%num_items)%part => part
end subroutine

function find_atomtype(atomtypetable, atom) result(part)
   type(atomtype_table), intent(in) :: atomtypetable
   type(atom_type), intent(in) :: atom
   type(part_node_t), pointer :: part
   integer :: i

   do i = 1, atomtypetable%num_items
      if (atomtypetable%items(i)%elnum == atom%elnum .and. &
          atomtypetable%items(i)%label == atom%label) then
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
   type(chain_root_t), pointer :: chain_root

   chain_root => eltypetree(mol1, mol2)
   call partition_to_partitionarray(chain_root%last_link, eltypes)
   call delete_chain(chain_root)
end subroutine

function eltypetree(mol1, mol2) result(chain_root)
! Partition atoms by atomic number and label
   type(mol_type), intent(in) :: mol1, mol2
   ! Local variables
   type(part_node_t), pointer :: root_part, child_part
   type(chain_root_t), pointer :: chain_root
   type(link_node_t), pointer :: new_link
   type(item_node_t), pointer :: item1, item2
   type(atomtype_table) :: atomtypetable
   integer :: i, num_atoms1, num_atoms2

   num_atoms1 = size(mol1%atoms)
   num_atoms2 = size(mol2%atoms)

   root_part => make_new_part()
   chain_root => make_chain_root(num_atoms1, num_atoms2)
   new_link => add_new_link(chain_root)
   call add_part(new_link, root_part)
   new_link => add_new_link(chain_root)

   allocate(atomtypetable%items(num_atoms1 + num_atoms2))
   atomtypetable%num_items = 0

   do i = 1, num_atoms1
      call add_new_item1(root_part, i)
   end do

   do i = 1, num_atoms2
      call add_new_item2(root_part, i)
   end do

   ! First molecule
   item1 => root_part%first_item1
   do while(associated(item1))
      child_part => find_atomtype(atomtypetable, mol1%atoms(item1%value))
      if (.not. associated(child_part)) then
         child_part => add_new_part(new_link, root_part)
         call add_atomtype(atomtypetable, mol1%atoms(item1%value), child_part)
      end if
      call add_new_item1(child_part, item1%value)
      item1 => item1%next_item
   end do

   ! Second molecule
   item2 => root_part%first_item2
   do while(associated(item2))
      child_part => find_atomtype(atomtypetable, mol2%atoms(item2%value))
      if (.not. associated(child_part)) then
         child_part => add_new_part(new_link, root_part)
         call add_atomtype(atomtypetable, mol2%atoms(item2%value), child_part)
      end if
      call add_new_item2(child_part, item2%value)
      item2 => item2%next_item
   end do
end function

end module
