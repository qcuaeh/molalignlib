module partitioning_ll
use parameters
use derived_types
use molecule
use lcrs_trees
implicit none
private
public collect_atomtypes_ll

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
      if (atomtypetable%items(i)%elnum == atom%elnum .and. &
          atomtypetable%items(i)%typeid == atom%typeid) then
         part => atomtypetable%items(i)%part
         return
      end if
   end do

   part => null()
end function

function collect_atomtypes_ll(atoms1, atoms2) result(chain_root)
! Partition atoms by atomic number and label
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   ! Local variables
   type(partree_node_t), pointer :: root_part, child_part
   type(assigntree_node_t), pointer :: chain_root
   type(chain_node_t), pointer :: new_link
   type(atomtype_table_t) :: atomtypetable
   integer :: num_atoms1, num_atoms2, i

   num_atoms1 = size(atoms1)
   num_atoms2 = size(atoms2)

   root_part => new_root_part()
   chain_root => new_root_chain(num_atoms1, num_atoms2)
   new_link => new_chain_link(chain_root)

   allocate(atomtypetable%items(num_atoms1 + num_atoms2))
   atomtypetable%num_items = 0

   ! First molecule
   do i = 1, num_atoms1
      child_part => find_atomtype(atomtypetable, atoms1(i))
      if (.not. associated(child_part)) then
         child_part => new_child_part(root_part)
         call link_part(new_link, child_part)
         call add_atomtype(atomtypetable, atoms1(i), child_part)
      end if
      call add_new_item1(child_part, i)
   end do

   ! Second molecule
   do i = 1, num_atoms1
      child_part => find_atomtype(atomtypetable, atoms2(i))
      if (.not. associated(child_part)) then
         child_part => new_child_part(root_part)
         call link_part(new_link, child_part)
         call add_atomtype(atomtypetable, atoms2(i), child_part)
      end if
      call add_new_item2(child_part, i)
   end do
end function

end module
