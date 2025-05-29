module eltype_compute
use parameters
use chemdata
use molecule
use lcrs_tree
implicit none

type :: atomtype_item
   integer :: elnum
   integer :: label
   integer :: part_index
end type

type :: atomtype_table
   integer :: num_items
   type(atomtype_item), dimension(:), allocatable :: items
end type

contains

subroutine add_atomtype(atomtypetable, atom, part_index)
   type(atomtype_table), intent(inout) :: atomtypetable
   type(atom_type), intent(in) :: atom
   integer, intent(in) :: part_index

   atomtypetable%num_items = atomtypetable%num_items + 1
   atomtypetable%items(atomtypetable%num_items)%elnum = atom%elnum
   atomtypetable%items(atomtypetable%num_items)%label = atom%label
   atomtypetable%items(atomtypetable%num_items)%part_index = part_index
end subroutine

function find_atomtype(atomtypetable, atom) result(part_index)
   type(atomtype_table), intent(in) :: atomtypetable
   type(atom_type), intent(in) :: atom
   integer :: part_index
   integer :: i

   do i = 1, atomtypetable%num_items
      if (atomtypetable%items(i)%elnum == atom%elnum .and. &
          atomtypetable%items(i)%label == atom%label) then
         part_index = atomtypetable%items(i)%part_index
         return
      end if
   end do

   part_index = 0  ! Not found
end function

subroutine set_eltypes(mol1, mol2, eltypes)
! Partition atoms by atomic number and label using arrays directly
   type(mol_type), intent(in) :: mol1, mol2
   type(partitionarray_t), intent(out) :: eltypes
   ! Local variables
   type(atomtype_table) :: atomtypetable
   integer :: i, num_atoms1, num_atoms2
   integer :: part_index, current_part
   integer :: max_parts
   ! Temporary arrays for building partitions
   integer, dimension(:), allocatable :: part_num_items1, part_num_items2
   integer, dimension(:,:), allocatable :: part_items1, part_items2
   integer, dimension(:), allocatable :: itemdir1_temp, itemdir2_temp

   num_atoms1 = size(mol1%atoms)
   num_atoms2 = size(mol2%atoms)
   max_parts = num_atoms1 + num_atoms2  ! Maximum possible partitions

   ! Allocate temporary arrays
   allocate(atomtypetable%items(max_parts))
   allocate(part_num_items1(max_parts))
   allocate(part_num_items2(max_parts))
   allocate(part_items1(num_atoms1, max_parts))
   allocate(part_items2(num_atoms2, max_parts))
   allocate(itemdir1_temp(num_atoms1))
   allocate(itemdir2_temp(num_atoms2))

   ! Initialize
   atomtypetable%num_items = 0
   part_num_items1 = 0
   part_num_items2 = 0
   current_part = 0

   ! First molecule
   do i = 1, num_atoms1
      part_index = find_atomtype(atomtypetable, mol1%atoms(i))

      if (part_index == 0) then
         ! Create new partition
         current_part = current_part + 1
         call add_atomtype(atomtypetable, mol1%atoms(i), current_part)
         part_index = current_part
      end if

      ! Add item to partition
      part_num_items1(part_index) = part_num_items1(part_index) + 1
      part_items1(part_num_items1(part_index), part_index) = i
      itemdir1_temp(i) = part_index
   end do

   ! Second molecule
   do i = 1, num_atoms2
      part_index = find_atomtype(atomtypetable, mol2%atoms(i))

      if (part_index == 0) then
         ! Create new partition
         current_part = current_part + 1
         call add_atomtype(atomtypetable, mol2%atoms(i), current_part)
         part_index = current_part
      end if

      ! Add item to partition
      part_num_items2(part_index) = part_num_items2(part_index) + 1
      part_items2(part_num_items2(part_index), part_index) = i
      itemdir2_temp(i) = part_index
   end do

   ! Now build the final partitionarray_t structure
   eltypes%num_parts = current_part
   allocate(eltypes%parts(eltypes%num_parts))
   allocate(eltypes%itemdir1(num_atoms1))
   allocate(eltypes%itemdir2(num_atoms2))

   ! Copy item directories
   eltypes%itemdir1 = itemdir1_temp
   eltypes%itemdir2 = itemdir2_temp

   ! Build each partition
   do i = 1, eltypes%num_parts
      ! Set sizes
      eltypes%parts(i)%num_items1 = part_num_items1(i)
      eltypes%parts(i)%num_items2 = part_num_items2(i)
      eltypes%parts(i)%num_children = 0

      ! Allocate and copy items
      if (part_num_items1(i) > 0) then
         allocate(eltypes%parts(i)%items1(part_num_items1(i)))
         eltypes%parts(i)%items1 = part_items1(1:part_num_items1(i), i)
      else
         allocate(eltypes%parts(i)%items1(0))
      end if

      if (part_num_items2(i) > 0) then
         allocate(eltypes%parts(i)%items2(part_num_items2(i)))
         eltypes%parts(i)%items2 = part_items2(1:part_num_items2(i), i)
      else
         allocate(eltypes%parts(i)%items2(0))
      end if

      ! Allocate empty arrays for neighbors and children
      allocate(eltypes%parts(i)%signature(0))
      allocate(eltypes%parts(i)%children(0))
   end do

   ! Clean up temporary arrays
   deallocate(part_num_items1)
   deallocate(part_num_items2)
   deallocate(part_items1)
   deallocate(part_items2)
   deallocate(itemdir1_temp)
   deallocate(itemdir2_temp)
end subroutine

end module
