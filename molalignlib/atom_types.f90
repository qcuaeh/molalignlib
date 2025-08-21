module atom_types
use options
use parameters
use derived_types
use chemistry
use molecule
use lcrs_tree
implicit none
private
public collect_atomtypes

type :: atomtype_item_t
   integer :: elnum
   integer :: typeid
   integer :: partidx
end type

type :: atomtype_table_t
   integer :: num_items
   type(atomtype_item_t), dimension(:), allocatable :: items
end type

contains

subroutine add_atomtype(atomtypetable, atom, partidx)
   type(atomtype_table_t), intent(inout) :: atomtypetable
   type(atom_t), intent(in) :: atom
   integer, intent(in) :: partidx

   atomtypetable%num_items = atomtypetable%num_items + 1
   atomtypetable%items(atomtypetable%num_items)%elnum = atom%elnum
   atomtypetable%items(atomtypetable%num_items)%typeid = atom%typeid
   atomtypetable%items(atomtypetable%num_items)%partidx = partidx
end subroutine

function find_atomtype(atomtypetable, atom) result(partidx)
   type(atomtype_table_t), intent(in) :: atomtypetable
   type(atom_t), intent(in) :: atom
   integer :: partidx
   integer :: i

   do i = 1, atomtypetable%num_items
      if (atomtypetable%items(i)%elnum == atom%elnum .and. &
          atomtypetable%items(i)%typeid == atom%typeid) then
         partidx = atomtypetable%items(i)%partidx
         return
      end if
   end do

   partidx = 0  ! Not found
end function

subroutine collect_atomtypes(atoms1, atoms2, atomtypes)
! Partition atoms by atomic number and label using arrays directly
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_t), intent(out) :: atomtypes
   ! Local variables
   type(atomtype_table_t) :: atomtypetable
   integer :: i, num_atoms1, num_atoms2
   integer :: partidx, current_part
   integer :: max_parts
   ! Temporary arrays for building partitions
   integer, dimension(:), allocatable :: part_num_items1, part_num_items2
   integer, dimension(:,:), allocatable :: part_items1, part_items2
   integer, dimension(:), allocatable :: itemdir1_temp, itemdir2_temp

   num_atoms1 = size(atoms1)
   num_atoms2 = size(atoms2)
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
      if (atoms1(i)%mask) then
         partidx = find_atomtype(atomtypetable, atoms1(i))
         if (partidx == 0) then
            ! Create new partition
            current_part = current_part + 1
            call add_atomtype(atomtypetable, atoms1(i), current_part)
            partidx = current_part
         end if
         ! Add item to partition
         part_num_items1(partidx) = part_num_items1(partidx) + 1
         part_items1(part_num_items1(partidx), partidx) = i
         itemdir1_temp(i) = partidx
      end if
   end do

   ! Second molecule
   do i = 1, num_atoms2
      if (atoms2(i)%mask) then
         partidx = find_atomtype(atomtypetable, atoms2(i))
         if (partidx == 0) then
            ! Create new partition
            current_part = current_part + 1
            call add_atomtype(atomtypetable, atoms2(i), current_part)
            partidx = current_part
         end if
         ! Add item to partition
         part_num_items2(partidx) = part_num_items2(partidx) + 1
         part_items2(part_num_items2(partidx), partidx) = i
         itemdir2_temp(i) = partidx
      end if
   end do

   ! Now build the final partition_t structure
   atomtypes%num_parts = current_part
   allocate(atomtypes%parts(atomtypes%num_parts))
   allocate(atomtypes%itemdir1(num_atoms1))
   allocate(atomtypes%itemdir2(num_atoms2))

   ! Copy item directories
   atomtypes%itemdir1 = itemdir1_temp
   atomtypes%itemdir2 = itemdir2_temp

   ! Build each partition
   do i = 1, atomtypes%num_parts
      ! Set sizes
      atomtypes%parts(i)%num_items1 = part_num_items1(i)
      atomtypes%parts(i)%num_items2 = part_num_items2(i)
      atomtypes%parts(i)%num_children = 0

      ! Allocate and copy items
      allocate(atomtypes%parts(i)%items1(part_num_items1(i)))
      allocate(atomtypes%parts(i)%items2(part_num_items2(i)))
      atomtypes%parts(i)%items1 = part_items1(1:part_num_items1(i), i)
      atomtypes%parts(i)%items2 = part_items2(1:part_num_items2(i), i)

      ! Allocate empty arrays for neighbors and children
      allocate(atomtypes%parts(i)%signature(0))
      allocate(atomtypes%parts(i)%children(0))
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
