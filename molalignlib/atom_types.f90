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
! Single-pass approach using reverse mapping - no large 2D arrays needed
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_t), intent(out) :: atomtypes

   ! Local variables
   type(atomtype_table_t) :: atomtypetable
   integer :: i, num_atoms1, num_atoms2
   integer :: partidx, current_part
   integer :: max_parts
   ! Small temporary arrays - only O(max_parts) size
   integer, dimension(:), allocatable :: part_count1, part_count2
   integer, dimension(:), allocatable :: part_fill1, part_fill2
   ! Item directories - O(num_atoms) size, unavoidable
   integer, dimension(:), allocatable :: itemdir1_temp, itemdir2_temp

   num_atoms1 = size(atoms1)
   num_atoms2 = size(atoms2)
   max_parts = num_atoms1 + num_atoms2  ! Maximum possible partitions

   ! Allocate temporary storage
   allocate(part_count1(max_parts))
   allocate(part_count2(max_parts))
   allocate(atomtypetable%items(max_parts))
   allocate(itemdir1_temp(num_atoms1))
   allocate(itemdir2_temp(num_atoms2))

   ! Initialize
   part_count1 = 0
   part_count2 = 0
   atomtypetable%num_items = 0
   current_part = 0

   ! SINGLE PASS: Process all atoms, build assignments AND count sizes
   ! First molecule
   do i = 1, num_atoms1
      if (atoms1(i)%mask) then
         partidx = find_atomtype(atomtypetable, atoms1(i))
         if (partidx == 0) then
            current_part = current_part + 1
            call add_atomtype(atomtypetable, atoms1(i), current_part)
            partidx = current_part
         end if
         itemdir1_temp(i) = partidx
         part_count1(partidx) = part_count1(partidx) + 1
      end if
   end do

   ! Second molecule
   do i = 1, num_atoms2
      if (atoms2(i)%mask) then
         partidx = find_atomtype(atomtypetable, atoms2(i))
         if (partidx == 0) then
            current_part = current_part + 1
            call add_atomtype(atomtypetable, atoms2(i), current_part)
            partidx = current_part
         end if
         itemdir2_temp(i) = partidx
         part_count2(partidx) = part_count2(partidx) + 1
      end if
   end do

   ! Now allocate final structure with exact sizes (no waste!)
   atomtypes%num_parts = current_part
   allocate(atomtypes%parts(atomtypes%num_parts))
   allocate(atomtypes%itemdir1(num_atoms1))
   allocate(atomtypes%itemdir2(num_atoms2))

   ! Copy item directories
   atomtypes%itemdir1 = itemdir1_temp
   atomtypes%itemdir2 = itemdir2_temp

   ! Allocate each partition with exact size
   do i = 1, atomtypes%num_parts
      atomtypes%parts(i)%num_items1 = part_count1(i)
      atomtypes%parts(i)%num_items2 = part_count2(i)
      atomtypes%parts(i)%num_children = 0

      allocate(atomtypes%parts(i)%items1(part_count1(i)))
      allocate(atomtypes%parts(i)%items2(part_count2(i)))
      allocate(atomtypes%parts(i)%signature(0))
      allocate(atomtypes%parts(i)%children(0))
   end do

   ! FAST FILL: Use assignments to populate final arrays efficiently
   allocate(part_fill1(current_part))
   allocate(part_fill2(current_part))
   part_fill1 = 0  ! Current fill position for each partition
   part_fill2 = 0

   ! Fill first molecule using reverse mapping
   do i = 1, num_atoms1
      if (atoms1(i)%mask) then
         partidx = itemdir1_temp(i)
         part_fill1(partidx) = part_fill1(partidx) + 1
         atomtypes%parts(partidx)%items1(part_fill1(partidx)) = i
      end if
   end do

   ! Fill second molecule using reverse mapping
   do i = 1, num_atoms2
      if (atoms2(i)%mask) then
         partidx = itemdir2_temp(i)
         part_fill2(partidx) = part_fill2(partidx) + 1
         atomtypes%parts(partidx)%items2(part_fill2(partidx)) = i
      end if
   end do

   ! Clean up small temporary arrays
   deallocate(part_count1, part_count2)
   deallocate(itemdir1_temp, itemdir2_temp)
   deallocate(part_fill1, part_fill2)
end subroutine

end module
