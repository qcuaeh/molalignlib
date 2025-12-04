! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module assorting
use parameters
use types_basic
use options
use chemistry
use molecule
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

abstract interface
   logical function compare_atoms_interface(item, elnum, typeid)
      import atomtype_item_t
      type(atomtype_item_t), intent(in) :: item
      integer, intent(in) :: elnum, typeid
   end function
end interface

procedure(compare_atoms_interface), pointer :: compare_atoms

contains

logical function compare_atoms_simple(item, elnum, typeid) result(equal)
   type(atomtype_item_t), intent(in) :: item
   integer, intent(in) :: elnum, typeid

   equal = item%elnum == elnum
end function

logical function compare_atoms_labeled(item, elnum, typeid) result(equal)
   type(atomtype_item_t), intent(in) :: item
   integer, intent(in) :: elnum, typeid

   if (item%elnum == elnum) then
      if (item%typeid == typeid) then
         equal = .TRUE.
         return
      end if
   end if

   equal = .FALSE.
end function

subroutine add_atomtype(atomtypetable, elnum, typeid, partidx)
   type(atomtype_table_t), intent(inout) :: atomtypetable
   integer, intent(in) :: elnum, typeid, partidx

   atomtypetable%num_items = atomtypetable%num_items + 1
   atomtypetable%items(atomtypetable%num_items)%elnum = elnum
   atomtypetable%items(atomtypetable%num_items)%typeid = typeid
   atomtypetable%items(atomtypetable%num_items)%partidx = partidx
end subroutine

function find_atomtype(atomtypetable, elnum, typeid) result(partidx)
   type(atomtype_table_t), intent(in) :: atomtypetable
   integer, intent(in) :: elnum, typeid
   integer :: partidx
   integer :: i

   do i = 1, atomtypetable%num_items
      if (compare_atoms(atomtypetable%items(i), elnum, typeid)) then
         partidx = atomtypetable%items(i)%partidx
         return
      end if
   end do

   partidx = 0  ! Not found
end function

subroutine collect_atomtypes(atomset1, atomset2, atoms1, atoms2, atomtypes)
! Single-pass approach using reverse mapping - no large 2D arrays needed
   integer, dimension(:), intent(in) :: atomset1, atomset2
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_t), intent(out) :: atomtypes

   ! Local variables
   type(atomtype_table_t) :: atomtypetable
   ! Small temporary arrays - only O(max_parts) size
   integer, dimension(:), allocatable :: part_count1, part_count2
   integer, dimension(:), allocatable :: part_fill1, part_fill2
   ! Item directories - O(n_atoms) size, unavoidable
   integer, dimension(:), allocatable :: itemdir1_temp, itemdir2_temp
   integer :: n_atoms1, n_atoms2
   integer :: current_part, max_parts
   integer :: atomidx, partidx, i

   if (label_flag) then
      compare_atoms => compare_atoms_labeled
   else
      compare_atoms => compare_atoms_simple
   end if

   n_atoms1 = size(atoms1)
   n_atoms2 = size(atoms2)
   max_parts = n_atoms1 + n_atoms2  ! Maximum possible partitions

   ! Allocate temporary storage
   allocate(part_count1(max_parts))
   allocate(part_count2(max_parts))
   allocate(atomtypetable%items(max_parts))
   allocate(itemdir1_temp(n_atoms1))
   allocate(itemdir2_temp(n_atoms2))

   ! Initialize
   part_count1 = 0
   part_count2 = 0
   atomtypetable%num_items = 0
   current_part = 0

   ! SINGLE PASS: Process all atoms, build assignments AND count sizes
   ! First molecule
   do i = 1, size(atomset1)
      atomidx = atomset1(i)
      partidx = find_atomtype(atomtypetable, atoms1(atomidx)%elnum, atoms1(atomidx)%typeid)
      if (partidx == 0) then
         current_part = current_part + 1
         call add_atomtype(atomtypetable, atoms1(atomidx)%elnum, atoms1(atomidx)%typeid, current_part)
         partidx = current_part
      end if
      itemdir1_temp(atomidx) = partidx
      part_count1(partidx) = part_count1(partidx) + 1
   end do

   ! Second molecule
   do i = 1, size(atomset2)
      atomidx = atomset2(i)
      partidx = find_atomtype(atomtypetable, atoms2(atomidx)%elnum, atoms2(atomidx)%typeid)
      if (partidx == 0) then
         current_part = current_part + 1
         call add_atomtype(atomtypetable, atoms2(atomidx)%elnum, atoms2(atomidx)%typeid, current_part)
         partidx = current_part
      end if
      itemdir2_temp(atomidx) = partidx
      part_count2(partidx) = part_count2(partidx) + 1
   end do

   ! Now allocate final structure with exact sizes (no waste!)
   atomtypes%num_parts = current_part
   allocate(atomtypes%parts(atomtypes%num_parts))
   allocate(atomtypes%itemdir1(n_atoms1))
   allocate(atomtypes%itemdir2(n_atoms2))

   ! Copy item directories
   atomtypes%itemdir1 = itemdir1_temp
   atomtypes%itemdir2 = itemdir2_temp

   ! Allocate each partition with exact size
   do i = 1, atomtypes%num_parts
      atomtypes%parts(i)%elnum = atomtypetable%items(i)%elnum
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
   do i = 1, n_atoms1
      if (any(atomset1 == i)) then
         partidx = itemdir1_temp(i)
         part_fill1(partidx) = part_fill1(partidx) + 1
         atomtypes%parts(partidx)%items1(part_fill1(partidx)) = i
      end if
   end do

   ! Fill second molecule using reverse mapping
   do i = 1, n_atoms2
      if (any(atomset2 == i)) then
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
