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
use common_types
use flags
use chemdata
use molecule
implicit none
private
public collect_atomtypes

type :: atomtype_item_t
   integer(ik) :: elnum
   integer(ik) :: group
   integer(ik) :: partidx
end type

type :: atomtype_table_t
   integer(ik) :: n_items
   type(atomtype_item_t), dimension(:), allocatable :: items
end type

abstract interface
   logical(lk) function compare_atoms_interface(item, elnum, group)
      use parameters
      import atomtype_item_t
      type(atomtype_item_t), intent(in) :: item
      integer(ik), intent(in) :: elnum, group
   end function
end interface

procedure(compare_atoms_interface), pointer :: compare_atoms

contains

logical(lk) function compare_atoms_simple(item, elnum, group) result(equal)
   type(atomtype_item_t), intent(in) :: item
   integer(ik), intent(in) :: elnum, group

   equal = item%elnum == elnum
end function

logical(lk) function compare_atoms_labeled(item, elnum, group) result(equal)
   type(atomtype_item_t), intent(in) :: item
   integer(ik), intent(in) :: elnum, group

   if (item%elnum == elnum) then
      if (item%group == group) then
         equal = .TRUE.
         return
      end if
   end if

   equal = .FALSE.
end function

subroutine add_atomtype(atomtypetable, elnum, group, partidx)
   type(atomtype_table_t), intent(inout) :: atomtypetable
   integer(ik), intent(in) :: elnum, group, partidx

   atomtypetable%n_items = atomtypetable%n_items + 1
   atomtypetable%items(atomtypetable%n_items)%elnum = elnum
   atomtypetable%items(atomtypetable%n_items)%group = group
   atomtypetable%items(atomtypetable%n_items)%partidx = partidx
end subroutine

function find_atomtype(atomtypetable, elnum, group) result(partidx)
   type(atomtype_table_t), intent(in) :: atomtypetable
   integer(ik), intent(in) :: elnum, group
   integer(ik) :: partidx
   integer(ik) :: i

   do i = 1, atomtypetable%n_items
      if (compare_atoms(atomtypetable%items(i), elnum, group)) then
         partidx = atomtypetable%items(i)%partidx
         return
      end if
   end do

   partidx = 0  ! Not found
end function

subroutine collect_atomtypes(atoms1, atoms2, atomtypes)
! Partition the atoms of both molecules by atom type. Every atom passed in
! is classified, so pass only the included atoms (e.g. atoms(atomset)):
! items and itemdirs then use the compact numbering of the included atoms.
! Single-pass approach using reverse mapping - no large 2D arrays needed
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_t), intent(out) :: atomtypes

   ! Local variables
   type(atomtype_table_t) :: atomtypetable
   ! Small temporary arrays - only O(max_parts) size
   integer(ik), dimension(:), allocatable :: part_count1, part_count2
   integer(ik), dimension(:), allocatable :: part_fill1, part_fill2
   integer(ik) :: n_atoms1, n_atoms2
   integer(ik) :: current_part, max_parts
   integer(ik) :: atomidx, partidx, i

   if (atomlabel_flag) then
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
   allocate(atomtypes%itemdir1(n_atoms1))
   allocate(atomtypes%itemdir2(n_atoms2))

   ! Initialize counters
   part_count1 = 0
   part_count2 = 0
   atomtypetable%n_items = 0
   current_part = 0

   ! SINGLE PASS: Process all atoms, build assignments AND count sizes
   ! First molecule
   do atomidx = 1, n_atoms1
      partidx = find_atomtype(atomtypetable, atoms1(atomidx)%elnum, atoms1(atomidx)%group)
      if (partidx == 0) then
         current_part = current_part + 1
         call add_atomtype(atomtypetable, atoms1(atomidx)%elnum, atoms1(atomidx)%group, current_part)
         partidx = current_part
      end if
      atomtypes%itemdir1(atomidx) = partidx
      part_count1(partidx) = part_count1(partidx) + 1
   end do

   ! Second molecule
   do atomidx = 1, n_atoms2
      partidx = find_atomtype(atomtypetable, atoms2(atomidx)%elnum, atoms2(atomidx)%group)
      if (partidx == 0) then
         current_part = current_part + 1
         call add_atomtype(atomtypetable, atoms2(atomidx)%elnum, atoms2(atomidx)%group, current_part)
         partidx = current_part
      end if
      atomtypes%itemdir2(atomidx) = partidx
      part_count2(partidx) = part_count2(partidx) + 1
   end do

   ! Now allocate final structure with exact sizes (no waste!)
   atomtypes%n_parts = current_part
   allocate(atomtypes%parts(atomtypes%n_parts))

   ! Allocate each partition with exact size
   do i = 1, atomtypes%n_parts
      atomtypes%parts(i)%elnum = atomtypetable%items(i)%elnum
      atomtypes%parts(i)%n_items1 = part_count1(i)
      atomtypes%parts(i)%n_items2 = part_count2(i)
      atomtypes%parts(i)%n_children = 0

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
      partidx = atomtypes%itemdir1(i)
      part_fill1(partidx) = part_fill1(partidx) + 1
      atomtypes%parts(partidx)%items1(part_fill1(partidx)) = i
   end do

   ! Fill second molecule using reverse mapping
   do i = 1, n_atoms2
      partidx = atomtypes%itemdir2(i)
      part_fill2(partidx) = part_fill2(partidx) + 1
      atomtypes%parts(partidx)%items2(part_fill2(partidx)) = i
   end do

   ! Clean up small temporary arrays
   deallocate(part_count1, part_count2)
   deallocate(part_fill1, part_fill2)
end subroutine

end module
