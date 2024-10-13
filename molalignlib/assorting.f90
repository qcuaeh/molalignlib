! MolAlignLib
! Copyright (C) 2022 José M. Vásquez

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
use stdio
use kinds
use partition
use permutation
use molecule
use flags
use bounds
use sorting
use chemdata

implicit none

contains

! Partition atoms by atomic number
subroutine compute_crosseltypes(mol0, mol1)
   type(molecule_type), intent(inout) :: mol0, mol1
   ! Local variables
   integer :: i, elnum, typeidx
   integer, allocatable :: typedict(:)
   type(partition_type) :: eltypes0, eltypes1

   call eltypes0%allocate(size(mol0%atoms))
   call eltypes1%allocate(size(mol1%atoms))

   allocate (typedict(nelem))
   typedict(:) = 0

   do i = 1, size(mol0%atoms)
      elnum = mol0%atoms(i)%elnum
      if (typedict(elnum) == 0) then
         typeidx = eltypes0%new_part()
         call eltypes1%add_part(typeidx)
         typedict(elnum) = typeidx
      end if
      call eltypes0%append_to_part(typedict(elnum), i)
   end do

   do i = 1, size(mol1%atoms)
      elnum = mol1%atoms(i)%elnum
      if (typedict(elnum) == 0) then
         typeidx = eltypes1%new_part()
         call eltypes0%add_part(typeidx)
         typedict(elnum) = typeidx
      end if
      call eltypes1%append_to_part(typedict(elnum), i)
   end do

   call mol0%set_eltypes(eltypes0)
   call mol1%set_eltypes(eltypes1)

end subroutine

! Iteratively compute MNA types at all levels
subroutine compute_mnatypes(mol)
   type(molecule_type), intent(inout) :: mol
   ! Local variables
   type(partition_type) :: intypes, types
!   integer h

   intypes = mol%eltypes

   do
      ! Compute MNA types at next level
      call compute_nextmnatypes(mol%atoms, intypes, types)
      ! Exit the loop if types are unchanged
      if (types == intypes) exit
      intypes = types
   end do

   call mol%set_mnatypes(types)
!   do h = 1 , size(types%parts)
!      print *, types%parts(h)%indices
!   end do

end subroutine

! Assort atoms by MNA type at next level
subroutine compute_nextmnatypes(atoms, intypes, types)
   type(atom_type), dimension(:), intent(in) :: atoms
   type(partition_type), intent(in) :: intypes
   type(partition_type), intent(out) :: types
   ! Local variables
   integer :: h, k, i
   integer :: atomidx, typeidx
   integer, allocatable :: adjlist(:)
   integer, allocatable :: atomintypes(:)
   integer, allocatable :: adjtypecounts(:)
   type(partitionlist_type) :: subtypelist
   type(countlist_type), allocatable :: adjtypeaccounting(:)

   call types%allocate(size(atoms))
   allocate (adjtypeaccounting(size(atoms)))
   atomintypes = intypes%get_item_types()

   do h = 1, size(intypes%parts)
      call subtypelist%allocate(size(atoms))
      outer: do i = 1, size(intypes%parts(h)%indices)
         atomidx = intypes%parts(h)%indices(i)
         adjlist = atoms(atomidx)%adjlist
         adjtypecounts = typecounts(atomintypes, adjlist)
         do k = 1, size(subtypelist%indices)
            if (all(adjtypecounts == adjtypeaccounting(subtypelist%indices(k))%counts)) then
               call types%append_to_part(subtypelist%indices(k), atomidx)
               cycle outer
            end if
         end do
         typeidx = types%new_part()
         call subtypelist%append(typeidx)
         call types%append_to_part(typeidx, atomidx)
         adjtypeaccounting(typeidx)%counts = typecounts(atomintypes, adjlist)
      end do outer
   end do

end subroutine

function typecounts(atomtypes, adjlist)
   integer, intent(in) :: atomtypes(:)
   integer, intent(in) :: adjlist(:)
   integer, allocatable :: typecounts(:)
   ! Local variables
   integer :: i

   allocate(typecounts(size(atomtypes)))
   typecounts(:) = 0

   do i = 1, size(adjlist)
      typecounts(atomtypes(adjlist(i))) = typecounts(atomtypes(adjlist(i))) + 1
   end do

end function

! Find next level cross MNA types between mol0 and mol1
subroutine compute_crossmnatypes(adjlists0, adjlists1, nintype, intypes0, intypes1, &
      ntype, types0, types1)
   type(indexlist_type), dimension(:), intent(in) :: adjlists0, adjlists1
   integer, intent(in) :: nintype
   integer, dimension(:), intent(in) :: intypes0, intypes1
   integer, intent(out) :: ntype
   integer, dimension(:), intent(out) :: types0, types1
   ! Local variables
   integer :: h, i, j
   logical :: untyped(size(adjlists0))
   integer :: archeatom(size(adjlists0))

   ntype = 0
   untyped(:) = .true.

   do i = 1, size(adjlists0)
      if (untyped(i)) then
         ntype = ntype + 1
         types0(i) = ntype
         archeatom(ntype) = i
         do j = i + 1, size(adjlists0)
            if (untyped(j)) then
               if (intypes0(j) == intypes0(i)) then
                  if (same_adjacency(nintype, intypes0, adjlists0(i)%indices, intypes0, adjlists0(j)%indices)) then
                     types0(j) = ntype
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

   untyped(:) = .true.

   do h = 1, ntype
      do i = 1, size(adjlists0)
         if (untyped(i)) then
            if (intypes1(i) == intypes0(archeatom(h))) then
               if (same_adjacency(nintype, intypes0, adjlists0(archeatom(h))%indices, intypes1, adjlists1(i)%indices)) then
                  types1(i) = h
                  untyped(i) = .false.
               end if
            end if
         end if
      end do
   end do

   do i = 1, size(adjlists0)
      if (untyped(i)) then
         ntype = ntype + 1
         types1(i) = ntype
         do j = i + 1, size(adjlists0)
            if (untyped(j)) then
               if (intypes1(j) == intypes1(i)) then
                  if (same_adjacency(nintype, intypes1, adjlists1(i)%indices, intypes1, adjlists1(j)%indices)) then
                     types1(j) = ntype
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

end subroutine

! Test if adjacent atoms to a pair of atoms in mol0 and mol1 are the same
function same_adjacency(neltype, atomtype0, adjlist0, atomtype1, adjlist1) result(sameadj)
   integer, intent(in) :: neltype
   integer, dimension(:), intent(in) :: adjlist0, adjlist1
   integer, dimension(:) :: atomtype0, atomtype1
   ! Result variable
   logical :: sameadj
   ! Local variables
   integer :: i0, i1
   integer, dimension(neltype) :: n0, n1

   sameadj = .true.

   if (size(adjlist0) /= size(adjlist1)) then
      sameadj = .false.
      return
   end if

!  Find if adjacent atoms are the same

   n0(:) = 0
   n1(:) = 0

   do i0 = 1, size(adjlist0)
      n0(atomtype0(adjlist0(i0))) = n0(atomtype0(adjlist0(i0))) + 1
   end do

   do i1 = 1, size(adjlist1)
      n1(atomtype1(adjlist1(i1))) = n1(atomtype1(adjlist1(i1))) + 1
   end do

   if (any(n0 /= n1)) then
      sameadj = .false.
      return
   end if

end function

end module
