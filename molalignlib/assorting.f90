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
use hashtable
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
   integer :: i, elnum
   type(partition_type) :: eltypes0, eltypes1
   type(pointertopart_type), dimension(:), allocatable :: typelist0, typelist1

   allocate (typelist0(nelem))
   allocate (typelist1(nelem))

   call eltypes0%init(size(mol0%atoms))
   call eltypes1%init(size(mol1%atoms))

   do i = 1, size(mol0%atoms)
      elnum = mol0%atoms(i)%elnum
      if (.not. associated(typelist0(elnum)%ptr)) then
         typelist0(elnum)%ptr => eltypes0%get_new_part()
         typelist1(elnum)%ptr => eltypes1%get_new_part()
      end if
      call typelist0(elnum)%ptr%add(i)
   end do

   do i = 1, size(mol1%atoms)
      elnum = mol1%atoms(i)%elnum
      if (.not. associated(typelist1(elnum)%ptr)) then
         typelist1(elnum)%ptr => eltypes1%get_new_part()
         typelist0(elnum)%ptr => eltypes0%get_new_part()
      end if
      call typelist1(elnum)%ptr%add(i)
   end do

!   call eltypes0%print_parts()
!   call eltypes1%print_parts()

   call mol0%set_eltypes(eltypes0)
   call mol1%set_eltypes(eltypes1)

end subroutine

! Level up cross MNA types
subroutine levelup_crossmnatypes(atoms0, atoms1, types0, types1, subtypes0, subtypes1)
   type(atom_type), dimension(:), intent(in) :: atoms0, atoms1
   type(partition_type), intent(in) :: types0, types1
   type(partition_type), intent(out) :: subtypes0, subtypes1
   ! Local variables
   integer :: h, i, index_i
   integer :: total_atoms, max_part_size
   type(hashtable_type) :: subtypedict
   type(pointertopart_type), dimension(:), allocatable :: subtypelist0, subtypelist1
   integer, allocatable :: neighborhood(:)

   total_atoms = size(atoms0) + size(atoms1)
   max_part_size = max(types0%max_part_size, types1%max_part_size)

   call subtypes0%init(total_atoms)
   call subtypes1%init(total_atoms)
   call subtypedict%init(max_part_size)

   allocate (subtypelist0(subtypedict%size))
   allocate (subtypelist1(subtypedict%size))

   do h = 1, types0%size
      do i = 1, types0%parts(h)%size
         index_i = types0%parts(h)%indices(i)
         neighborhood = types0%index_part_map(atoms0(index_i)%adjlist)
         if (.not. subtypedict%has_index(neighborhood)) then
            subtypelist0(subtypedict%get_new_index(neighborhood))%ptr => subtypes0%get_new_part()
            subtypelist1(subtypedict%get_new_index(neighborhood))%ptr => subtypes1%get_new_part()
         end if
         call subtypelist0(subtypedict%get_index(neighborhood))%ptr%add(index_i)
      end do
      do i = 1, types1%parts(h)%size
         index_i = types1%parts(h)%indices(i)
         neighborhood = types1%index_part_map(atoms1(index_i)%adjlist)
         if (.not. subtypedict%has_index(neighborhood)) then
            subtypelist1(subtypedict%get_new_index(neighborhood))%ptr => subtypes1%get_new_part()
            subtypelist0(subtypedict%get_new_index(neighborhood))%ptr => subtypes0%get_new_part()
         end if
         call subtypelist1(subtypedict%get_index(neighborhood))%ptr%add(index_i)
      end do
      call subtypedict%reset()
   end do

end subroutine

! Iteratively level up cross MNA types
subroutine compute_crossmnatypes2(mol0, mol1)
   type(molecule_type), intent(inout) :: mol0, mol1
   ! Local variables
   type(partition_type) :: mnatypes0, mnatypes1
   type(partition_type) :: mnasubtypes0, mnasubtypes1

   mnatypes0 = mol0%eltypes
   mnatypes1 = mol1%eltypes

   do
      ! Compute next level MNA types
      call levelup_crossmnatypes(mol0%atoms, mol1%atoms, mnatypes0, mnatypes1, mnasubtypes0, mnasubtypes1)
      ! Exit the loop if types are unchanged
      if (mnasubtypes0 == mnatypes0 .and. mnasubtypes1 == mnatypes1) exit
      mnatypes0 = mnasubtypes0
      mnatypes1 = mnasubtypes1
   end do

   call mnatypes0%print_parts()
   call mnatypes1%print_parts()

!   call mol0%set_mnatypes(mnatypes0)
!   call mol1%set_mnatypes(mnatypes1)

end subroutine

! Level up MNA types
subroutine levelup_mnatypes(atoms, types, subtypes)
   type(atom_type), dimension(:), intent(in) :: atoms
   type(partition_type), intent(in) :: types
   type(partition_type), intent(out) :: subtypes
   ! Local variables
   integer :: h, i, index_i
   type(hashtable_type) :: subtypedict
   type(pointertopart_type), allocatable :: subtypelist(:)
   integer, allocatable :: neighborhood(:)

   call subtypes%init(size(atoms))
   call subtypedict%init(types%max_part_size)
   allocate (subtypelist(subtypedict%size))

   do h = 1, types%size
      do i = 1, types%parts(h)%size
         index_i = types%parts(h)%indices(i)
         neighborhood = types%index_part_map(atoms(index_i)%adjlist)
         if (.not. subtypedict%has_index(neighborhood)) then
            subtypelist(subtypedict%get_new_index(neighborhood))%ptr => subtypes%get_new_part()
         end if
         call subtypelist(subtypedict%get_index(neighborhood))%ptr%add(index_i)
      end do
      call subtypedict%reset()
   end do

end subroutine

! Iteratively level up MNA types
subroutine compute_mnatypes(mol)
   type(molecule_type), intent(inout) :: mol
   ! Local variables
   type(partition_type) :: mnatypes, mnasubtypes

   mnatypes = mol%eltypes

   do
      ! Compute next level MNA mnasubtypes
      call levelup_mnatypes(mol%atoms, mnatypes, mnasubtypes)
      ! Exit the loop if mnasubtypes are unchanged
      if (mnasubtypes == mnatypes) exit
      mnatypes = mnasubtypes
   end do

   call mnatypes%print_parts()
!   call mol%set_mnatypes(mnatypes)

end subroutine

! Compute next level cross MNA types between mol0 and mol1
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
