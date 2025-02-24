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

module bipartitioning
use parameters
use globals
use sorting
use molecule
use tupledict
use bipartition
use partitiondict
use permutation

implicit none

contains

! Partition atoms by atomic number and label
subroutine compute_crosseltypes(mol1, mol2, eltypes)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(out) :: eltypes
   ! Local variables
   integer :: eltype(2)
   type(tupledict_type) :: typedict
   type(bipartpointer_type), allocatable :: typelist(:)
   integer :: i, num_atoms1, num_atoms2

   num_atoms1 = size(mol1%atoms)
   num_atoms2 = size(mol2%atoms)
   call eltypes%initialize(num_atoms1, num_atoms2)
   call typedict%initialize(num_atoms1 + num_atoms2, 'ordered')
   allocate (typelist(typedict%num_slots))

   do i = 1, num_atoms1
      eltype(1) = mol1%atoms(i)%elnum
      eltype(2) = mol1%atoms(i)%label
      if (.not. (eltype .in. typedict)) then
         typelist(typedict%new_index(eltype))%ptr => &
            eltypes%new_part(num_atoms1, num_atoms2)
      end if
      call typelist(typedict%get_index(eltype))%ptr%add1(i)
   end do

   do i = 1, num_atoms2
      eltype(1) = mol2%atoms(i)%elnum
      eltype(2) = mol2%atoms(i)%label
      if (.not. (eltype .in. typedict)) then
         typelist(typedict%new_index(eltype))%ptr => &
            eltypes%new_part(num_atoms1, num_atoms2)
      end if
      call typelist(typedict%get_index(eltype))%ptr%add2(i)
   end do

end subroutine

! Level up MNA types
subroutine levelup_crossmnatypes(mol1, mol2, mnatypes, subtypes)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(in) :: mnatypes
   type(bipartition_type), intent(out) :: subtypes
   ! Local variables
   integer :: h, i, iatom
   type(tupledict_type) :: typedict
   type(bipartpointer_type), allocatable :: typelist(:)
   integer, allocatable :: neighborhood(:)

   call subtypes%initialize(mnatypes%num_items1, mnatypes%num_items2)
   call typedict%initialize(mnatypes%largest_part_size, 'unordered')
   allocate (typelist(typedict%num_slots))

   do h = 1, mnatypes%num_parts

      do i = 1, mnatypes%parts(h)%part_size1
         iatom = mnatypes%parts(h)%items1(i)
         neighborhood = mnatypes%idcs1(mol1%atoms(iatom)%adjlist)
         if (.not. (neighborhood .in. typedict)) then
            typelist(typedict%new_index(neighborhood))%ptr => &
               subtypes%new_part(mnatypes%parts(h)%part_size1, mnatypes%parts(h)%part_size2)
         end if
         call typelist(typedict%get_index(neighborhood))%ptr%add1(iatom)
      end do

      do i = 1, mnatypes%parts(h)%part_size2
         iatom = mnatypes%parts(h)%items2(i)
         neighborhood = mnatypes%idcs2(mol2%atoms(iatom)%adjlist)
         if (.not. (neighborhood .in. typedict)) then
            typelist(typedict%new_index(neighborhood))%ptr => &
               subtypes%new_part(mnatypes%parts(h)%part_size1, mnatypes%parts(h)%part_size2)
         end if
         call typelist(typedict%get_index(neighborhood))%ptr%add2(iatom)
      end do

      call typedict%reset()

   end do

end subroutine

subroutine compute_crossmnatypes(mol1, mol2, mnatypes)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(inout) :: mnatypes
   ! Local variables
   type(bipartition_type) :: subtypes

   do

      ! Compute MNA upper level types
      call levelup_crossmnatypes(mol1, mol2, mnatypes, subtypes)

      ! Exit loop if types did not change
      if (subtypes == mnatypes) then
         mnatypes = subtypes
         exit
      end if

      ! Update mnatypes
      mnatypes = subtypes

   end do

end subroutine

subroutine split_crossmnatypes(h, perm, mnatypes)
   integer, intent(in) :: h
   integer, intent(in) :: perm(:)
   type(bipartition_type), intent(inout) :: mnatypes
   ! Local variables
   integer :: i, k
   type(bipart_type), pointer :: newtype
   type(bipartition_type) :: subtypes

   if (mnatypes%parts(h)%part_size1 /= mnatypes%parts(h)%part_size2) then
      call mnatypes%print_parts()
      write (stderr, *) mnatypes%parts(h)%part_size1, mnatypes%parts(h)%part_size2
      error stop 'part_size1 /= part_size2'
   end if

   call subtypes%initialize(mnatypes%num_items1, mnatypes%num_items2)

   do k = 1, h - 1
      call subtypes%add_part(mnatypes%parts(k))
   end do

   do i = 1, mnatypes%parts(h)%part_size1
      newtype => subtypes%new_part(mnatypes%parts(h)%part_size1, mnatypes%parts(h)%part_size2)
      call newtype%add1(mnatypes%parts(h)%items1(i))
      call newtype%add2(mnatypes%parts(h)%items2(perm(i)))
   end do

   do k = h + 1, mnatypes%num_parts
      call subtypes%add_part(mnatypes%parts(k))
   end do

   mnatypes = subtypes

end subroutine

end module
