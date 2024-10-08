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

module reactivity
use kinds
use sorting
use permutation
use alignment
use partition
use molecule

implicit none

contains

subroutine remove_reactive_bonds(mol0, mol1, mapping)

   type(molecule_type), intent(inout) :: mol0, mol1
   integer, dimension(:), intent(in) :: mapping

   integer :: i, j, k, j_, k_
   integer :: natom0, natom1
   type(atompartition_type) :: mnatypes0, mnatypes1
   integer, allocatable, dimension(:) :: elnums0, elnums1
   integer, allocatable, dimension(:) :: atommnatypes0, atommnatypes1
   type(atomlist_type), allocatable, dimension(:) :: adjlists0, adjlists1
   type(atomlist_type), allocatable, dimension(:) :: molfragparts0, molfragparts1
   logical, allocatable, dimension(:, :) :: adjmat0, adjmat1
   integer, allocatable, dimension(:) :: mnatypepartidcs0, mnatypepartidcs1
   integer, allocatable, dimension(:) :: unmapping
   real(rk) :: rotquat(4)

   ! Align coordinates

   rotquat = leastrotquat(mol0%natom, mol0%get_weights(), mol0%get_coords(), mol1%get_coords(), mapping)
   call mol1%rotate_coords(rotquat)

!   write (stderr, *) sqrt(squaredist(mol0%natom, mol0%get_weights(), mol0%get_coords(), mol1%get_coords(), mapping) &
!         /sum(mol0%get_weights()))

   ! Initialization

   natom0 = mol0%get_natom()
   natom1 = mol1%get_natom()
   adjmat0 = mol0%get_adjmatrix()
   adjmat1 = mol1%get_adjmatrix()
   elnums0 = mol0%get_elnums()
   elnums1 = mol1%get_elnums()
   adjlists0 = mol0%get_adjlists()
   adjlists1 = mol1%get_adjlists()
   mnatypes0 = mol0%mnatypes
   mnatypes1 = mol1%mnatypes
   atommnatypes0 = mol0%mnatypes%get_atomtypes()
   atommnatypes1 = mol1%mnatypes%get_atomtypes()
   molfragparts0 = mol0%get_molfrags()
   molfragparts1 = mol1%get_molfrags()
   unmapping = inverse_permutation(mapping)

   ! Remove mismatched bonds

   do i = 1, natom0
      do j_ = 1, size(adjlists0(i)%atomidcs)
         j = adjlists0(i)%atomidcs(j_)
         if (.not. adjmat1(mapping(i), mapping(j))) then
            mnatypepartidcs0 = mnatypes0%subsets(atommnatypes0(j))%atomidcs
            do k_ = 1, size(mnatypepartidcs0)
               k = mnatypepartidcs0(k_)
               call mol0%remove_bond(i, k)
               call mol1%remove_bond(mapping(i), mapping(k))
            end do
         end if
      end do
   end do

   do i = 1, natom1
      do j_ = 1, size(adjlists1(i)%atomidcs)
         j = adjlists1(i)%atomidcs(j_)
         if (.not. adjmat0(unmapping(i), unmapping(j))) then
            mnatypepartidcs1 = mnatypes1%subsets(atommnatypes1(j))%atomidcs
            do k_ = 1, size(mnatypepartidcs1)
               k = mnatypepartidcs1(k_)
               call mol1%remove_bond(i, k)
               call mol0%remove_bond(unmapping(i), unmapping(k))
            end do
         end if
      end do
   end do

   ! Remove water bonds

   do i = 1, size(molfragparts0)
      if (all(sorted(elnums0(molfragparts0(i)%atomidcs)) == [1, 1, 8])) then
         do j_ = 1, size(molfragparts0(i)%atomidcs)
            j = molfragparts0(i)%atomidcs(j_)
            do k_ = 1, size(adjlists0(j)%atomidcs)
               k = adjlists0(j)%atomidcs(k_)
               call mol0%remove_bond(j, k)
            end do
         end do
      end if
   end do

   do i = 1, size(molfragparts1)
      if (all(sorted(elnums1(molfragparts1(i)%atomidcs)) == [1, 1, 8])) then
         do j_ = 1, size(molfragparts1(i)%atomidcs)
            j = molfragparts1(i)%atomidcs(j_)
            do k_ = 1, size(adjlists1(j)%atomidcs)
               k = adjlists1(j)%atomidcs(k_)
               call mol1%remove_bond(j, k)
            end do
         end do
      end if
   end do

end subroutine

end module
