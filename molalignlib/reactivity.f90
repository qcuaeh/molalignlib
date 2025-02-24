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
use parameters
use sorting
use permutation
use alignment
use molecule
use tracking
use common_types
use lcrs_tree

implicit none

contains

subroutine remove_reactive_bonds(mol1, mol2, molfrags1, molfrags2, mnatypes, atomperm)
!
! Remove mismatched and reactive bonds
!

   ! Arguments
   type(mol_type), intent(inout) :: mol1, mol2
   type(intlist_type), dimension(:) :: molfrags1, molfrags2
   type(bipartition_container), intent(in) :: mnatypes
   integer, dimension(:), intent(in) :: atomperm

   ! Local variables
   type(atom_type), dimension(:), allocatable :: atoms1, atoms2
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
!   integer, dimension(:), allocatable :: adjidcs1, adjidcs2
   integer, dimension(:), allocatable :: invatomperm
   integer :: iatom
   integer :: j, jatom
!   integer :: k, katom

   ! Initialization

   atoms1 = mol1%atoms
   atoms2 = mol2%atoms
   adjmat1 = mol1%adjmat
   adjmat2 = mol2%adjmat
   invatomperm = inverse_perm(atomperm)

   ! Remove mismatched bonds

   do iatom = 1, size(atoms1)
      do j = 1, size(atoms1(iatom)%adjlist)
         jatom = atoms1(iatom)%adjlist(j)
         if (.not. adjmat2(atomperm(iatom), atomperm(jatom))) then
!            write (stderr, *) 'remove mol1 bond:', iatom, jatom
            call mol1%remove_bond(iatom, jatom)
!            adjidcs1 = mnatypes%parts(mnatypes%itemdir1(jatom))%indices1
!            do k = 1, size(adjidcs1)
!               katom = adjidcs1(k)
!               call mol1%remove_bond(iatom, katom)
!               call mol2%remove_bond(atomperm(iatom), atomperm(katom))
!            end do
         end if
      end do
   end do

   do iatom = 1, size(atoms2)
      do j = 1, size(atoms2(iatom)%adjlist)
         jatom = atoms2(iatom)%adjlist(j)
         if (.not. adjmat1(invatomperm(iatom), invatomperm(jatom))) then
!            write (stderr, *) 'remove mol2 bond:', iatom, jatom
            call mol2%remove_bond(iatom, jatom)
!            adjidcs2 = mnatypes%parts(mnatypes%itemdir2(jatom))%indices2
!            do k = 1, size(adjidcs2)
!               katom = adjidcs2(k)
!               call mol2%remove_bond(iatom, katom)
!               call mol1%remove_bond(invatomperm(iatom), invatomperm(katom))
!            end do
         end if
      end do
   end do

   ! Dissociate water molecules
!
!   do iatom = 1, size(molfrags1)
!      if (all(sorted(atoms1(molfrags1(iatom)%n)%elnum) == [1, 1, 8])) then
!         do j = 1, size(molfrags1(iatom)%n)
!            jatom = molfrags1(iatom)%n(j)
!            do k = 1, size(atoms1(jatom)%adjlist)
!               katom = atoms1(jatom)%adjlist(k)
!               call mol1%remove_bond(jatom, katom)
!            end do
!         end do
!      end if
!   end do
!
!   do iatom = 1, size(molfrags2)
!      if (all(sorted(atoms2(molfrags2(iatom)%n)%elnum) == [1, 1, 8])) then
!         do j = 1, size(molfrags2(iatom)%n)
!            jatom = molfrags2(iatom)%n(j)
!            do k = 1, size(atoms2(jatom)%adjlist)
!               katom = atoms2(jatom)%adjlist(k)
!               call mol2%remove_bond(jatom, katom)
!            end do
!         end do
!      end if
!   end do

end subroutine

end module
