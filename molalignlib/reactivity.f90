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
use common_types
use globals
use random
use rigid_body
use rotation
use adjacency
use permutation
use alignment
use assignment
use backtracking
use biasing
use molecule
use tracking
use lcrs_tree
use partitioning
use registry

implicit none

contains

subroutine remove_reactive_bonds( mol1, mol2, eltypes, atomperm)
   type(mol_type), intent(inout) :: mol1, mol2
   type(bipartition_container), intent(in) :: eltypes
   integer, dimension(:), intent(out) :: atomperm

   ! Local variables
   type(registry_type), target :: results
   type(tree_node), pointer :: mnatree
   type(bipartition_container) :: mnatypes
   type(intlist_type), dimension(:), allocatable :: molfrags1, molfrags2
   type(intmatrix_type), allocatable :: biases(:)
   integer, dimension(:), allocatable :: auxperm, invatomperm
   real(rk), dimension(:,:), allocatable :: coords1, coords2
   real(rk) :: eigquat(4), totquat(4)
   integer :: adjd, num_atoms1, num_trials, num_steps
   integer :: j, iatom, jatom
   integer, pointer :: lead_count
!   integer, dimension(:), allocatable :: indices1, indices2
!   integer :: k, katom

   num_atoms1 = size(mol1%atoms)
   coords1 = mol1%get_weighted_coords()
   coords2 = mol2%get_weighted_coords()
   call results%initialize(max_records)
   allocate (auxperm(num_atoms1))

   ! Compute MNA types
   call tree_from_partition( eltypes, mnatree)
   call compute_consistent_mnatypes( mol1, mol2, mnatree)
   call partition_from_tree( mnatree, mnatypes)

   ! Find molecular fragments
   call find_molfrags( mol1, first_partition(eltypes), molfrags1)
   call find_molfrags( mol2, second_partition(eltypes), molfrags2)

   ! Mirror coordinates
   if (mirror_flag) then
      call mirror_coords( coords2)
   end if

   ! Translate atoms to their centroids
   call translate_coords( coords1, -centroid(coords1))
   call translate_coords( coords2, -centroid(coords2))

   ! Find unfeasible assignments
   call compute_mna_biases( mol1, mol2, eltypes, biases)

   ! Initialize random number generator
   call random_initialize()

   ! Find reactive bonds

   num_trials = 0
   lead_count => results%records(1)%count

   do while (lead_count < max_count .and. num_trials < max_trials)

      num_trials = num_trials + 1

      ! Aply a random rotation to coords2
      call rotate_coords(coords2, randrotquat())

      ! Assign atoms with current orientation
      call assign_atoms_biased(eltypes, coords1, coords2, biases, atomperm)
      totquat = leasteigquat(atomperm, coords1, coords2)
      call rotate_coords(coords2, totquat)
      num_steps = 1

      do while (iter_flag)
         call assign_atoms_biased(eltypes, coords1, coords2, biases, auxperm)
         if (all(auxperm == atomperm)) exit
         atomperm = auxperm
         eigquat = leasteigquat(atomperm, coords1, coords2)
         call rotate_coords(coords2, eigquat)
         totquat = quatmul(eigquat, totquat)
         num_steps = num_steps + 1
      end do

      call minadjdiff( eltypes, mnatypes, molfrags1, mol1, mol2, coords1, coords2, atomperm)

      ! Update results
      adjd = adjacencydiff(atomperm, mol1%adjmat, mol2%adjmat)
      call results%push_adjd(atomperm, num_steps, angle(totquat), adjd)

   end do

   write (stderr, *) 'adjd before', adjacencydiff( identity_perm(num_atoms1), mol1%adjmat, mol2%adjmat)
   write (stderr, *) 'adjd after', adjacencydiff( results%records(1)%atomperm, mol1%adjmat, mol2%adjmat)

   ! Remove reactive bonds

   atomperm = results%records(1)%atomperm
   invatomperm = inverse_perm(atomperm)

   ! Remove mismatched bonds
   do iatom = 1, size(mol1%atoms)
      do j = 1, size(mol1%atoms(iatom)%adjlist)
         jatom = mol1%atoms(iatom)%adjlist(j)
         if (.not. mol2%adjmat(atomperm(iatom), atomperm(jatom))) then
!            write (stderr, *) 'remove mol1 bond:', iatom, jatom
            call mol1%remove_bond(iatom, jatom)
!            indices1 = mnatypes%parts(mnatypes%itemdir1(jatom))%indices1
!            do k = 1, size(indices1)
!               katom = indices1(k)
!               call mol1%remove_bond(iatom, katom)
!               call mol2%remove_bond(atomperm(iatom), atomperm(katom))
!            end do
         end if
      end do
   end do

   do iatom = 1, size(mol2%atoms)
      do j = 1, size(mol2%atoms(iatom)%adjlist)
         jatom = mol2%atoms(iatom)%adjlist(j)
         if (.not. mol1%adjmat(invatomperm(iatom), invatomperm(jatom))) then
!            write (stderr, *) 'remove mol2 bond:', iatom, jatom
            call mol2%remove_bond(iatom, jatom)
!            indices2 = mnatypes%parts(mnatypes%itemdir2(jatom))%indices2
!            do k = 1, size(indices2)
!               katom = indices2(k)
!               call mol2%remove_bond(iatom, katom)
!               call mol1%remove_bond(invatomperm(iatom), invatomperm(katom))
!            end do
         end if
      end do
   end do

   ! Dissociate water molecules
!
!   do iatom = 1, size(molfrags1)
!      if (all(sorted(mol1%atoms(molfrags1(iatom)%n)%elnum) == [1, 1, 8])) then
!         do j = 1, size(molfrags1(iatom)%n)
!            jatom = molfrags1(iatom)%n(j)
!            do k = 1, size(mol1%atoms(jatom)%adjlist)
!               katom = mol1%atoms(jatom)%adjlist(k)
!               call mol1%remove_bond(jatom, katom)
!            end do
!         end do
!      end if
!   end do
!
!   do iatom = 1, size(molfrags2)
!      if (all(sorted(mol2%atoms(molfrags2(iatom)%n)%elnum) == [1, 1, 8])) then
!         do j = 1, size(molfrags2(iatom)%n)
!            jatom = molfrags2(iatom)%n(j)
!            do k = 1, size(mol2%atoms(jatom)%adjlist)
!               katom = mol2%atoms(jatom)%adjlist(k)
!               call mol2%remove_bond(jatom, katom)
!            end do
!         end do
!      end if
!   end do

end subroutine

end module
