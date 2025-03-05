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
use spatial
use adjacency
use permutation
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
   type(adjd_registry), target :: results
   type(tree_node), pointer :: mnatree
   type(bipartition_container) :: mnatypes
   type(intlist_type), dimension(:), allocatable :: molfrags1, molfrags2
   type(intmatrix_type), allocatable :: biases(:)
   integer, dimension(:), allocatable :: auxperm, invatomperm
   integer, dimension(:), allocatable :: elnums1, elnums2
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   real(rk), dimension(:,:), allocatable :: coords1, coords2
   real(rk) :: step_rotation(4), total_rotation(4)
   integer :: num_steps
   integer :: j, iatom, jatom
!   integer, dimension(:), allocatable :: indices1, indices2
!   integer :: k, katom

   allocate (auxperm, mold=atomperm)

   elnums1 = mol1%atoms%elnum
   elnums2 = mol2%atoms%elnum
   coords1 = get_coords(mol1)
   coords2 = get_coords(mol2)
   adjmat1 = get_adjmat(mol1)
   adjmat2 = get_adjmat(mol2)

   ! Mirror coordinates
   if (mirror_flag) then
      call mirror_coords( coords2)
   end if

   ! Mass weight coordinates
   call weight_coords(coords1, atomic_weights(elnums1))
   call weight_coords(coords2, atomic_weights(elnums2))

   ! Translate atoms to their centroids
   call translate_coords( coords1, -centroid(coords1))
   call translate_coords( coords2, -centroid(coords2))

   ! Compute MNA types
   call tree_from_partition( eltypes, mnatree)
   call compute_consistent_mnatypes( mol1, mol2, mnatree)
   call partition_from_tree( mnatree, mnatypes)

   ! Find molecular fragments
   call find_molfrags( mol1, first_partition(eltypes), molfrags1)
   call find_molfrags( mol2, second_partition(eltypes), molfrags2)

   ! Find unfeasible assignments
   call compute_mna_biases( mol1, mol2, eltypes, biases)

   ! Initialize random number generator
   call random_initialize()

   ! Initialize local minima registry
   call registry_init(results, max_records, coords1, adjmat1)

   ! Optimize atom permutation
   do while (results%records(1)%count < max_count .and. results%num_trials < max_trials)

      ! Aply a random rotation to coords2
      call rotate_coords(coords2, randrotquat())

      ! Assign atoms with current orientation
      call assign_atoms_biased(eltypes, coords1, coords2, biases, atomperm)
      total_rotation = optimal_rotation(atomperm, coords1, coords2)
      call rotate_coords(coords2, total_rotation)
      num_steps = 1

      do while (iter_flag)
         call assign_atoms_biased(eltypes, coords1, coords2, biases, auxperm)
         if (all(auxperm == atomperm)) exit
         atomperm = auxperm
         step_rotation = optimal_rotation(atomperm, coords1, coords2)
         call rotate_coords(coords2, step_rotation)
         total_rotation = quatmul(step_rotation, total_rotation)
         num_steps = num_steps + 1
      end do

      call minadjdiff( eltypes, mnatypes, molfrags1, mol1, mol2, coords1, coords2, atomperm)

      ! Update results
      call registry_push( results, coords2, adjmat2, atomperm, num_steps, total_rotation)

   end do

   ! Remove reactive bonds

   atomperm = results%records(1)%atomperm
   invatomperm = inverse_perm(atomperm)

   ! Remove mismatched bonds
   do iatom = 1, size(mol1%atoms)
      do j = 1, size(mol1%atoms(iatom)%adjlist)
         jatom = mol1%atoms(iatom)%adjlist(j)
         if (.not. adjmat2(atomperm(iatom), atomperm(jatom))) then
!            write (stderr, *) 'remove mol1 bond:', iatom, jatom
            call remove_bond(mol1, iatom, jatom)
!            indices1 = mnatypes%parts(mnatypes%itemdir1(jatom))%indices1
!            do k = 1, size(indices1)
!               katom = indices1(k)
!               call remove_bond(mol1, iatom, katom)
!               call remove_bond(mol2, atomperm(iatom), atomperm(katom))
!            end do
         end if
      end do
   end do

   do iatom = 1, size(mol2%atoms)
      do j = 1, size(mol2%atoms(iatom)%adjlist)
         jatom = mol2%atoms(iatom)%adjlist(j)
         if (.not. adjmat1(invatomperm(iatom), invatomperm(jatom))) then
!            write (stderr, *) 'remove mol2 bond:', iatom, jatom
            call remove_bond(mol2, iatom, jatom)
!            indices2 = mnatypes%parts(mnatypes%itemdir2(jatom))%indices2
!            do k = 1, size(indices2)
!               katom = indices2(k)
!               call remove_bond(mol1, invatomperm(iatom), invatomperm(katom))
!               call remove_bond(mol2, iatom, katom)
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
!               call remove_bond(mol1, jatom, katom)
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
!               call remove_bond(mol2, jatom, katom)
!            end do
!         end do
!      end if
!   end do

end subroutine

end module
