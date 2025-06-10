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
use basetypes
use globals
use random
use spatial
use adjacency
use permutation
use assignment
use biasing
use molecule
use tracking
use lcrs_tree
use eltype_compute
use mna_compute
use registry

implicit none

contains

subroutine find_reactive_bonds( mol1, mol2, eltypes, atomperm)
   type(mol_type), intent(in) :: mol1, mol2
   type(partition_t), intent(in) :: eltypes
   integer, dimension(:), intent(out) :: atomperm

   ! Local variables
   type(partition_t) :: mnatypes
   type(part_node_t), pointer :: root_part
   type(chain_node_t), pointer :: mnachain
   type(int_list), dimension(:), allocatable :: molfrags1, molfrags2
   type(int_matrix), dimension(:), allocatable :: biases
   type(topoatomperm_registry), target :: results
   integer, dimension(:), allocatable :: auxperm
   integer, dimension(:), allocatable :: elnums1, elnums2
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   real(rk), dimension(:,:), allocatable :: coords1, coords2
   real(rk) :: center1(3), center2(3)

   allocate (auxperm, mold=atomperm)

   elnums1 = mol1%atoms%elnum
   elnums2 = mol2%atoms%elnum
   coords1 = get_coords( mol1)
   coords2 = get_coords( mol2)
   adjmat1 = get_adjmat( mol1)
   adjmat2 = get_adjmat( mol2)

   ! Mirror coordinates
   if (mirror_flag) then
      call mirror_coords( coords2)
   end if

   ! Mass weight coordinates
   call weight_coords( coords1, atomic_weights(elnums1))
   call weight_coords( coords2, atomic_weights(elnums2))

   ! Translate atoms to their centroids
   center1 = centroid( coords1)
   center2 = centroid( coords2)
   call translate_coords( coords2, -center2)
   call translate_coords( coords2, center1)

   ! Compute MNA types
   call init_chain_from_partition( eltypes, mnachain, root_part)
   call compute_consistent_mnas( mol1, mol2, mnachain)
   call link_to_partition( mnachain%last_link, mnatypes)

   ! Find molecular fragments
   call find_molfrags( mol1, first_partition(eltypes), molfrags1)
   call find_molfrags( mol2, second_partition(eltypes), molfrags2)

   ! Find unfeasible assignments
   call compute_mna_biases( mol1, mol2, eltypes, biases)

   ! Initialize random number generator
   call random_initialize()

   ! Initialize local minima registry
   call registry_init( results, max_records, adjmat1)

   ! Optimize atom permutation
   do while (results%records(1)%count < max_count .and. results%num_trials < max_trials)

      ! Assign atoms with current orientation
      call assign_atoms_biased( eltypes, coords1, coords2, biases, atomperm)
      ! Reassign mismatches
      call minadjdiff( eltypes, mnatypes, molfrags1, mol1, mol2, coords1, coords2, atomperm)
      ! Update results
      call registry_push( results, adjmat2, atomperm)

   end do

   atomperm = results%records(1)%atomperm
end subroutine

subroutine remove_reactive_bonds( mol1, mol2, eltypes, atomperm)
   ! Remove reactive bonds
   type(mol_type), intent(inout) :: mol1, mol2
   type(partition_t), intent(in) :: eltypes
   integer, dimension(:), intent(in) :: atomperm
   ! Local variables
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   integer, dimension(:), allocatable :: invatomperm
   integer :: j, iatom, jatom
!   integer, dimension(:), allocatable :: items1, items2
!   integer :: k, katom

   adjmat1 = get_adjmat(mol1)
   adjmat2 = get_adjmat(mol2)
   invatomperm = inverse_perm(atomperm)

   ! Remove mismatched bonds
   do iatom = 1, size(mol1%atoms)
      do j = 1, size(mol1%atoms(iatom)%adjlist)
         jatom = mol1%atoms(iatom)%adjlist(j)
         if (.not. adjmat2(atomperm(iatom), atomperm(jatom))) then
!            write (stderr, *) 'remove mol1 bond:', iatom, jatom
            call remove_bond(mol1, iatom, jatom)
!            items1 = mnatypes%parts(mnatypes%itemdir1(jatom))%items1
!            do k = 1, size(items1)
!               katom = items1(k)
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
!            items2 = mnatypes%parts(mnatypes%itemdir2(jatom))%items2
!            do k = 1, size(items2)
!               katom = items2(k)
!               call remove_bond(mol1, invatomperm(iatom), invatomperm(katom))
!               call remove_bond(mol2, iatom, katom)
!            end do
         end if
      end do
   end do

   ! Dissociate water molecules
!
!   do iatom = 1, size(molfrags1)
!      if (all(sorted(mol1%atoms(molfrags1(iatom)%e)%elnum) == [1, 1, 8])) then
!         do j = 1, size(molfrags1(iatom)%e)
!            jatom = molfrags1(iatom)%e(j)
!            do k = 1, size(mol1%atoms(jatom)%adjlist)
!               katom = mol1%atoms(jatom)%adjlist(k)
!               call remove_bond(mol1, jatom, katom)
!            end do
!         end do
!      end if
!   end do
!
!   do iatom = 1, size(molfrags2)
!      if (all(sorted(mol2%atoms(molfrags2(iatom)%e)%elnum) == [1, 1, 8])) then
!         do j = 1, size(molfrags2(iatom)%e)
!            jatom = molfrags2(iatom)%e(j)
!            do k = 1, size(mol2%atoms(jatom)%adjlist)
!               katom = mol2%atoms(jatom)%adjlist(k)
!               call remove_bond(mol2, jatom, katom)
!            end do
!         end do
!      end if
!   end do
end subroutine

end module
