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

module alignment_reaction
use parameters
use options
use random
use molecule
use utils
use chemistry
use permutation
use euclidean
use assignment_cluster
use adjacency
use biasing
use pruning
use lcrs_trees
use lcrs_arrays
use hna
use assignment_tree
use assignment_conformer
use registration
use alignment_conformer

implicit none

contains

subroutine find_reactive_bonds( atoms1, atoms2, atomtypes, atomperm)
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_t), intent(in) :: atomtypes
   integer, dimension(:), intent(out) :: atomperm

   ! Local variables
   type(partition_t) :: scnatypes
   type(assigntree_node_t), pointer :: hnachain
   type(int_list), dimension(:), allocatable :: molfrags1, molfrags2
   type(int_matrix), dimension(:), allocatable :: biases
   type(real_matrix), dimension(:), allocatable :: minbiases
   type(registry_t), target :: registry
   integer :: adjd
   integer, dimension(:), allocatable :: elnums1, elnums2
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   real(rk), dimension(:), allocatable :: weights1, weights2
   real(rk), dimension(:,:), allocatable :: wcoords1, wcoords2
   real(rk) :: center1(3), center2(3)

   elnums1 = atoms1%elnum
   elnums2 = atoms2%elnum
   weights1 = atomic_weights(elnums1)/sum(atomic_weights(elnums1))
   weights2 = atomic_weights(elnums2)/sum(atomic_weights(elnums2))
   wcoords1 = get_coords( atoms1)
   wcoords2 = get_coords( atoms2)
   adjmat1 = get_adjmat( atoms1)
   adjmat2 = get_adjmat( atoms2)

   ! Mirror coordinates
   if (mirror_flag) then
      call mirror_coords( wcoords2)
   end if

   ! Mass weight coordinates
   call weight_coords( wcoords1, weights1)
   call weight_coords( wcoords2, weights2)

   ! Calculate centroids
   center1 = centroid( wcoords1, weights1)
   center2 = centroid( wcoords2, weights2)

   ! Translate atoms to their centroids
   call translate_coords( wcoords2, -center2)
   call translate_coords( wcoords2, center1)

   ! Compute HNA types
   call compute_hna_partition( atoms1, atoms2, atomtypes, hnachain)
   call link_to_partition( hnachain%last_link, scnatypes)

   ! Find molecular fragments
   call find_molfrags( atoms1, first_partition(atomtypes), molfrags1)
   call find_molfrags( atoms2, second_partition(atomtypes), molfrags2)

   ! Find unfeasible assignments
   call compute_hna_biases( atoms1, atoms2, atomtypes, biases)

   call build_minbiases(atomtypes, atoms1, atoms2, biases, minbiases)
!CZGC: calcular min_rmsd_matrix usando 'minimum_rmsd' con vecinos en spatial.f90
!CZGC: sumar min_rmsd_matrix a biases

   ! Initialize random number generator
   call random_initialize()

   ! Initialize local minima registry
   call init_adjd_registry( registry, max_records)

   ! Optimize atom permutation
   do while (registry%records(1)%count < count_thres .and. registry%num_trials < max_trials)

      ! Assign atoms with current orientation
      call assign_atoms_biased( atomtypes, wcoords1, wcoords2, biases, atomperm)
      ! Reassign mismatches
      call minadjdiff( atomtypes, scnatypes, molfrags1, atoms1, atoms2, wcoords1, wcoords2, atomperm)
      ! Update results
      adjd = adjacencydiff( atomperm, adjmat1, adjmat2)
      call push_record( registry, atomperm, 1, adjd=adjd)

   end do

   atomperm = registry%records(1)%atomperm
end subroutine

subroutine remove_reactive_bonds( atoms1, atoms2, atomtypes, atomperm)
   ! Remove reactive bonds
   type(atom_t), dimension(:), intent(inout) :: atoms1, atoms2
   type(partition_t), intent(in) :: atomtypes
   integer, dimension(:), intent(in) :: atomperm
   ! Local variables
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   integer, dimension(:), allocatable :: invperm
   integer :: j, iatom, jatom
!   integer, dimension(:), allocatable :: items1, items2
!   integer :: k, katom

   adjmat1 = get_adjmat(atoms1)
   adjmat2 = get_adjmat(atoms2)
   invperm = inverse_permutation(atomperm)

   ! Remove mismatched bonds
   do iatom = 1, size(atoms1)
      do j = 1, size(atoms1(iatom)%adjlist)
         jatom = atoms1(iatom)%adjlist(j)
         if (.not. adjmat2(atomperm(iatom), atomperm(jatom))) then
!            write (stderr, *) 'remove atoms1 bond:', iatom, jatom
            call remove_bond(atoms1, iatom, jatom)
!            items1 = scnatypes%parts(scnatypes%itemdir1(jatom))%items1
!            do k = 1, size(items1)
!               katom = items1(k)
!               call remove_bond(atoms1, iatom, katom)
!               call remove_bond(atoms2, atomperm(iatom), atomperm(katom))
!            end do
         end if
      end do
   end do

   do iatom = 1, size(atoms2)
      do j = 1, size(atoms2(iatom)%adjlist)
         jatom = atoms2(iatom)%adjlist(j)
         if (.not. adjmat1(invperm(iatom), invperm(jatom))) then
!            write (stderr, *) 'remove atoms2 bond:', iatom, jatom
            call remove_bond(atoms2, iatom, jatom)
!            items2 = scnatypes%parts(scnatypes%itemdir2(jatom))%items2
!            do k = 1, size(items2)
!               katom = items2(k)
!               call remove_bond(atoms1, invperm(iatom), invperm(katom))
!               call remove_bond(atoms2, iatom, katom)
!            end do
         end if
      end do
   end do

   ! Dissociate water molecules
!
!   do iatom = 1, size(molfrags1)
!      if (all(sorted(atoms1(molfrags1(iatom)%e)%elnum) == [1, 1, 8])) then
!         do j = 1, size(molfrags1(iatom)%e)
!            jatom = molfrags1(iatom)%e(j)
!            do k = 1, size(atoms1(jatom)%adjlist)
!               katom = atoms1(jatom)%adjlist(k)
!               call remove_bond(atoms1, jatom, katom)
!            end do
!         end do
!      end if
!   end do
!
!   do iatom = 1, size(molfrags2)
!      if (all(sorted(atoms2(molfrags2(iatom)%e)%elnum) == [1, 1, 8])) then
!         do j = 1, size(molfrags2(iatom)%e)
!            jatom = molfrags2(iatom)%e(j)
!            do k = 1, size(atoms2(jatom)%adjlist)
!               katom = atoms2(jatom)%adjlist(k)
!               call remove_bond(atoms2, jatom, katom)
!            end do
!         end do
!      end if
!   end do
end subroutine

subroutine optimize_atomperm_isomer( atoms1, atoms2, registry)
   type(atom_t), dimension(:), intent(inout) :: atoms1, atoms2
   type(registry_t), target, intent(out) :: registry

   ! Local variables
   type(partition_t) :: atomtypes
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   integer, dimension(:), allocatable :: elnums1, elnums2
   real(rk), dimension(:,:), allocatable :: wcoords1, wcoords2
   integer, dimension(:), allocatable :: atomperm

   ! Abort if molecules have different number of atoms
   if (size(atoms1) /= size(atoms2)) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if molecules are not isomers
   if (any(sorted(atoms1%elnum) /= sorted(atoms2%elnum))) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Compute atomic types
   call collect_atomtypes( atoms1, atoms2, atomtypes)

   ! Abort if there are conflicting atomic types
   if (any(atomtypes%parts%num_items1 /= atomtypes%parts%num_items2)) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      stop
   end if

   allocate (atomperm(size(atoms1)))

   elnums1 = atoms1%elnum
   elnums2 = atoms2%elnum
   wcoords1 = get_coords(atoms1)
   wcoords2 = get_coords(atoms2)
   adjmat1 = get_adjmat(atoms1)
   adjmat2 = get_adjmat(atoms2)

   call find_reactive_bonds( atoms1, atoms2, atomtypes, atomperm)
   write (stderr, *) 'before', adjacencydiff( atomperm, adjmat1, adjmat2)
   call remove_reactive_bonds( atoms1, atoms2, atomtypes, atomperm)
   adjmat1 = get_adjmat( atoms1)
   adjmat2 = get_adjmat( atoms2)
   write (stderr, *) 'after ', adjacencydiff( atomperm, adjmat1, adjmat2)
   call optimize_atomperm_conform( atoms1, atoms2, atomtypes, registry)
end subroutine

end module
