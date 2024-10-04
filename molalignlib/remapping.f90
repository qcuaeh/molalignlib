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

module remapping
use kinds
use molecule
use flags
use bounds
use strutils
use permutation
use translation
use rotation
use alignment
use assignment
use adjacency
use assorting
use biasing
use backtracking
use printing

implicit none

private
public optimize_mapping
public remove_reactive_bonds

contains

subroutine optimize_mapping(mol0, mol1, maplist, countlist, nrec)
   type(t_mol), intent(in) :: mol0, mol1
   integer, intent(out) :: maplist(:, :)
   integer, intent(out) :: countlist(:)
   integer, intent(out) :: nrec

   integer :: i
   integer :: natom
   integer :: neltype
   integer :: nmnatype0, nmnatype1
   real(rk), dimension(:), allocatable :: weights
   real(rk), dimension(:, :), allocatable :: coords0, coords1
   logical, dimension(:, :), allocatable :: adjmat0, adjmat1
   logical :: visited, overflow
   integer :: nfrag0, nfrag1
   integer :: irec, mrec, ntrial, nstep, steps
   integer :: adjd, recadjd(maxrec)
   integer :: nadjs0(mol0%natom), adjlists0(maxcoord, mol0%natom)
   integer :: nadjs1(mol1%natom), adjlists1(maxcoord, mol1%natom)
   integer, dimension(mol0%natom) :: eltypepartlens
   integer, dimension(mol0%natom, maxlevel) :: nadjmna0
   integer, dimension(mol1%natom, maxlevel) :: nadjmna1
   integer, dimension(maxcoord, mol0%natom, maxlevel) :: adjmnalen0
   integer, dimension(maxcoord, mol1%natom, maxlevel) :: adjmnalen1
   integer, dimension(maxcoord, mol0%natom, maxlevel) :: adjmnalist0
   integer, dimension(maxcoord, mol1%natom, maxlevel) :: adjmnalist1
   integer, dimension(mol0%natom) :: nadjmnatypes0
   integer, dimension(mol1%natom) :: nadjmnatypes1
   integer, dimension(maxcoord, mol0%natom) :: adjmnatypepartlens0
   integer, dimension(maxcoord, mol1%natom) :: adjmnatypepartlens1
   integer :: equivmat(mol0%natom, mol1%natom)
   integer, dimension(mol0%natom) :: mapping, newmapping
   real(rk) :: rmsd, totalrot
   real(rk), dimension(4) :: rotquat, prodquat
   real(rk), dimension(maxrec) :: recrmsd, avgsteps, avgtotalrot, avgrealrot
   logical :: prunemat(mol0%natom, mol1%natom)
   real(rk) :: biasmat(mol0%natom, mol1%natom)
   real(rk) :: workcoords1(3, mol1%natom)

   integer, allocatable, dimension(:) :: fragroots0, fragroots1
   integer, allocatable, dimension(:) :: mnatypepartlens0, mnatypepartlens1

   natom = mol0%natom
   weights = mol0%get_sorted_weights()
   coords0 = mol0%get_sorted_coords()
   coords1 = mol1%get_sorted_coords()
   adjmat0 = mol0%get_sorted_adjmat()
   adjmat1 = mol1%get_sorted_adjmat()

   neltype = mol0%get_neltype()
   eltypepartlens = mol0%get_eltypepartlens()

   nadjs0 = mol0%get_sorted_nadjs()
   nadjs1 = mol1%get_sorted_nadjs()
   adjlists0 = mol0%get_sorted_oldadjlists()
   adjlists1 = mol1%get_sorted_oldadjlists()

   nmnatype0 = mol0%get_nmnatype()
   nmnatype1 = mol1%get_nmnatype()
   mnatypepartlens0 = mol0%get_mnatypepartlens()
   mnatypepartlens1 = mol1%get_mnatypepartlens()

   nadjmnatypes0 = mol0%get_sorted_nadjmnatypes()
   nadjmnatypes1 = mol1%get_sorted_nadjmnatypes()
   adjmnatypepartlens0 = mol0%get_sorted_adjmnatypepartlens()
   adjmnatypepartlens1 = mol1%get_sorted_adjmnatypepartlens()

   fragroots0 = mol0%get_sorted_fragroots()
   fragroots1 = mol1%get_sorted_fragroots()

   nfrag0 = size(fragroots0)
   nfrag1 = size(fragroots1)

   ! Calculate MNA equivalence matrix

   call fill_equivmat(mol0, mol1, nadjmna0, adjmnalen0, adjmnalist0, nadjmna1, &
         adjmnalen1, adjmnalist1, equivmat)

   ! Calculate bias matrix

   call cross_function(natom, neltype, eltypepartlens, coords0, coords1, equivmat, prunemat, biasmat)

   ! Initialize loop variables

   nrec = 0
   nstep = 0
   ntrial = 0
   countlist(1) = 0
   overflow = .false.

   ! Loop for map searching

   do while (countlist(1) < maxcount .and. (.not. trial_flag .or. ntrial < maxtrials))

      ntrial = ntrial + 1

      ! Work with a copy of coords1

      workcoords1 = coords1

      ! Aply a random rotation to workcoords1

      call rotate(natom, workcoords1, randrotquat())

      ! Minimize the euclidean distance

      call assign_atoms_function(natom, neltype, eltypepartlens, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
         nadjmna1, adjmnalen1, adjmnalist1, workcoords1, weights, equivmat, prunemat, biasmat, mapping)
      rotquat = leastrotquat(natom, weights, coords0, workcoords1, mapping)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate(natom, workcoords1, rotquat)

      steps = 1

      do while (iter_flag)
         call assign_atoms_function(natom, neltype, eltypepartlens, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
            nadjmna1, adjmnalen1, adjmnalist1, workcoords1, weights, equivmat, prunemat, biasmat, newmapping)
         if (all(newmapping == mapping)) exit
         rotquat = leastrotquat(natom, weights, coords0, workcoords1, newmapping)
         prodquat = quatmul(rotquat, prodquat)
         call rotate(natom, workcoords1, rotquat)
         totalrot = totalrot + rotangle(rotquat)
         mapping = newmapping
         steps = steps + 1
      end do

      nstep = nstep + steps

      if (back_flag) then

         call minadjdiff(natom, weights, neltype, eltypepartlens, coords0, nadjs0, adjlists0, adjmat0, nmnatype0, &
            mnatypepartlens0, workcoords1, nadjs1, adjlists1, adjmat1, nmnatype1, mnatypepartlens1, mapping, nfrag0, fragroots0)

         call eqvatomperm(natom, weights, coords0, adjmat0, adjlists0, nmnatype0, mnatypepartlens0, &
            nadjmnatypes0, adjmnatypepartlens0, workcoords1, adjmat1, mapping, nfrag0, fragroots0)

      end if

      adjd = adjacencydiff(natom, adjmat0, adjmat1, mapping)
      rmsd = sqrt(leastsquaredist(natom, weights, coords0, coords1, mapping)/sum(weights))

      ! Check for new best maplist

      visited = .false.

      do irec = 1, nrec
         if (all(mapping == maplist(:, irec))) then
            countlist(irec) = countlist(irec) + 1
            avgsteps(irec) = avgsteps(irec) + (steps - avgsteps(irec))/countlist(irec)
            avgrealrot(irec) = avgrealrot(irec) + (rotangle(prodquat) - avgrealrot(irec))/countlist(irec)
            avgtotalrot(irec) = avgtotalrot(irec) + (totalrot - avgtotalrot(irec))/countlist(irec)
            visited = .true.
            exit
         end if
      end do

      if (.not. visited) then
         mrec = nrec + 1
         do irec = nrec, 1, -1
            if (adjd < recadjd(irec) .or. (adjd == recadjd(irec) .and. rmsd < recrmsd(irec))) then
               mrec = irec
            else
               exit
            end if
         end do
         if (nrec < maxrec) then
            nrec = nrec + 1
         else
            overflow = .true.
         end if
         if (mrec <= maxrec) then
            do irec = nrec, mrec + 1, -1
               countlist(irec) = countlist(irec - 1)
               recrmsd(irec) = recrmsd(irec - 1)
               recadjd(irec) = recadjd(irec - 1)
               avgsteps(irec) = avgsteps(irec - 1)
               avgrealrot(irec) = avgrealrot(irec - 1)
               avgtotalrot(irec) = avgtotalrot(irec - 1)
               maplist(:, irec) = maplist(:, irec - 1)
            end do
            countlist(mrec) = 1
            recrmsd(mrec) = rmsd
            recadjd(mrec) = adjd
            avgsteps(mrec) = steps
            avgrealrot(mrec) = rotangle(prodquat)
            avgtotalrot(mrec) = totalrot
            maplist(:, mrec) = mapping
         end if
      end if

   end do

   ! Reorder back to default atom ordering

   do i = 1, nrec
      maplist(:, i) = mol1%mnatypepartition%atom_order(maplist(mol0%mnatypepartition%atom_mapping, i))
   end do

   ! Print stats if requested

   if (stats_flag) then
      call print_stats(nrec, countlist, avgsteps, avgtotalrot, avgrealrot, recadjd, recrmsd)
      call print_final_stats(overflow, maxrec, nrec, ntrial, nstep)
   end if

end subroutine

subroutine remove_reactive_bonds(mol0, mol1, mapping)

   type(t_mol), intent(inout) :: mol0, mol1
   integer, dimension(:), intent(in) :: mapping

   integer :: i, j, k, l, j_, k_, l_
   integer :: natom0, natom1
   integer, allocatable, dimension(:) :: elnums0, elnums1
   integer, allocatable, dimension(:) :: mnatypes0, mnatypes1
   logical, allocatable, dimension(:, :) :: adjmat0, adjmat1
   type(t_atomlist), allocatable, dimension(:) :: adjlists0, adjlists1
   type(t_atomlist), allocatable, dimension(:) :: molfragparts0, molfragparts1
   type(t_atomlist), allocatable, dimension(:) :: mnatypeparts0, mnatypeparts1
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
   adjmat0 = mol0%get_adjmat()
   adjmat1 = mol1%get_adjmat()
   elnums0 = mol0%get_elnums()
   elnums1 = mol1%get_elnums()
   adjlists0 = mol0%get_newadjlists()
   adjlists1 = mol1%get_newadjlists()
   mnatypes0 = mol0%get_mnatypes()
   mnatypes1 = mol1%get_mnatypes()
   mnatypeparts0 = mol0%get_mnatypeparts()
   mnatypeparts1 = mol1%get_mnatypeparts()
   molfragparts0 = mol0%get_molfragparts()
   molfragparts1 = mol1%get_molfragparts()
   unmapping = inverse_permutation(mapping)

   ! Remove mismatched bonds

   do i = 1, natom0
      do j_ = 1, size(adjlists0(i)%atomidcs)
         j = adjlists0(i)%atomidcs(j_)
         if (.not. adjmat1(mapping(i), mapping(j))) then
            mnatypepartidcs0 = mnatypeparts0(mnatypes0(j))%atomidcs
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
            mnatypepartidcs1 = mnatypeparts1(mnatypes1(j))%atomidcs
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
