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
use pruning
use printing
!use backtracking

implicit none

contains

subroutine remap_atoms(mol0, mol1, maplist, countlist, nrec)
   type(molecule_type), intent(in) :: mol0, mol1
   integer, intent(out) :: maplist(:, :)
   integer, intent(out) :: countlist(:)
   integer, intent(out) :: nrec

   ! Local variables

   integer :: natom
   real(rk), dimension(:), allocatable :: weights
   real(rk), dimension(:, :), allocatable :: coords0, coords1
   logical, dimension(:, :), allocatable :: adjmat0, adjmat1
   type(atomlist_type), allocatable, dimension(:) :: adjlists0, adjlists1
   type(atompartition_type) :: eltypes0, eltypes1

   logical :: visited, overflow
   integer :: irec, krec, ntrial, nstep, steps
   integer :: adjd, recadjd(maxrec)
   integer, dimension(mol0%natom) :: mapping, newmapping
   real(rk) :: rmsd, totalrot
   real(rk), dimension(4) :: rotquat, prodquat
   real(rk), dimension(maxrec) :: recrmsd, avgsteps, avgtotalrot, avgrealrot
   logical :: prunemat(mol0%natom, mol1%natom)
   real(rk) :: biasmat(mol0%natom, mol1%natom)
   real(rk) :: workcoords1(3, mol1%natom)

   integer h

   natom = mol0%get_natom()
   weights = mol0%get_weights()
   coords0 = mol0%get_coords()
   coords1 = mol1%get_coords()
   eltypes0 = mol0%get_eltypes()
   eltypes1 = mol1%get_eltypes()
   adjlists0 = mol0%get_adjlists()
   adjlists1 = mol1%get_adjlists()
   adjmat0 = mol0%get_adjmatrix()
   adjmat1 = mol1%get_adjmatrix()

   ! Calculate prune matrix

   call prune_procedure(eltypes0, eltypes1, coords0, coords1, prunemat)

   ! Calculate bias matrix

   call bias_procedure(eltypes0, eltypes1, adjlists0, adjlists1, biasmat)

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

      call assign_atoms(eltypes0, eltypes1, coords0, workcoords1, prunemat, biasmat, mapping)
      rotquat = leastrotquat(natom, weights, coords0, workcoords1, mapping)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate(natom, workcoords1, rotquat)
!      print *, sqrt(leastsquaredist(natom, weights, coords0, coords1, mapping)/sum(weights))
!      stop

      steps = 1

      do while (iter_flag)
         call assign_atoms(eltypes0, eltypes1, coords0, workcoords1, prunemat, biasmat, newmapping)
         if (all(newmapping == mapping)) exit
         rotquat = leastrotquat(natom, weights, coords0, workcoords1, newmapping)
         prodquat = quatmul(rotquat, prodquat)
         call rotate(natom, workcoords1, rotquat)
         totalrot = totalrot + rotangle(rotquat)
         mapping = newmapping
         steps = steps + 1
      end do

      nstep = nstep + steps

!      if (back_flag) then
!         call minadjdiff(mol0, mol1, mapping)
!         call eqvatomperm(mol0, mol1, workcoords1, mapping)
!      end if

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
         krec = nrec + 1
         do irec = nrec, 1, -1
            if (adjd < recadjd(irec) .or. (adjd == recadjd(irec) .and. rmsd < recrmsd(irec))) then
               krec = irec
            else
               exit
            end if
         end do
         if (nrec < maxrec) then
            nrec = nrec + 1
         else
            overflow = .true.
         end if
         if (krec <= maxrec) then
            do irec = nrec, krec + 1, -1
               countlist(irec) = countlist(irec - 1)
               recrmsd(irec) = recrmsd(irec - 1)
               recadjd(irec) = recadjd(irec - 1)
               avgsteps(irec) = avgsteps(irec - 1)
               avgrealrot(irec) = avgrealrot(irec - 1)
               avgtotalrot(irec) = avgtotalrot(irec - 1)
               maplist(:, irec) = maplist(:, irec - 1)
            end do
            countlist(krec) = 1
            recrmsd(krec) = rmsd
            recadjd(krec) = adjd
            avgsteps(krec) = steps
            avgrealrot(krec) = rotangle(prodquat)
            avgtotalrot(krec) = totalrot
            maplist(:, krec) = mapping
         end if
      end if

   end do

   ! Reorder back to default atom ordering

!   do irec = 1, nrec
!      maplist(:, irec) = mol1%atomorder(maplist(mol0%atomordermap, irec))
!   end do

   ! Print stats if requested

   if (stats_flag) then
      call print_stats(nrec, countlist, avgsteps, avgtotalrot, avgrealrot, recadjd, recrmsd)
      call print_final_stats(overflow, maxrec, nrec, ntrial, nstep)
   end if

end subroutine

end module
