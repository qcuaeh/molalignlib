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
use types
use flags
use bounds
use random
use strutils
use adjacency
use rotation
use translation
use tracking
use backtracking
use alignment
use assorting
use printing
use biasing
use assignment

implicit none

private
public optimize_mapping
public remove_reactive_bonds

contains

subroutine optimize_mapping(mol0, mol1, mnaord0, mnaord1, maplist, countlist, nrec)
   type(cMol), intent(in) :: mol0, mol1
   integer, intent(in), dimension(:) :: mnaord0, mnaord1
   integer, intent(out) :: maplist(:, :)
   integer, intent(out) :: countlist(:)
   integer, intent(out) :: nrec

   integer :: i
   integer :: natom
   integer :: neltype
   integer :: nmnatype0, nmnatype1
   real(wp), dimension(:), allocatable :: weights
   real(wp), dimension(:, :), allocatable :: coords0, coords1
   logical, dimension(:, :), allocatable :: adjmat0, adjmat1
   logical :: visited, overflow
   integer :: nfrag0, nfrag1
   integer irec, mrec, ntrial, nstep, steps
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
   integer, dimension(mol0%natom) :: natomneimnatypes0
   integer, dimension(mol1%natom) :: natomneimnatypes1
   integer, dimension(maxcoord, mol0%natom) :: atomneimnatypepartlens0
   integer, dimension(maxcoord, mol1%natom) :: atomneimnatypepartlens1
   integer :: equivmat(mol0%natom, mol1%natom)
   integer, dimension(mol0%natom) :: mapping, newmapping
   real(wp) :: rmsd, totalrot
   real(wp), dimension(4) :: rotquat, prodquat
   real(wp), dimension(maxrec) :: recrmsd, avgsteps, avgtotalrot, avgrealrot
   real(wp) :: biasmat(mol0%natom, mol1%natom)
   real(wp) :: workcoords1(3, mol1%natom)

   integer, allocatable, dimension(:) :: fragroots0, fragroots1
   integer, allocatable, dimension(:) :: mnatypepartlens0, mnatypepartlens1
   integer, allocatable, dimension(:) :: unord0, unord1

   natom = mol0%natom
   weights = mol0%get_atomweights(mnaord0)
   coords0 = mol0%get_coords(mnaord0)
   coords1 = mol1%get_coords(mnaord1)
   adjmat0 = mol0%get_adjmat(mnaord0)
   adjmat1 = mol1%get_adjmat(mnaord1)

   neltype = mol0%get_neltype()
   eltypepartlens = mol0%get_eltypepartlens()

   nadjs0 = mol0%get_nadjs(mnaord0)
   nadjs1 = mol1%get_nadjs(mnaord1)
   adjlists0 = mol0%get_adjlists(mnaord0)
   adjlists1 = mol1%get_adjlists(mnaord1)

   nmnatype0 = mol0%get_nmnatype()
   nmnatype1 = mol1%get_nmnatype()
   mnatypepartlens0 = mol0%get_mnatypepartlens()
   mnatypepartlens1 = mol1%get_mnatypepartlens()

   natomneimnatypes0 = mol0%get_natomneimnatypes(mnaord0)
   natomneimnatypes1 = mol1%get_natomneimnatypes(mnaord1)
   atomneimnatypepartlens0 = mol0%get_atomneimnatypepartlens(mnaord0)
   atomneimnatypepartlens1 = mol1%get_atomneimnatypepartlens(mnaord1)

   fragroots0 = mol0%get_fragroots(mnaord0)
   fragroots1 = mol1%get_fragroots(mnaord1)

   nfrag0 = size(fragroots0)
   nfrag1 = size(fragroots1)

   ! Calculate MNA equivalence matrix

   call calcequivmat(mol0, mol1, mnaord0, mnaord1, nadjmna0, adjmnalen0, adjmnalist0, nadjmna1, &
         adjmnalen1, adjmnalist1, equivmat)

   ! Calculate bias matrix

   call setcrossbias(natom, neltype, eltypepartlens, coords0, coords1, equivmat, biasmat)

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

      call rotate(natom, workcoords1, genrotquat(randarray(3)))

      ! Minimize the euclidean distance

      call mapatoms(natom, neltype, eltypepartlens, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
         nadjmna1, adjmnalen1, adjmnalist1, workcoords1, weights, equivmat, biasmat, mapping)
      rotquat = leastrotquat(natom, weights, coords0, workcoords1, mapping)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate(natom, workcoords1, rotquat)

      steps = 1

      do while (iter_flag)
         call mapatoms(natom, neltype, eltypepartlens, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
            nadjmna1, adjmnalen1, adjmnalist1, workcoords1, weights, equivmat, biasmat, newmapping)
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
            natomneimnatypes0, atomneimnatypepartlens0, workcoords1, adjmat1, mapping, nfrag0, fragroots0)

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

   unord0 = inverse_mapping(mnaord0)
   unord1 = inverse_mapping(mnaord1)

   do i = 1, nrec
      maplist(:, i) = mnaord1(maplist(unord0(:), i))
   end do

   ! Print stats if requested

   if (stats_flag) then
      call print_stats(nrec, countlist, avgsteps, avgtotalrot, avgrealrot, recadjd, recrmsd)
      call print_final_stats(overflow, maxrec, nrec, ntrial, nstep)
   end if

end subroutine

subroutine remove_reactive_bonds(mol0, mol1, mapping)

   type(cMol), intent(inout) :: mol0, mol1
   integer, dimension(:), intent(in) :: mapping

   integer :: natom
   integer :: i, j, k
   type(cAtom), allocatable, dimension(:) :: atoms0, atoms1
   logical, allocatable, dimension(:, :) :: adjmat0, adjmat1
   integer, allocatable, dimension(:) :: atomeltypes0, atomeltypes1
   type(cPart), allocatable, dimension(:) :: eltypeparts0, eltypeparts1
   integer, allocatable, dimension(:) :: eltypepartidcs0j, eltypepartidcs1j
   real(wp) :: atomdist, rotquat(4)

   ! Align coordinates

   rotquat = leastrotquat(mol0%natom, mol0%get_atomweights(), mol0%get_coords(), mol1%get_coords(), mapping)
   call mol1%rotate_coords(rotquat)

!   write (stderr, *) sqrt(squaredist(mol0%natom, mol0%get_atomweights(), mol0%get_coords(), mol1%get_coords(), mapping) &
!         /sum(mol0%get_atomweights()))

   ! Initialization

   atoms0 = mol0%get_atoms()
   atoms1 = mol1%get_atoms()
   adjmat0 = mol0%get_adjmat()
   adjmat1 = mol1%get_adjmat()
   atomeltypes0 = mol0%get_atomeltypes()
   atomeltypes1 = mol1%get_atomeltypes()
   eltypeparts0 = mol0%get_eltypeparts()
   eltypeparts1 = mol1%get_eltypeparts()

   ! Remove reactive bonds

   do i = 1, mol0%natom
      do j = 1, mol0%natom
         if (adjmat0(i, j) .neqv. adjmat1(mapping(i), mapping(j))) then
            eltypepartidcs0j = eltypeparts0(atomeltypes0(j))%atomidcs
            eltypepartidcs1j = eltypeparts1(atomeltypes1(j))%atomidcs
!            write (stderr, '(i3,2x,i3)') i, j
            call remove_reactive_bond(i, j, mol0, mol1, mapping)
            do k = 1, size(eltypepartidcs1j)
               atomdist = sum((atoms0(j)%coords - atoms1(mapping(eltypepartidcs1j(k)))%coords)**2)
!               write (stderr, '(a3,2x,i3,2x,i3,2x,f8.4)') '<<<', i, eltypepartidcs1j(k), atomdist
               if (atomdist < 2.0) then
                  call remove_reactive_bond(i, eltypepartidcs1j(k), mol0, mol1, mapping)
               end if
            end do
            do k = 1, size(eltypepartidcs0j)
               atomdist = sum((atoms0(eltypepartidcs0j(k))%coords - atoms1(mapping(j))%coords)**2)
!               write (stderr, '(a3,2x,i3,2x,i3,2x,f8.4)') '>>>', i, eltypepartidcs0j(k), atomdist
               if (atomdist < 2.0) then
                  call remove_reactive_bond(i, eltypepartidcs0j(k), mol0, mol1, mapping)
               end if
            end do
         end if
      end do
   end do

end subroutine

subroutine remove_reactive_bond(i, j, mol0, mol1, mapping)
! Purpose: Remove reactive bonds
   integer, intent(in) :: i, j
   integer, dimension(:), intent(in) :: mapping
   type(cMol), intent(inout) :: mol0, mol1

   integer :: k
   integer, allocatable, dimension(:) :: atomnums0, atomnums1
   integer, allocatable, dimension(:) :: nadjs0, nadjs1
   integer, allocatable, dimension(:, :) :: adjlists0, adjlists1

   atomnums0 = mol0%get_elnums()
   atomnums1 = mol1%get_elnums()
   nadjs0 = mol0%get_nadjs()
   nadjs1 = mol1%get_nadjs()
   adjlists0 = mol0%get_adjlists()
   adjlists1 = mol1%get_adjlists()

   call mol0%remove_bond(i, j)
   call mol1%remove_bond(mapping(i), mapping(j))

   if (atomnums0(i) == 1) then
      do k = 1, nadjs0(i)
         if (any([7, 8] == atomnums0(adjlists0(k, i)))) then
            call mol0%remove_bond(i, adjlists0(k, i))
         end if
      end do
   end if

   if (atomnums1(i) == 1) then
      do k = 1, nadjs1(i)
         if (any([7, 8] == atomnums1(adjlists1(k, i)))) then
            call mol1%remove_bond(mapping(i), mapping(adjlists1(k, i)))
         end if
      end do
   end if

   if (atomnums0(j) == 1) then
      do k = 1, nadjs0(j)
         if (any([7, 8] == atomnums0(adjlists0(k, j)))) then
            call mol0%remove_bond(j, adjlists0(k, j))
         end if
      end do
   end if

   if (atomnums1(j) == 1) then
      do k = 1, nadjs1(j)
         if (any([7, 8] == atomnums1(adjlists1(k, j)))) then
            call mol1%remove_bond(mapping(j), mapping(adjlists1(k, j)))
         end if
      end do
   end if

end subroutine

end module
