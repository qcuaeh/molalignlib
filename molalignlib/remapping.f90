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
public find_reactive_sites

contains

subroutine optimize_mapping(mol0, mol1, maplist, countlist, nrec)
!subroutine optimize_mapping( &
!   mol0, &
!   mol1, &
!   natomtype, &
!   atomtypelenlist, &
!   natomequiv0, &
!   atomequivlenlist0, &
!   natomequiv1, &
!   atomequivlenlist1, &
!   maplist, &
!   countlist, &
!   nrec)

   type(Molecule), intent(inout) :: mol0, mol1
   integer :: natomtype
   integer, dimension(mol0%natom) :: atomtypelenlist
   integer :: natomequiv0, natomequiv1
   integer, dimension(mol0%natom) :: atomequivlenlist0
   integer, dimension(mol1%natom) :: atomequivlenlist1
!   integer, intent(in) :: natomtype
!   integer, dimension(:), intent(in) :: atomtypelenlist
!   integer, intent(in) :: natomequiv0, natomequiv1
!   integer, dimension(:), intent(in) :: atomequivlenlist0, atomequivlenlist1
   integer, intent(out) :: maplist(:, :)
   integer, intent(out) :: countlist(:)
   integer, intent(out) :: nrec

   integer :: natom
   real(wp), dimension(:), allocatable :: weights
   real(wp), dimension(:, :), allocatable :: coords0, coords1
   logical, dimension(:, :), allocatable :: adjmat0, adjmat1
   logical :: visited, overflow
   integer :: nfrag0, nfrag1
   integer irec, mrec, ntrial, nstep, steps
   integer :: adjd, recadjd(maxrec)
   integer :: nadjs0(mol0%natom), adjlists0(maxcoord, mol0%natom)
   integer :: nadjs1(mol1%natom), adjlists1(maxcoord, mol1%natom)
   integer, dimension(mol0%natom, maxlevel) :: nadjmna0
   integer, dimension(mol1%natom, maxlevel) :: nadjmna1
   integer, dimension(maxcoord, mol0%natom, maxlevel) :: adjmnalen0
   integer, dimension(maxcoord, mol1%natom, maxlevel) :: adjmnalen1
   integer, dimension(maxcoord, mol0%natom, maxlevel) :: adjmnalist0
   integer, dimension(maxcoord, mol1%natom, maxlevel) :: adjmnalist1
   integer, dimension(mol0%natom) :: nadjequivs0
   integer, dimension(mol1%natom) :: nadjequivs1
   integer, dimension(maxcoord, mol0%natom) :: adjequivlenlists0
   integer, dimension(maxcoord, mol1%natom) :: adjequivlenlists1
   integer, dimension(mol0%natom) :: fragroot0
   integer, dimension(mol1%natom) :: fragroot1
   integer :: equivmat(mol0%natom, mol1%natom)
   integer, dimension(mol0%natom) :: mapping, newmapping
   real(wp) :: rmsd, totalrot
   real(wp), dimension(4) :: rotquat, prodquat
   real(wp), dimension(maxrec) :: recrmsd, avgsteps, avgtotalrot, avgrealrot
   real(wp) :: biasmat(mol0%natom, mol1%natom)
   real(wp) :: workcoords1(3, mol1%natom)

   integer i
!   integer h, i, n, offset
!   real(wp), dimension(3, natom) :: randcoords0, randcoords1
!   integer :: votes(natom, natom)

!  Temporary variable copies

   natom = mol0%natom
   weights = mol0%get_weights()
   coords0 = mol0%get_atomcoords()
   coords1 = mol1%get_atomcoords()
   adjmat0 = mol0%get_adjmat()
   adjmat1 = mol1%get_adjmat()

   natomtype = mol0%get_natomtype()
   atomtypelenlist = mol0%get_atomtypelenlist()

   nadjs0 = mol0%get_nadjs()
   nadjs1 = mol1%get_nadjs()
   adjlists0 = mol0%get_adjlists()
   adjlists1 = mol1%get_adjlists()

   natomequiv0 = mol0%get_natomequiv()
   natomequiv1 = mol1%get_natomequiv()
   atomequivlenlist0 = mol0%get_atomequivlenlist()
   atomequivlenlist1 = mol1%get_atomequivlenlist()

   nadjequivs0 = mol0%get_nadjequivs()
   nadjequivs1 = mol1%get_nadjequivs()
   adjequivlenlists0 = mol0%get_adjequivlenlists()
   adjequivlenlists1 = mol1%get_adjequivlenlists()

   ! Detect fragments and root atoms

!CZGC: new call
!   call findmolfrags(mol0, nadjs0, adjlists0, nfrag0, fragroot0)
!   call findmolfrags(mol1, nadjs1, adjlists1, nfrag1, fragroot1)
!CZGC: optionally
!   call findmolfrags(mol0, nfrag0, fragroot0)
!   call findmolfrags(mol1, nfrag1, fragroot1)
!CZGC: optionally 2 (nota: sólo pasan todos los tests si se llama findmolfrags
!                          aquí y no al inicio en fileio.f90)
   call findmolfrags(mol0)
   call findmolfrags(mol1)
   nfrag0 = mol0%nfrag
   fragroot0 = mol0%get_fragroot()
   nfrag1 = mol1%nfrag
   fragroot1 = mol1%get_fragroot()
!CZGC: old call
!   call findmolfrags(natom, nadjs0, adjlists0, natomtype, atomtypelenlist, natomequiv0, atomequivlenlist0, nfrag0, fragroot0)
!   call findmolfrags(natom, nadjs1, adjlists1, natomtype, atomtypelenlist, natomequiv1, atomequivlenlist1, nfrag1, fragroot1)

   ! Calculate MNA equivalence matrix

!CZGC: new call
   call calcequivmat(mol0, mol1, nadjmna0, adjmnalen0, adjmnalist0, nadjmna1, adjmnalen1, adjmnalist1, equivmat)
!CZGC: old call
!   call calcequivmat(natom, natomtype, atomtypelenlist, nadjs0, adjlists0, nadjmna0, adjmnalen0, adjmnalist0, &
!      nadjs1, adjlists1, nadjmna1, adjmnalen1, adjmnalist1, equivmat)

   ! Print equivalence matrix

!   offset = 0
!   do h = 1, natomtype
!      print *, h
!      do i = offset + 1, offset + atomtypelenlist(h)
!         print '(100(1x, i3))', equivmat(offset+1:offset+atomtypelenlist(h), i)
!      end do
!      offset = offset + atomtypelenlist(h)
!      print *
!   end do

   ! Calculate bias matrix


!CZGC: new call
!   call setcrossbias(mol0, mol1, equivmat, biasmat)
!CZGC: old call
   call setcrossbias(natom, natomtype, atomtypelenlist, coords0, coords1, equivmat, biasmat)

   ! Print biases

!   offset = 0
!   do h = 1, blockcount0
!      do i = offset+1, offset+atomtypelenlist(h)
!         write (stdout, '(i0,":")', advance='no') i
!         do j = offset+1, offset+atomtypelenlist(h)
!            if (biasmat(i, j) == 0) then
!               write (stdout, '(1x,i0)', advance='no') j
!            end if
!         end do
!         print *
!      end do
!      offset = offset + atomtypelenlist(h)
!   end do

   ! Calculate assignment frequencies

!   votes(:, :) = 0
!   do n = 1, 100
!      do i = 1, natom
!         randcoords0(:, i) = randarray(3)
!         randcoords1(:, i) = randarray(3)
!      end do
!      call mapatoms(natom, randcoords0, randcoords1, natomtype, atomtypelenlist, biasmat, mapping)
!      print *, adjacencydiff(natom, adjmat0, adjmat1, mapping)
!      do i = 1, natom
!         votes(mapping(i), i) = votes(mapping(i), i) + 1
!      end do
!   end do
!   print *
!   offset = 0
!   do h = 1, natomtype
!      print *, h
!      do i = offset + 1, offset + atomtypelenlist(h)
!         print '(100(1x, i3))', votes(offset+1:offset+atomtypelenlist(h), i)
!      end do
!      offset = offset + atomtypelenlist(h)
!      print *
!   end do

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

      call mapatoms(natom, natomtype, atomtypelenlist, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
         nadjmna1, adjmnalen1, adjmnalist1, workcoords1, weights, equivmat, biasmat, mapping)
      rotquat = leastrotquat(natom, weights, coords0, workcoords1, mapping)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate(natom, workcoords1, rotquat)

      steps = 1

      do while (iter_flag)
         call mapatoms(natom, natomtype, atomtypelenlist, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
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

         call minadjdiff(natom, weights, natomtype, atomtypelenlist, coords0, nadjs0, adjlists0, adjmat0, natomequiv0, &
            atomequivlenlist0, workcoords1, nadjs1, adjlists1, adjmat1, natomequiv1, atomequivlenlist1, mapping, nfrag0, fragroot0)

         call eqvatomperm(natom, weights, coords0, adjmat0, adjlists0, natomequiv0, atomequivlenlist0, &
            nadjequivs0, adjequivlenlists0, workcoords1, adjmat1, mapping, nfrag0, fragroot0)

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

   if (stats_flag) then
      call print_stats(nrec, countlist, avgsteps, avgtotalrot, avgrealrot, recadjd, recrmsd)
      call print_final_stats(overflow, maxrec, nrec, ntrial, nstep)
   end if

end subroutine

subroutine find_reactive_sites(mol0, mol1, mapping)

   type(Molecule), intent(inout) :: mol0, mol1
   integer, dimension(:), intent(in) :: mapping

   integer :: natom
   integer :: h, i, j, k, m, n
   integer, dimension(:), allocatable :: atomnums0, atomnums1
   real(wp), dimension(:, :), allocatable :: coords0, coords1
   logical, dimension(:, :), allocatable :: adjmat0, adjmat1
   integer, dimension(:), allocatable :: atomtypeidcs0, atomtypeidcs1
   integer, dimension(:), allocatable :: atomtypelenlist0, atomtypelenlist1
   type(Part), dimension(:), allocatable :: atomtypeblocks0, atomtypeblocks1
   real(wp) :: rotquat(4)

   ! Align coordinates

   rotquat = leastrotquat(mol0%natom, mol0%get_weights(), mol0%get_atomcoords(), mol1%get_atomcoords(), mapping)
   call mol1%rotate_coords(rotquat)

   ! Initialization

   coords0 = mol0%get_atomcoords()
   coords1 = mol1%get_atomcoords()
   adjmat0 = mol0%get_adjmat()
   adjmat1 = mol1%get_adjmat()
   atomtypeidcs0 = mol0%get_atomtypeidcs()
   atomtypeidcs1 = mol1%get_atomtypeidcs()
   atomtypeblocks0 = mol0%get_atomtypeblocks()
   atomtypeblocks1 = mol1%get_atomtypeblocks()
   atomtypelenlist0 = mol0%get_atomtypelenlist()
   atomtypelenlist1 = mol1%get_atomtypelenlist()

   ! Remove reactive bonds

   do i = 1, mol0%natom
      m = atomtypeidcs0(i)
      do j = i + 1, mol0%natom
         n = atomtypeidcs0(j)
         if (adjmat0(i, j) .neqv. adjmat1(mapping(i), mapping(j))) then
            call remove_reactive_bond(i, j, mol0, mol1, mapping)
            do k = 1, atomtypelenlist0(m)
               if (sum((coords0(:, i) - coords1(:, mapping(atomtypeblocks1(m)%atomidcs(k))))**2) < 2.0) then
                  call remove_reactive_bond(atomtypeblocks1(m)%atomidcs(k), j, mol0, mol1, mapping)
               end if
               if (sum((coords0(:, atomtypeblocks0(m)%atomidcs(k)) - coords1(:, mapping(i)))**2) < 2.0) then
                  call remove_reactive_bond(atomtypeblocks0(m)%atomidcs(k), j, mol0, mol1, mapping)
               end if
            end do
            do k = 1, atomtypelenlist1(n)
               if (sum((coords0(:, j) - coords1(:, mapping(atomtypeblocks1(n)%atomidcs(k))))**2) < 2.0) then
                  call remove_reactive_bond(i, atomtypeblocks1(n)%atomidcs(k), mol0, mol1, mapping)
               end if
               if (sum((coords0(:, atomtypeblocks0(n)%atomidcs(k)) - coords1(:, mapping(j)))**2) < 2.0) then
                  call remove_reactive_bond(i, atomtypeblocks0(n)%atomidcs(k), mol0, mol1, mapping)
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
   type(Molecule), intent(inout) :: mol0, mol1

   integer :: k
   integer, allocatable, dimension(:) :: atomnums0, atomnums1
   integer :: nadjs0(mol0%natom), adjlists0(maxcoord, mol0%natom)
   integer :: nadjs1(mol1%natom), adjlists1(maxcoord, mol1%natom)

   nadjs0 = mol0%get_nadjs()
   nadjs1 = mol1%get_nadjs()
   adjlists0 = mol0%get_adjlists()
   adjlists1 = mol1%get_adjlists()
   atomnums0 = mol0%get_atomnums()
   atomnums1 = mol1%get_atomnums()

   if (mol0%bonded(i, j)) then
      call mol0%remove_bond(i, j)
   end if

   if (mol1%bonded(mapping(i), mapping(j))) then
      call mol1%remove_bond(mapping(i), mapping(j))
   end if

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
            if (mol1%bonded(mapping(i), mapping(adjlists1(k, i)))) then
               call mol1%remove_bond(mapping(i), mapping(adjlists1(k, i)))
            end if
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
            if (mol1%bonded(mapping(j), mapping(adjlists1(k, j)))) then
               call mol1%remove_bond(mapping(j), mapping(adjlists1(k, j)))
            end if
         end if
      end do
   end if

end subroutine

end module
