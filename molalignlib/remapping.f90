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
public optimize_assignment

contains

subroutine optimize_assignment( &
   natom, &
   nblk, &
   blklen, &
   neqv0, &
   eqvlen0, &
   coords0, &
   adjmat0, &
   neqv1, &
   eqvlen1, &
   coords1, &
   adjmat1, &
   weights, &
   maplist, &
   countlist, &
   nrec)

   integer, intent(in) :: natom, nblk, neqv0, neqv1
   integer, dimension(:), intent(in) :: blklen
   integer, dimension(:), intent(in) :: eqvlen0, eqvlen1
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   logical, dimension(:, :), intent(in) :: adjmat0, adjmat1
   real(wp), dimension(:), intent(in) :: weights
   integer, intent(out) :: maplist(:, :)
   integer, intent(out) :: countlist(:)
   integer, intent(out) :: nrec

   logical visited, overflow
   integer :: nfrag0, nfrag1
   integer irec, mrec, ntrial, nstep, steps
   integer, dimension(natom) :: mapping, newmapping
   integer :: adjd, recadjd(maxrec)
   integer, dimension(natom) :: nadj0, nadj1
   integer, dimension(maxcoord, natom) :: adjlist0, adjlist1
   integer, dimension(natom, maxlevel) :: nadjmna0, nadjmna1
   integer, dimension(maxcoord, natom, maxlevel) :: adjmnalen0, adjmnalen1
   integer, dimension(maxcoord, natom, maxlevel) :: adjmnalist0, adjmnalist1
   integer, dimension(natom) :: nadjeqv0, nadjeqv1
   integer, dimension(maxcoord, natom) :: adjeqvlen0, adjeqvlen1
   integer, dimension(natom) :: fragroot0, fragroot1
   integer :: equivmat(natom, natom)
   real(wp) :: rmsd, totalrot
   real(wp), dimension(4) :: rotquat, prodquat
   real(wp), dimension(maxrec) :: recrmsd, avgsteps, avgtotalrot, avgrealrot
   real(wp) :: biasmat(natom, natom)
   real(wp) :: workcoords1(3, natom)

!   integer h, i, n, offset
!   real(wp), dimension(3, natom) :: randcoords0, randcoords1
!   integer :: votes(natom, natom)

   ! Recalculate adjacency lists

   call adjmat2list(natom, adjmat0, nadj0, adjlist0)
   call adjmat2list(natom, adjmat1, nadj1, adjlist1)

   ! Group neighbors by type

!   call groupneighbors(natom, nblk, blklen, nadj0, adjlist0, nadjmna0, adjmnalen0)
!   call groupneighbors(natom, nblk, blklen, nadj1, adjlist1, nadjmna1, adjmnalen1)

   ! Group neighbors by MNA

   call groupneighbors(natom, neqv0, eqvlen0, nadj0, adjlist0, nadjeqv0, adjeqvlen0)
   call groupneighbors(natom, neqv1, eqvlen1, nadj1, adjlist1, nadjeqv1, adjeqvlen1)

   ! Print adjacency lists

!   do i = 1, natom
!      offset = 0
!      do h = 1, nadjmna0(i)
!         do k = offset + 1, offset + adjmnalen0(h, i)
!            print *, i, h, adjlist0(k, i)
!         end do
!         offset = offset + adjmnalen0(h, i)
!      end do
!   end do
!
!   do i = 1, natom
!      offset = 0
!      do h = 1, nadjmna1(i)
!         do k = offset + 1, offset + adjmnalen1(h, i)
!            print *, i, h, adjlist1(k, i)
!         end do
!         offset = offset + adjmnalen1(h, i)
!      end do
!   end do

   ! Detect fagments and starting atoms

   call findmolfrags(natom, nadj0, adjlist0, nblk, blklen, neqv0, eqvlen0, nfrag0, fragroot0)
   call findmolfrags(natom, nadj1, adjlist1, nblk, blklen, neqv1, eqvlen1, nfrag1, fragroot1)

   ! Calculate MNA equivalence matrix

   call calcequivmat(natom, nblk, blklen, nadj0, adjlist0, nadjmna0, adjmnalen0, adjmnalist0, &
      nadj1, adjlist1, nadjmna1, adjmnalen1, adjmnalist1, equivmat)

   ! Print equivalence matrix

!   offset = 0
!   do h = 1, nblk
!      print *, h
!      do i = offset + 1, offset + blklen(h)
!         print '(100(1x, i3))', equivmat(offset+1:offset+blklen(h), i)
!      end do
!      offset = offset + blklen(h)
!      print *
!   end do

   ! Calculate bias matrix

   call setcrossbias(natom, nblk, blklen, coords0, coords1, equivmat, biasmat)

   ! Print biases

!   offset = 0
!   do h = 1, nblk0
!      do i = offset+1, offset+blklen(h)
!         write (stdout, '(i0,":")', advance='no') i
!         do j = offset+1, offset+blklen(h)
!            if (biasmat(i, j) == 0) then
!               write (stdout, '(1x,i0)', advance='no') j
!            end if
!         end do
!         print *
!      end do
!      offset = offset + blklen(h)
!   end do

   ! Calculate assignment frequencies

!   votes(:, :) = 0
!   do n = 1, 100
!      do i = 1, natom
!         randcoords0(:, i) = randarray(3)
!         randcoords1(:, i) = randarray(3)
!      end do
!      call mapatoms(natom, randcoords0, randcoords1, nblk, blklen, biasmat, mapping)
!      print *, adjacencydiff(natom, adjmat0, adjmat1, mapping)
!      do i = 1, natom
!         votes(mapping(i), i) = votes(mapping(i), i) + 1
!      end do
!   end do
!   print *
!   offset = 0
!   do h = 1, nblk
!      print *, h
!      do i = offset + 1, offset + blklen(h)
!         print '(100(1x, i3))', votes(offset+1:offset+blklen(h), i)
!      end do
!      offset = offset + blklen(h)
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

      call mapatoms(natom, nblk, blklen, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
         nadjmna1, adjmnalen1, adjmnalist1, workcoords1, weights, equivmat, biasmat, mapping)
      rotquat = leastrotquat(natom, weights, coords0, workcoords1, mapping)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate(natom, workcoords1, rotquat)

      steps = 1

      do while (iter_flag)
         call mapatoms(natom, nblk, blklen, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
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

         call minadjdiff(natom, weights, nblk, blklen, coords0, nadj0, adjlist0, adjmat0, neqv0, &
            eqvlen0, workcoords1, nadj1, adjlist1, adjmat1, neqv1, eqvlen1, mapping, nfrag0, fragroot0)

         call eqvatomperm(natom, weights, coords0, adjmat0, adjlist0, neqv0, eqvlen0, &
            nadjeqv0, adjeqvlen0, workcoords1, adjmat1, mapping, nfrag0, fragroot0)

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

end module
