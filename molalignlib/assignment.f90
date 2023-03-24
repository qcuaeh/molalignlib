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

module assignment
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
use lap

implicit none

private
public optimize_assignment

contains

subroutine optimize_assignment( &
   natom, &
   nblk, &
   blklen, &
   blkwgt, &
   neqv0, &
   eqvlen0, &
   coords0, &
   adjmat0, &
   neqv1, &
   eqvlen1, &
   coords1, &
   adjmat1, &
   maplist, &
   countlist, &
   nrec)

   integer, intent(in) :: natom, nblk, neqv0, neqv1
   integer, dimension(:), intent(in) :: blklen
   integer, dimension(:), intent(in) :: eqvlen0, eqvlen1
   real(wp), dimension(:), intent(in) :: blkwgt
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   logical, dimension(:, :), intent(in) :: adjmat0, adjmat1
   integer, intent(out) :: maplist(:, :)
   integer, intent(out) :: countlist(:)
   integer, intent(out) :: nrec

   integer :: h, offset
   logical visited, overflow
   integer :: nfrag0, nfrag1
   integer irec, mrec, ntrial, nstep, steps
   integer, dimension(natom) :: mapping, auxperm
   integer :: adjd, recadjd(maxrec)
   integer, dimension(natom) :: nadj0, nadj1
   integer, dimension(maxcoord, natom) :: adjlist0, adjlist1
   integer, dimension(natom) :: nadjblk0, nadjblk1
   integer, dimension(maxcoord, natom) :: adjblklen0, adjblklen1
   integer, dimension(natom) :: neqvnei0, neqvnei1
   integer, dimension(maxcoord, natom) :: eqvneilen0, eqvneilen1
   integer, dimension(natom) :: fragroot0, fragroot1
   integer :: equivmat(natom, natom)
   real(wp) :: rmsd, olddist, newdist, totalrot
   real(wp), dimension(4) :: rotquat, prodquat
   real(wp), dimension(maxrec) :: recrmsd, avgsteps, avgtotalrot, avgrealrot
   real(wp) :: biasmat(natom, natom)
   real(wp) :: workcoords1(3, natom)
   real(wp) :: weights(natom)

   integer i, n
   real(wp), dimension(3, natom) :: randcoords0, randcoords1
   integer :: votes(natom, natom)

   ! Assign id's and weights to atoms

   offset = 0
   do h = 1, nblk
      weights(offset+1:offset+blklen(h)) = blkwgt(h)
      offset = offset + blklen(h)
   end do

   ! Recalculate adjacency lists

   call adjmat2list(natom, adjmat0, nadj0, adjlist0)
   call adjmat2list(natom, adjmat1, nadj1, adjlist1)

   ! Group neighbors by type

   call groupneighbors(natom, nblk, blklen, nadj0, adjlist0, nadjblk0, adjblklen0)
   call groupneighbors(natom, nblk, blklen, nadj1, adjlist1, nadjblk1, adjblklen1)

   ! Group neighbors by MNA

   call groupneighbors(natom, neqv0, eqvlen0, nadj0, adjlist0, neqvnei0, eqvneilen0)
   call groupneighbors(natom, neqv1, eqvlen1, nadj1, adjlist1, neqvnei1, eqvneilen1)

   ! Detect fagments and starting atoms

   call findmolfrags(natom, nadj0, adjlist0, nblk, blklen, neqv0, eqvlen0, nfrag0, fragroot0)
   call findmolfrags(natom, nadj1, adjlist1, nblk, blklen, neqv1, eqvlen1, nfrag1, fragroot1)

   ! Calculate MNA equivalence matrix

   call calcequivmat(natom, nblk, blklen, nadj0, adjlist0, nadj1, adjlist1, equivmat)

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

   call calcbiasmat(natom, nblk, blklen, coords0, coords1, equivmat, biasmat)

   ! Print biases

!   offset = 0
!   do h = 1, nblk0
!      do i = offset+1, offset+blklen(h)
!         write (output_unit, '(i0,":")', advance='no') i
!         do j = offset+1, offset+blklen(h)
!            if (biasmat(i, j) == 0) then
!               write (output_unit, '(1x,i0)', advance='no') j
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
!      call minatomperm(natom, randcoords0, randcoords1, nblk, blklen, blkwgt, biasmat, mapping, olddist)
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

      call minatomperm(natom, coords0, workcoords1, nblk, blklen, blkwgt, biasmat, mapping, olddist)
!      call mapatoms(natom, nblk, blklen, blkwgt, nadjblk0, adjblklen0, adjlist0, coords0, adjlist1, &
!         coords1, equivmat, mapping, olddist)
      rotquat = leastrotquat(natom, weights, coords0, workcoords1, mapping)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate(natom, workcoords1, rotquat)

      steps = 1

      do while (iter_flag)
         olddist = squaredist(natom, weights, coords0, workcoords1, mapping) + biasdist(natom, weights, biasmat, mapping)
         call minatomperm(natom, coords0, workcoords1, nblk, blklen, blkwgt, biasmat, auxperm, newdist)
!         call mapatoms(natom, nblk, blklen, blkwgt, nadjblk0, adjblklen0, adjlist0, coords0, adjlist1, &
!            coords1, equivmat, mapping, newdist)
         if (all(auxperm == mapping)) exit
         if (newdist > olddist) then
            write (error_unit, '(a)') 'newdist is larger than olddist!'
!            print *, olddist, newdist
         end if
         rotquat = leastrotquat(natom, weights, coords0, workcoords1, auxperm)
         prodquat = quatmul(rotquat, prodquat)
         call rotate(natom, workcoords1, rotquat)
         totalrot = totalrot + rotangle(rotquat)
         mapping = auxperm
         steps = steps + 1
      end do

      nstep = nstep + steps

      if (back_flag) then

         call minadjdiff(natom, weights, nblk, blklen, coords0, nadj0, adjlist0, adjmat0, neqv0, &
            eqvlen0, workcoords1, nadj1, adjlist1, adjmat1, neqv1, eqvlen1, mapping, nfrag0, fragroot0)

         call eqvatomperm(natom, weights, coords0, adjmat0, adjlist0, neqv0, eqvlen0, &
            neqvnei0, eqvneilen0, workcoords1, adjmat1, mapping, nfrag0, fragroot0)

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
