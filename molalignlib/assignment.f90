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
   blksz, &
   blkwt, &
   neqv0, &
   eqvsz0, &
   coords0, &
   adjmat0, &
   neqv1, &
   eqvsz1, &
   coords1, &
   adjmat1, &
   permlist, &
   countlist, &
   nrec)

   integer, intent(in) :: natom, nblk, neqv0, neqv1
   integer, dimension(:), intent(in) :: blksz
   integer, dimension(:), intent(in) :: eqvsz0, eqvsz1
   real(wp), dimension(:), intent(in) :: blkwt
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   logical, dimension(:, :), intent(in) :: adjmat0, adjmat1
   integer, intent(out) :: permlist(:, :)
   integer, intent(out) :: countlist(:)
   integer, intent(out) :: nrec

   logical visited, overflow
   integer :: nfrag0, nfrag1
   integer irec, mrec, ntrial, nstep, steps
   integer, dimension(natom) :: atomperm, auxperm
   integer :: adjdiff, recadjdiff(maxrec)
   integer, dimension(natom) :: nadj0, nadj1
   integer, dimension(maxcoord, natom) :: adjlist0, adjlist1
   integer, dimension(natom) :: neqvnei0, neqvnei1
   integer, dimension(maxcoord, natom) :: eqvneisz0, eqvneisz1
   integer, dimension(natom) :: blkid
   integer, dimension(natom) :: fragid0, fragid1
   integer, dimension(natom) :: fragroot0, fragroot1
   integer :: h, offset
   real(wp) :: dist2, olddist, newdist, totalrot
   real(wp), dimension(4) :: rotquat, prodquat
   real(wp), dimension(maxrec) :: recdist2, avgsteps, avgtotalrot, avgrealrot
   real(wp) :: biasmat(natom, natom)
   real(wp) :: workcoords1(3, natom)
   real(wp) :: weights(natom)

   ! Assign id's and weights to atoms

   offset = 0
   do h = 1, nblk
      blkid(offset+1:offset+blksz(h)) = h
      weights(offset+1:offset+blksz(h)) = blkwt(h)
      offset = offset + blksz(h)
   end do

   ! Recalculate adjacency lists

   call adjmat2list(natom, adjmat0, nadj0, adjlist0)
   call adjmat2list(natom, adjmat1, nadj1, adjlist1)

   ! Group equivalent neighbors

   call groupeqvnei(natom, neqv0, eqvsz0, nadj0, adjlist0, neqvnei0, eqvneisz0)
   call groupeqvnei(natom, neqv1, eqvsz1, nadj1, adjlist1, neqvnei1, eqvneisz1)

   ! Detect fagments and starting atoms

   call getmolfrags(natom, nadj0, adjlist0, blksz, blkid, neqv0, eqvsz0, nfrag0, fragroot0, fragid0)
   call getmolfrags(natom, nadj1, adjlist1, blksz, blkid, neqv1, eqvsz1, nfrag1, fragroot1, fragid1)

   ! Calculate biases

   call bias_func(natom, nblk, blksz, nadj0, adjlist0, nadj1, adjlist1, coords0, coords1, biasmat)

   ! Print biases

!   offset = 0
!   do h = 1, nblk0
!      do i = offset+1, offset+blksz(h)
!         write (output_unit, '(i0,":")', advance='no') i
!         do j = offset+1, offset+blksz(h)
!            if (biasmat(i, j) == 0) then
!               write (output_unit, '(1x,i0)', advance='no') j
!            end if
!         end do
!         print *
!      end do
!      offset = offset + blksz(h)
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

      call rotate(natom, workcoords1, genrotquat(rand3()))

   ! Minimize the euclidean distance

      call minatomperm(natom, coords0, workcoords1, nblk, blksz, blkwt, biasmat, atomperm, olddist)
      rotquat = leastrotquat(natom, weights, coords0, workcoords1, atomperm)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate(natom, workcoords1, rotquat)

      steps = 1

      do while (iter_flag)
         olddist = squaredist(natom, weights, coords0, workcoords1, atomperm) + biasdist(natom, weights, biasmat, atomperm)
         call minatomperm(natom, coords0, workcoords1, nblk, blksz, blkwt, biasmat, auxperm, newdist)
         if (all(auxperm == atomperm)) exit
         if (newdist > olddist) then
            write (error_unit, '(a)') 'newdist is larger than olddist!'
!            print *, olddist, newdist
         end if
         rotquat = leastrotquat(natom, weights, coords0, workcoords1, auxperm)
         prodquat = quatmul(rotquat, prodquat)
         call rotate(natom, workcoords1, rotquat)
         totalrot = totalrot + rotangle(rotquat)
         atomperm = auxperm
         steps = steps + 1
      end do

      nstep = nstep + steps

!      dist2 = squaredist(natom, weights, coords0, workcoords1, atomperm)

      call minadjdiff(natom, weights, blkid, coords0, nadj0, adjlist0, adjmat0, &
         coords1, nadj1, adjlist1, adjmat1, atomperm, nfrag0, fragroot0)

      call eqvatomperm(natom, weights, coords0, adjmat0, adjlist0, neqv0, eqvsz0, &
         neqvnei0, eqvneisz0, workcoords1, adjmat1, atomperm, nfrag0, fragroot0)

      dist2 = leastsquaredist(natom, weights, coords0, coords1, atomperm)
      adjdiff = adjacencydiff(natom, adjmat0, adjmat1, atomperm)

      ! Check for new best permlist

      visited = .false.

      do irec = 1, nrec
         if (all(atomperm == permlist(:, irec))) then
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
            if (adjdiff < recadjdiff(irec) .or. (adjdiff == recadjdiff(irec) .and. dist2 < recdist2(irec))) then
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
               recdist2(irec) = recdist2(irec - 1)
               recadjdiff(irec) = recadjdiff(irec - 1)
               avgsteps(irec) = avgsteps(irec - 1)
               avgrealrot(irec) = avgrealrot(irec - 1)
               avgtotalrot(irec) = avgtotalrot(irec - 1)
               permlist(:, irec) = permlist(:, irec - 1)
            end do
            countlist(mrec) = 1
            recdist2(mrec) = dist2
            recadjdiff(mrec) = adjdiff
            avgsteps(mrec) = steps
            avgrealrot(mrec) = rotangle(prodquat)
            avgtotalrot(mrec) = totalrot
            permlist(:, mrec) = atomperm
         end if
      end if

   end do

   if (stats_flag) then
      call print_stats(nrec, weights, countlist, avgsteps, avgtotalrot, avgrealrot, recadjdiff, recdist2)
      call print_final_stats(overflow, maxrec, nrec, ntrial, nstep)
   end if

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine minatomperm(natom, coords0, coords1, nblk, blksz, blkwt, biasmat, atomperm, totdist)

! nblk: Number of block atoms
! blksz: Number of atoms in each block
! atomperm: Map between correspondent points in the adjmat
! offset: First element of current block

   integer, intent(in) :: natom, nblk
   integer, dimension(:), intent(in) :: blksz
   real(wp), dimension(:), intent(in) :: blkwt
   real(wp), dimension(:, :), intent(in) :: coords0
   real(wp), dimension(:, :), intent(in) :: coords1
   real(wp), dimension(:, :), intent(in) :: biasmat
   integer, dimension(:), intent(out) :: atomperm
   real(wp), intent(out) :: totdist

   integer :: h, offset
   integer, dimension(natom) :: perm
   real(wp) :: dist

   ! Fill distance matrix for each block

   offset = 0
   totdist = 0

   do h = 1, nblk
      call minperm(blksz(h), coords0(:, offset+1:offset+blksz(h)), coords1(:, offset+1:offset+blksz(h)), &
         biasmat(offset+1:offset+blksz(h), offset+1:offset+blksz(h)), perm, dist)
      atomperm(offset+1:offset+blksz(h)) = perm(:blksz(h)) + offset
      totdist = totdist + blkwt(h)*dist
      offset = offset + blksz(h)
   end do

end subroutine

end module
