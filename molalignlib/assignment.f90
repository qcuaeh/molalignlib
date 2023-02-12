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
use lap
use translation
use alignment
use rotation
use biasing
use printing

implicit none

private
public optimize_assignment

contains

subroutine optimize_assignment( &
   natom, &
   nblk, &
   blksz, &
   weights, &
   coords0, &
   coords1, &
   permlist, &
   countlist, &
   nrec)

   integer, intent(in) :: natom, nblk
   integer, dimension(:), intent(in) :: blksz
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   real(wp), dimension(:), intent(in) :: weights
   integer, intent(out) :: permlist(:, :)
   integer, intent(out) :: countlist(:)
   integer, intent(out) :: nrec

   logical visited, overflow
   integer irec, ntrial, nstep, steps
   integer, dimension(natom) :: atomperm, auxmap
   real(wp) :: dist2, olddist, newdist, totalrot
   real(wp), dimension(4) :: rotquat, prodquat
   real(wp), dimension(maxrec) :: dist2rec, avgsteps, avgtotalrot, avgrealrot
   real(wp), dimension(natom, natom) :: biasmat
   real(wp) :: workcoords1(3, natom)

! Calculate biases

   if (bias_flag) then
      call setsdnbias(natom, nblk, blksz, coords0, coords1, biasmat)
   else
      biasmat(:, :) = 0
   end if

! Print header and initial stats

   if (stats_flag .and. live_flag) then
      write (error_unit, '(a)', advance='no') achar(27) // '[1H' // achar(27) // '[J'
      call print_header()
   end if

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

      call minatomperm(natom, coords0, workcoords1, nblk, blksz, biasmat, weights, atomperm, olddist)
      rotquat = leastrotquat(natom, weights, coords0, workcoords1, atomperm)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate(natom, workcoords1, rotquat)

      steps = 1

      do while (iter_flag)
         olddist = biasdist(natom, atomperm, weights, coords0, workcoords1, biasmat)
         call minatomperm(natom, coords0, workcoords1, nblk, blksz, biasmat, weights, auxmap, newdist)
         if (all(auxmap == atomperm)) exit
         if (newdist > olddist) then
            write (error_unit, '(a)') 'newdist is larger than olddist!'
!            print *, olddist, newdist
         end if
         rotquat = leastrotquat(natom, weights, coords0, workcoords1, auxmap)
         prodquat = quatmul(rotquat, prodquat)
         call rotate(natom, workcoords1, rotquat)
         totalrot = totalrot + rotangle(rotquat)
         atomperm = auxmap
         steps = steps + 1
      end do

      nstep = nstep + steps

      dist2 = squaredist(natom, weights, coords0, workcoords1, atomperm)

! Check for new best permlist

      visited = .false.

      do irec = 1, nrec
         if (all(atomperm == permlist(:, irec))) then
            countlist(irec) = countlist(irec) + 1
            avgsteps(irec) = avgsteps(irec) + (steps - avgsteps(irec))/countlist(irec)
            avgrealrot(irec) = avgrealrot(irec) + (rotangle(prodquat) - avgrealrot(irec))/countlist(irec)
            avgtotalrot(irec) = avgtotalrot(irec) + (totalrot - avgtotalrot(irec))/countlist(irec)
            if (stats_flag .and. live_flag) then
               write (error_unit, '(a)', advance='no') achar(27) // '[' // intstr(irec + 2) // 'H'
               call print_body(irec, countlist(irec), avgsteps(irec), avgtotalrot(irec), &
                  avgrealrot(irec), dist2rec(irec))
            end if
            visited = .true.
            exit
         end if
      end do

      if (.not. visited) then
         if (nrec < maxrec) then
            nrec = nrec + 1
         else
            overflow = .true.
            if (dist2 > dist2rec(nrec)) cycle
         end if
         do irec = nrec, 2, -1
            if (dist2 > dist2rec(irec - 1)) exit
            permlist(:, irec) = permlist(:, irec - 1)
            countlist(irec) = countlist(irec - 1)
            dist2rec(irec) = dist2rec(irec - 1)
            avgsteps(irec) = avgsteps(irec - 1)
            avgrealrot(irec) = avgrealrot(irec - 1)
            avgtotalrot(irec) = avgtotalrot(irec - 1)
            if (stats_flag .and. live_flag) then
               write (error_unit, '(a)', advance='no') achar(27) // '[' // intstr(irec + 2) // 'H'
               call print_body(irec, countlist(irec), avgsteps(irec), avgtotalrot(irec), &
                  avgrealrot(irec), dist2rec(irec))
            end if
         end do
         permlist(:, irec) = atomperm
         countlist(irec) = 1
         dist2rec(irec) = dist2
         avgsteps(irec) = steps
         avgrealrot(irec) = rotangle(prodquat)
         avgtotalrot(irec) = totalrot
         if (stats_flag .and. live_flag) then
            write (error_unit, '(a)', advance='no') achar(27) // '[' // intstr(irec + 2) // 'H'
            call print_body(irec, countlist(irec), avgsteps(irec), avgtotalrot(irec), &
               avgrealrot(irec), dist2rec(irec))
         end if
      end if

      if (stats_flag .and. live_flag) then
         write (error_unit, '(a)', advance='no') achar(27) // '[' // intstr(nrec + 3) // 'H'
         call print_footer()
      end if

   end do

   if (stats_flag .and. .not. live_flag) then
      call print_header()
      do irec = 1, nrec
         call print_body(irec, countlist(irec), avgsteps(irec), avgtotalrot(irec), &
            avgrealrot(irec), dist2rec(irec))
      end do
      call print_footer()
   end if

   if (stats_flag) then
      call print_stats(overflow, maxrec, nrec, ntrial, nstep)
   end if

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine minatomperm(natom, coords0, coords1, nblk, blksz, biasmat, weights, atomperm, totdist)

! nblk: Number of block atoms
! blksz: Number of atoms in each block
! atomperm: Map between correspondent points in the adjmat
! offset: First element of current block

   integer, intent(in) :: natom, nblk
   integer, dimension(:), intent(in) :: blksz
   real(wp), dimension(:, :), intent(in) :: coords0
   real(wp), dimension(:, :), intent(in) :: coords1
   real(wp), dimension(:, :), intent(in) :: biasmat
   real(wp), dimension(:), intent(in) :: weights
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
      totdist = totdist + weights(h)*dist
      offset = offset + blksz(h)
   end do

end subroutine

end module
