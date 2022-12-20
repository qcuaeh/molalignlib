! MolAlign
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
use parameters
use settings
use math
use random
use strutils
use translation
use alignment
use rotation
use biasing
use printing
use lap

implicit none

private
public optimize_assignment

contains

subroutine optimize_assignment( &
   natom, &
   nblock, &
   bsize, &
   weights, &
   coords0, &
   coords1, &
   nrec, &
   nmap, &
   maplist, &
   countlist, &
   dist2list &
)

   integer, intent(in) :: natom, nblock, nrec
   integer, dimension(:), intent(in) :: bsize
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   real(wp), dimension(:), intent(in) :: weights
   integer, intent(out) :: nmap
   integer, intent(out) :: maplist(:, :)
   integer, intent(out) :: countlist(:)
   real(wp), intent(out) :: dist2list(:)

   logical visited, overflow
   integer imap, ntrial, nstep, steps
   integer, dimension(natom) :: atomap, auxmap
   real(wp) :: dist2, olddist, newdist, totalrot
   real(wp), dimension(4) :: rotquat, prodquat
   real(wp), dimension(nrec) :: avgsteps, avgtotalrot, avgrealrot
   real(wp) :: bias(natom, natom)
   real(wp) :: workcoords1(3, natom)

! Set bias for non equivalent atoms 

   call setadjbias(natom, nblock, bsize, coords0, coords1, bias)

! Print header and initial stats

   if (stats_flag .and. live_flag) then
      write (output_unit, '(a)', advance='no') achar(27)//'[1H'//achar(27)//'[J'
      call print_header()
   end if

! Initialize loop variables

   nmap = 0
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

      call rotate(natom, workcoords1, torotquat(rand3()))

! Minimize the euclidean distance

      call minatomap(natom, coords0, workcoords1, nblock, bsize, bias, weights, atomap, olddist)
      rotquat = leastrotquat(natom, weights, coords0, workcoords1, atomap)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate(natom, workcoords1, rotquat)

      steps = 1

      do while (iter_flag)
         olddist = biasedist(natom, atomap, weights, coords0, workcoords1, bias)
         call minatomap(natom, coords0, workcoords1, nblock, bsize, bias, weights, auxmap, newdist)
         if (all(auxmap == atomap)) exit
         if (newdist > olddist) then
            write (output_unit, '(a)') 'newdist is larger than olddist!'
!            print *, olddist, newdist
         end if
         rotquat = leastrotquat(natom, weights, coords0, workcoords1, auxmap)
         prodquat = quatmul(rotquat, prodquat)
         call rotate(natom, workcoords1, rotquat)
         totalrot = totalrot + rotangle(rotquat)
         atomap = auxmap
         steps = steps + 1
      end do

      nstep = nstep + steps

      dist2 = squaredist(natom, weights, coords0, workcoords1, atomap)

! Check for new best maplist

      visited = .false.

      do imap = 1, nmap
         if (all(atomap == maplist(:, imap))) then
            countlist(imap) = countlist(imap) + 1
            avgsteps(imap) = avgsteps(imap) + (steps - avgsteps(imap))/countlist(imap)
            avgrealrot(imap) = avgrealrot(imap) + (rotangle(prodquat) - avgrealrot(imap))/countlist(imap)
            avgtotalrot(imap) = avgtotalrot(imap) + (totalrot - avgtotalrot(imap))/countlist(imap)
            if (stats_flag .and. live_flag) then
               write (output_unit, '(a)', advance='no') achar(27)//'['//str(imap+2)//'H'
               call print_body(imap, countlist(imap), avgsteps(imap), avgtotalrot(imap), &
                  avgrealrot(imap), dist2list(imap))
            end if
            visited = .true.
            exit
         end if
      end do

      if (.not. visited) then
         if (nmap < nrec) then
            nmap = nmap + 1
         else
            overflow = .true.
            if (dist2 > dist2list(nmap)) cycle
         end if
         do imap = nmap, 2, -1
            if (dist2 > dist2list(imap - 1)) exit
            maplist(:, imap) = maplist(:, imap - 1)
            countlist(imap) = countlist(imap - 1)
            dist2list(imap) = dist2list(imap - 1)
            avgsteps(imap) = avgsteps(imap - 1)
            avgrealrot(imap) = avgrealrot(imap - 1)
            avgtotalrot(imap) = avgtotalrot(imap - 1)
            if (stats_flag .and. live_flag) then
               write (output_unit, '(a)', advance='no') achar(27)//'['//str(imap + 2)//'H'
               call print_body(imap, countlist(imap), avgsteps(imap), avgtotalrot(imap), &
                  avgrealrot(imap), dist2list(imap))
            end if
         end do
         maplist(:, imap) = atomap
         countlist(imap) = 1
         dist2list(imap) = dist2
         avgsteps(imap) = steps
         avgrealrot(imap) = rotangle(prodquat)
         avgtotalrot(imap) = totalrot
         if (stats_flag .and. live_flag) then
            write (output_unit, '(a)', advance='no') achar(27)//'['//str(imap + 2)//'H'
            call print_body(imap, countlist(imap), avgsteps(imap), avgtotalrot(imap), &
               avgrealrot(imap), dist2list(imap))
         end if
      end if

      if (stats_flag .and. live_flag) then
         write (output_unit, '(a)', advance='no') achar(27)//'['//str(nmap + 3)//'H'
         call print_footer()
      end if

   end do

   if (stats_flag .and. .not. live_flag) then
      call print_header()
      do imap = 1, nmap
         call print_body(imap, countlist(imap), avgsteps(imap), avgtotalrot(imap), &
            avgrealrot(imap), dist2list(imap))
      end do
      call print_footer()
   end if

   if (stats_flag) then
      call print_stats(overflow, nrec, nmap, ntrial, nstep)
   end if

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine minatomap(natom, coords0, coords1, nblock, bsize, bias, weights, atomap, totdist)

! nblock: Number of block atoms
! bsize: Number of atoms in each block
! atomap: Map between correspondent points in the adjmat
! offset: First element of current block

   integer, intent(in) :: natom, nblock
   integer, dimension(:), intent(in) :: bsize
   real(wp), dimension(:, :), intent(in) :: coords0
   real(wp), dimension(:, :), intent(in) :: coords1
   real(wp), dimension(:, :), intent(in) :: bias
   real(wp), dimension(:), intent(in) :: weights
   integer, dimension(:), intent(out) :: atomap
   real(wp), intent(out) :: totdist

   integer :: h, offset
   integer, dimension(natom) :: perm
   real(wp) :: dist

! Fill distance matrix for each block

   offset = 0
   totdist = 0

   do h = 1, nblock
      call minperm(bsize(h), coords0(:, offset+1:offset+bsize(h)), coords1(:, offset+1:offset+bsize(h)), &
         bias(offset+1:offset+bsize(h), offset+1:offset+bsize(h)), perm, dist)
      atomap(offset+1:offset+bsize(h)) = perm(:bsize(h)) + offset
      totdist = totdist + weights(h)*dist
      offset = offset + bsize(h)
   end do

end subroutine

! Calculate total bias
real(wp) function biasedist(natom, atomap, weights, coords0, coords1, bias) result(dist)

   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: atomap
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords0
   real(wp), dimension(:, :), intent(in) :: coords1
   real(wp), dimension(:, :), intent(in) :: bias
   integer :: i

   dist = 0.

   do i = 1, natom
      dist = dist + weights(i)*(sum((coords0(:, i) - coords1(:, atomap(i)))**2) + bias(i, atomap(i)))
   end do

end function

end module
