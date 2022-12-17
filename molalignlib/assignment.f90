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
   blocksize, &
   weights, &
   coords0, &
   coords1, &
   nrec, &
   nmap, &
   mapping, &
   mapcount, &
   mapdist2 &
)

   integer, intent(in) :: natom, nblock, nrec
   integer, dimension(:), intent(in) :: blocksize
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   real(wp), dimension(:), intent(in) :: weights
   integer, intent(out) :: nmap
   integer, intent(out) :: mapping(:, :)
   integer, intent(out) :: mapcount(:)
   real(wp), intent(out) :: mapdist2(:)

   logical found, overflow
   integer imap, jmap, ncount, ntrial, nstep, steps
   integer, dimension(natom) :: atomap, auxmap
   real(wp) :: dist2, olddist, newdist, totalrot
   real(wp), dimension(4) :: rotquat, prodquat
   real(wp), dimension(nrec) :: avgsteps, avgtotalrot, avgrealrot
   real(wp) :: bias(natom, natom)
   real(wp) :: workcoords1(3, natom)

! Set bias for non equivalent atoms 

   call setadjbias(natom, nblock, blocksize, coords0, coords1, bias)

! Print header and initial stats

   if (live_flag) then
      write (output_unit, '(a)', advance='no') achar(27)//'[1H'//achar(27)//'[J'
      call print_header()
   end if

! Initialize loop variables

   nmap = 0
   nstep = 0
   ntrial = 0
   ncount = 0
   overflow = .false.

! Loop for map searching

   do while (ncount < maxcount .and. (.not. trial_flag .or. ntrial < maxtrials))

      ntrial = ntrial + 1

! Work with a copy of coords1

      workcoords1 = coords1

! Aply a random rotation to workcoords1

      call rotate(natom, workcoords1, torotquat(rand3()))

! Minimize the euclidean distance

      call minatomap(natom, coords0, workcoords1, nblock, blocksize, bias, weights, atomap, olddist)
      rotquat = leastrotquat(natom, weights, coords0, workcoords1, atomap)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate(natom, workcoords1, rotquat)

      steps = 1

      do while (iter_flag)
         olddist = biasedist(natom, atomap, weights, coords0, workcoords1, bias)
         call minatomap(natom, coords0, workcoords1, nblock, blocksize, bias, weights, auxmap, newdist)
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

! Check for new best mapping

      found = .false.

      do imap = 1, nmap
         if (all(atomap == mapping(:, imap))) then
            if (imap == 1) ncount = ncount + 1
            mapcount(imap) = mapcount(imap) + 1
            avgsteps(imap) = avgsteps(imap) + (steps - avgsteps(imap))/mapcount(imap)
            avgrealrot(imap) = avgrealrot(imap) + (rotangle(prodquat) - avgrealrot(imap))/mapcount(imap)
            avgtotalrot(imap) = avgtotalrot(imap) + (totalrot - avgtotalrot(imap))/mapcount(imap)
            if (live_flag) then
               write (output_unit, '(a)', advance='no') achar(27)//'['//str(imap + 2)//'H'
               call print_body(imap, mapcount(imap), avgsteps(imap), avgtotalrot(imap), &
                  avgrealrot(imap), mapdist2(imap))
            end if
            found = .true.
            exit
         end if
      end do

      if (.not. found) then
         if (nmap >= nrec) then
            overflow = .true.
         end if
         do imap = 1, nrec
            if (imap > nmap .or. dist2 < mapdist2(imap)) then
               if (imap == 1) ncount = 1
               if (nmap < nrec) nmap = nmap + 1
               do jmap = nmap, imap + 1, -1
                  mapping(:, jmap) = mapping(:, jmap - 1)
                  mapcount(jmap) = mapcount(jmap - 1)
                  mapdist2(jmap) = mapdist2(jmap - 1)
                  avgsteps(jmap) = avgsteps(jmap - 1)
                  avgrealrot(jmap) = avgrealrot(jmap - 1)
                  avgtotalrot(jmap) = avgtotalrot(jmap - 1)
                  if (live_flag) then
                     write (output_unit, '(a)', advance='no') achar(27)//'['//str(jmap + 2)//'H'
                     call print_body(jmap, mapcount(jmap), avgsteps(jmap), avgtotalrot(jmap), &
                        avgrealrot(jmap), mapdist2(jmap))
                  end if
               end do
               mapping(:, imap) = atomap
               mapcount(imap) = 1
               mapdist2(imap) = dist2
               avgsteps(imap) = steps
               avgrealrot(imap) = rotangle(prodquat)
               avgtotalrot(imap) = totalrot
               if (live_flag) then
                  write (output_unit, '(a)', advance='no') achar(27)//'['//str(imap + 2)//'H'
                  call print_body(imap, mapcount(imap), avgsteps(imap), avgtotalrot(imap), &
                     avgrealrot(imap), mapdist2(imap))
               end if
               exit
            end if
         end do
      end if

      if (live_flag) then
         write (output_unit, '(a)', advance='no') achar(27)//'['//str(nmap + 3)//'H'
         call print_footer()
      end if

   end do

   if (.not. live_flag) then
      call print_header()
      do imap = 1, nmap
         call print_body(imap, mapcount(imap), avgsteps(imap), avgtotalrot(imap), &
            avgrealrot(imap), mapdist2(imap))
      end do
      call print_footer()
   end if

   call print_stats(overflow, nrec, nmap, ntrial, nstep)

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine minatomap(natom, coords0, coords1, nblock, blocksize, bias, weights, atomap, totdist)

! nblock: Number of block atoms
! blocksize: Number of atoms in each block
! atomap: Map between correspondent points in the adjmat
! offset: First element of current block

   integer, intent(in) :: natom, nblock
   integer, dimension(:), intent(in) :: blocksize
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
      call minperm(blocksize(h), offset, coords0, coords1, bias, perm, dist)
      atomap(offset+1:offset+blocksize(h)) = perm(:blocksize(h)) + offset
      totdist = totdist + weights(h)*dist
      offset = offset + blocksize(h)
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
