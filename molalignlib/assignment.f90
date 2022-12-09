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
   mapind, &
   mapcount, &
   mapdist &
)

   integer, intent(in) :: natom, nblock, nrec
   integer, dimension(:), intent(in) :: blocksize
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   real(wp), dimension(:), intent(in) :: weights
   integer, intent(out) :: nmap
   integer, intent(out) :: mapind(:, :)
   integer, intent(out) :: mapcount(:)
   real(wp), intent(out) :: mapdist(:)

   logical found, overflow
   integer imap, jmap, ntrial, nmatch, cycles
   integer, dimension(natom) :: atomap, auxmap
   real(wp) :: dist, biased_dist, new_biased_dist, totalrot
   real(wp), dimension(4) :: rotquat, prodquat
   real(wp), dimension(nrec) :: avgiter, avgtotalrot, avgrealrot
   real(wp) :: bias(natom, natom)
   real(wp) :: auxcoords1(3, natom)

! Set bias for non equivalent atoms 

   call setadjbias(natom, nblock, blocksize, coords0, coords1, bias)

! Print header and initial stats

   if (live_flag) then
      write (output_unit, '(a)', advance='no') achar(27)//'[1H'//achar(27)//'[J'
      call print_header()
   end if

! Initialize loop variables

   nmap = 0
   ntrial = 0
   nmatch = 0
   overflow = .false.

! Loop for map searching

   do while (nmatch < max_count .and. (free_flag .or. ntrial < max_trials))

      ntrial = ntrial + 1

! Work with a copy of coords1

      auxcoords1 = coords1

! Aply a random rotation to auxcoords1

      call rotate(natom, auxcoords1, torotquat(rand3()))

! Minimize the euclidean distance

      call getatomap(natom, coords0, auxcoords1, nblock, blocksize, bias, atomap)
      rotquat = leastrotquat(natom, weights, coords0, auxcoords1, atomap)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate(natom, auxcoords1, rotquat)
      cycles = 1

      do while (iter_flag)
         biased_dist = squadist(natom, weights, coords0, auxcoords1, atomap) &
               + totalbias(natom, weights, bias, atomap)
         call getatomap(natom, coords0, auxcoords1, nblock, blocksize, bias, auxmap)
         if (all(auxmap == atomap)) exit
         new_biased_dist = squadist(natom, weights, coords0, auxcoords1, auxmap) &
               + totalbias(natom, weights, bias, auxmap)
         if (new_biased_dist > biased_dist) then
            write (output_unit, '(a)') 'new_biased_dist is larger than biased_dist!'
!                print *, biased_dist, new_biased_dist
         end if
         rotquat = leastrotquat(natom, weights, coords0, auxcoords1, auxmap)
         prodquat = quatmul(rotquat, prodquat)
         call rotate(natom, auxcoords1, rotquat)
         cycles = cycles + 1
         totalrot = totalrot + rotangle(rotquat)
         atomap = auxmap
      end do

      dist = squadist(natom, weights, coords0, auxcoords1, atomap)

! Check for new best mapping

      found = .false.

      do imap = 1, nmap
         if (all(atomap == mapind(:, imap))) then
            if (imap == 1) nmatch = nmatch + 1
            mapcount(imap) = mapcount(imap) + 1
            avgiter(imap) = avgiter(imap) + (cycles - avgiter(imap))/mapcount(imap)
            avgrealrot(imap) = avgrealrot(imap) + (rotangle(prodquat) - avgrealrot(imap))/mapcount(imap)
            avgtotalrot(imap) = avgtotalrot(imap) + (totalrot - avgtotalrot(imap))/mapcount(imap)
            if (live_flag) then
               write (output_unit, '(a)', advance='no') achar(27)//'['//str(imap + 2)//'H'
               call print_stats(imap, mapcount(imap), avgiter(imap), avgtotalrot(imap), &
                  avgrealrot(imap), mapdist(imap))
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
            if (imap > nmap .or. dist < mapdist(imap)) then
               if (imap == 1) nmatch = 1
               if (nmap < nrec) nmap = nmap + 1
               do jmap = nmap, imap + 1, -1
                  mapind(:, jmap) = mapind(:, jmap - 1)
                  mapcount(jmap) = mapcount(jmap - 1)
                  mapdist(jmap) = mapdist(jmap - 1)
                  avgiter(jmap) = avgiter(jmap - 1)
                  avgrealrot(jmap) = avgrealrot(jmap - 1)
                  avgtotalrot(jmap) = avgtotalrot(jmap - 1)
                  if (live_flag) then
                     write (output_unit, '(a)', advance='no') achar(27)//'['//str(jmap + 2)//'H'
                     call print_stats(jmap, mapcount(jmap), avgiter(jmap), avgtotalrot(jmap), &
                        avgrealrot(jmap), mapdist(jmap))
                  end if
               end do
               mapind(:, imap) = atomap
               mapcount(imap) = 1
               mapdist(imap) = dist
               avgiter(imap) = cycles
               avgrealrot(imap) = rotangle(prodquat)
               avgtotalrot(imap) = totalrot
               if (live_flag) then
                  write (output_unit, '(a)', advance='no') achar(27)//'['//str(imap + 2)//'H'
                  call print_stats(imap, mapcount(imap), avgiter(imap), avgtotalrot(imap), &
                     avgrealrot(imap), mapdist(imap))
               end if
               exit
            end if
         end do
      end if

      if (live_flag) then
         write (output_unit, '(a)', advance='no') achar(27)//'['//str(nmap + 3)//'H'
         call print_footer(overflow, nmap, ntrial)
      end if

   end do

   if (.not. live_flag) then
      call print_header()
      do imap = 1, nmap
         call print_stats(imap, mapcount(imap), avgiter(imap), avgtotalrot(imap), &
            avgrealrot(imap), mapdist(imap))
      end do
      call print_footer(overflow, nmap, ntrial)
   end if

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine getatomap(natom, coords0, coords1, nblock, blocksize, bias, atomap)

! nblock: Number of block atoms
! blocksize: Number of atoms in each block
! atomap: Map between correspondent points in the adjmat
! offset: First element of current block

   integer, intent(in) :: natom, nblock
   integer, dimension(:), intent(in) :: blocksize
   real(wp), dimension(:, :), intent(in) :: coords0
   real(wp), dimension(:, :), intent(in) :: bias
   real(wp), dimension(:, :), intent(inout) :: coords1
   integer, dimension(:), intent(out) :: atomap

   integer :: h, offset
   integer, dimension(natom) :: perm
   real(wp) :: dist

! Fill distance matrix for each block

   offset = 0

   do h = 1, nblock
      call minperm(blocksize(h), offset, coords0, coords1, bias, perm, dist)
      atomap(offset+1:offset+blocksize(h)) = perm(:blocksize(h)) + offset
      offset = offset + blocksize(h)
   end do

end subroutine

! Calculate total bias
real(wp) function totalbias(natom, weights, bias, map) result(total)

   integer, intent(in) :: natom
   real(wp), dimension(:), intent(in) :: weights
   integer, dimension(:), intent(in) :: map
   real(wp), dimension(:, :), intent(in) :: bias
   integer :: i

   total = 0.

   do i = 1, natom
      total = total + weights(i)*bias(i, map(i))
   end do

end function

end module
