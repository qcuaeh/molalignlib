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

module assorting
use kinds
use discrete
use sorting
use chemdata

implicit none

private
public getblocks

contains

subroutine getblocks(natom, znums, types, nblk, blksz, blkid, atomorder)
! Purpose: Group atoms by atomic numbers and types
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: znums, types
   integer, intent(out) :: nblk
   integer, intent(out) :: atomorder(:)
   integer, dimension(:), intent(out) :: blksz, blkid

   integer :: i, j
   logical :: remaining(natom)
   integer :: blkznum(natom)
   integer :: blktype(natom)
   integer :: blkorder(natom)

! Initialization

   nblk = 0
   remaining = .true.

! Create block list

   do i = 1, natom
      if (remaining(i)) then
         nblk = nblk + 1
         blksz(nblk) = 1
         blkznum(nblk) = znums(i)
         blktype(nblk) = types(i)
         blkid(i) = nblk
         do j = i + 1, natom
            if (remaining(i)) then
               if (znums(j) == znums(i) .and. types(j) == types(i)) then
                  blkid(j) = nblk
                  blksz(nblk) = blksz(nblk) + 1
                  remaining(j) = .false.
               end if
            end if
         end do
      end if
   end do

! Order blocks by atomic type

   blkorder(:nblk) = order(blktype, nblk)
   blksz(:nblk) = blksz(blkorder(:nblk))
   blkorder(:nblk) = inverseperm(blkorder(:nblk))
   blkid = blkorder(blkid)

! Order blocks by atomic number

   blkorder(:nblk) = order(blkznum, nblk)
   blksz(:nblk) = blksz(blkorder(:nblk))
   blkorder(:nblk) = inverseperm(blkorder(:nblk))
   blkid = blkorder(blkid)

! Order blocks by block size

!    blkorder(:nblk) = order(blksz, nblk)
!    blksz(:nblk) = blksz(blkorder(:nblk))
!    blkorder(:nblk) = inverseperm(blkorder(:nblk))
!    blkid = blkorder(blkid)

! Save atom order

   atomorder = order(blkid, natom)

end subroutine

end module
