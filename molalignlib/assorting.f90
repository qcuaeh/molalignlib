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
use parameters
use settings
use discrete
use sorting
use chemdata

implicit none

private
public getblocks

contains

subroutine getblocks(natom, znums, types, nblock, bsize, blockindex, atomorder)
! Purpose: Group atoms by atomic numbers and types
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: znums, types
   integer, intent(out) :: nblock
   integer, intent(out) :: atomorder(:)
   integer, dimension(:), intent(out) :: bsize, blockindex

   integer :: i, j
   logical :: remaining(natom)
   integer :: blockznum(natom)
   integer :: blocktype(natom)
   integer :: blockorder(natom)

! Initialization

   nblock = 0
   remaining = .true.

! Create block list

   do i = 1, natom
      if (remaining(i)) then
         nblock = nblock + 1
         bsize(nblock) = 1
         blockznum(nblock) = znums(i)
         blocktype(nblock) = types(i)
         blockindex(i) = nblock
         do j = i + 1, natom
            if (remaining(i)) then
               if (znums(j) == znums(i) .and. types(j) == types(i)) then
                  blockindex(j) = nblock
                  bsize(nblock) = bsize(nblock) + 1
                  remaining(j) = .false.
               end if
            end if
         end do
      end if
   end do

! Order blocks by atomic type

   blockorder(:nblock) = order(blocktype, nblock)
   bsize(:nblock) = bsize(blockorder(:nblock))
   blockorder(:nblock) = inverseperm(blockorder(:nblock))
   blockindex = blockorder(blockindex)

! Order blocks by atomic number

   blockorder(:nblock) = order(blockznum, nblock)
   bsize(:nblock) = bsize(blockorder(:nblock))
   blockorder(:nblock) = inverseperm(blockorder(:nblock))
   blockindex = blockorder(blockindex)

! Order blocks by block size

!    blockorder(:nblock) = order(bsize, nblock)
!    bsize(:nblock) = bsize(blockorder(:nblock))
!    blockorder(:nblock) = inverseperm(blockorder(:nblock))
!    blockindex = blockorder(blockindex)

! Save atom order

   atomorder = order(blockindex, natom)

end subroutine

end module
