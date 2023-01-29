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

module library
use parameters
use settings

implicit none

contains

! Assign atoms0 and atoms1
subroutine assign_atoms( &
   natom0, &
   znums0,  &
   types0, &
   coords0, &
   weights0, &
   natom1, &
   znums1, &
   types1, &
   coords1, &
   weights1, &
   maxrec, &
   nrec, &
   permlist, &
   countlist, &
   error)

   use random
   use discrete
   use sorting
   use assorting
   use translation
   use assignment

   integer, intent(in) :: natom0, natom1, maxrec
   integer, dimension(natom0), intent(in) :: znums0, types0
   integer, dimension(natom1), intent(in) :: znums1, types1
   real(wp), intent(in) :: coords0(3, natom0)
   real(wp), intent(in) :: coords1(3, natom1)
   real(wp), intent(in) :: weights0(natom0)
   real(wp), intent(in) :: weights1(natom1)
   integer, intent(out) :: nrec, error
   integer, intent(out) :: permlist(natom0, maxrec)
   integer, intent(out) :: countlist(maxrec)
   integer :: i, h, offset
   integer :: nblock0, nblock1
   integer, dimension(:), allocatable :: atomorder0, atomorder1
   integer, dimension(:), allocatable :: atomunorder0, atomunorder1
   integer, dimension(:), allocatable :: blockidx0, blockidx1
   integer, dimension(:), allocatable :: blocksize0, blocksize1
   real(wp) :: center0(3), center1(3), totalweight

   ! Set error code to 0 by default

   error = 0

   ! Abort if clusters have different number of atoms

   if (natom0 /= natom1) then
      write (error_unit, '(a)') 'Error: The clusters have different number of atoms'
      error = 1
      return
   end if

   ! Allocate arrays

   allocate(atomorder0(natom0), atomorder1(natom1))
   allocate(atomunorder0(natom0), atomunorder1(natom1))
   allocate(blockidx0(natom0), blockidx1(natom1))
   allocate(blocksize0(natom0), blocksize1(natom1))

   ! Group atoms by label

   call getblocks(natom0, znums0, types0, nblock0, blocksize0, blockidx0, atomorder0)
   call getblocks(natom1, znums1, types1, nblock1, blocksize1, blockidx1, atomorder1)

   ! Get inverse atom ordering

   atomunorder0 = inverseperm(atomorder0)
   atomunorder1 = inverseperm(atomorder1)

   ! Abort if clusters are not isomers

   if (any(znums0(atomorder0) /= znums1(atomorder1))) then
      write (error_unit, '(a)') 'Error: The clusters are not isomers'
      error = 1
      return
   end if

   ! Abort if there are conflicting types

   if (any(types0(atomorder0) /= types1(atomorder1))) then
      write (error_unit, '(a)') 'Error: There are conflicting atom types'
      error = 1
      return
   end if

   ! Abort if there are conflicting weights

   if (any(weights0(atomorder0) /= weights1(atomorder1))) then
      write (error_unit, '(a)') 'Error: There are conflicting weights'
      error = 1
      return
   end if

   ! Abort if there are inconsistent weights

   offset = 0
   do h = 1, nblock0
      do i = 2, blocksize0(h)
         if (weights0(atomorder0(offset+i)) /= weights0(atomorder0(offset+1))) then
            write (error_unit, '(a)') 'Error: Atoms of the same type must weight the same'
            error = 1
            return
         end if
      end do
      offset = offset + blocksize0(h)
   end do

   ! Calculate total weight

   totalweight = sum(weights0)

   ! Calculate centroids

   center0 = centroid(natom0, weights0, coords0)
   center1 = centroid(natom1, weights1, coords1)

   ! Initialize random number generator

   call initialize_random()

   ! Remap atoms to minimize distance and difference

   call optimize_assignment( &
      natom0, nblock0, blocksize0, &
      weights0(atomorder0)/totalweight, &
      centered(natom0, coords0(:, atomorder0), center0), &
      centered(natom1, coords1(:, atomorder1), center1), &
      maxrec, nrec, permlist, countlist)

   ! Reorder back to original atom ordering

   do i = 1, nrec
      permlist(:, i) = atomorder1(permlist(atomunorder0, i))
   end do

end subroutine

! Align atoms0 and atoms1
subroutine align_atoms( &
   natom0, &
   znums0, &
   types0, &
   coords0, &
   weights0, &
   natom1, &
   znums1, &
   types1, &
   coords1, &
   weights1, &
   travec, &
   rotmat, &
   error)

   use discrete
   use sorting
   use translation
   use rotation
   use alignment

   integer, intent(in) :: natom0, natom1
   integer, dimension(natom0), intent(in) :: znums0, types0
   integer, dimension(natom1), intent(in) :: znums1, types1
   real(wp), intent(in) :: weights0(natom0)
   real(wp), intent(in) :: weights1(natom1)
   real(wp), intent(in) :: coords0(3, natom0)
   real(wp), intent(in) :: coords1(3, natom1)
   real(wp), intent(out) :: travec(3), rotmat(3, 3)
   integer, intent(out) :: error
   real(wp) :: totalweight
   real(wp) :: center0(3), center1(3)
   real(wp) :: aligned1(3, natom1)

   ! Set error code to 0 by default

   error = 0

   ! Abort if clusters have different number of atoms

   if (natom0 /= natom1) then
      write (error_unit, '(a)') 'Error: The clusters have different number of atoms'
      error = 1
      return
   end if

   ! Abort if clusters are not isomers

   if (any(sorted(znums0, natom0) /= sorted(znums1, natom1))) then
      write (error_unit, '(a)') 'Error: The clusters are not isomers'
      error = 1
      return
   end if

   ! Abort if atoms are not ordered

   if (any(znums0 /= znums1)) then
      write (error_unit, '(a)') 'Error: The atoms are not in the same order'
      error = 1
      return
   end if

   ! Abort if there are conflicting types

   if (any(sorted(types0, natom0) /= sorted(types1, natom1))) then
      write (error_unit, '(a)') 'Error: There are conflicting atom types'
      error = 1
      return
   end if

   ! Abort if types are not ordered

   if (any(types0 /= types1)) then
      write (error_unit, '(a)') 'Error: The atom types are not in the same order'
      error = 1
      return
   end if

   ! Calculate total weight

   totalweight = sum(weights0)

   ! Calculate centroids

   center0 = centroid(natom0, weights0, coords0)
   center1 = centroid(natom1, weights1, coords1)

   ! Calculate optimal rotation matrix

   rotmat = rotquat2rotmat(leastrotquat( &
      natom0, weights0/totalweight, &
      centered(natom0, coords0, center0), &
      centered(natom1, coords1, center1), &
      identityperm(natom0)))

   ! Calculate optimal translation vector

   travec = center0 - matmul(rotmat, center1)

   ! Calculate RMSD

   aligned1 = translated(natom1, rotated(natom1, coords1, rotmat), travec)

end subroutine

end module
