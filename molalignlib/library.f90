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
use stdio
use kinds
use bounds

implicit none

contains

! Assign atoms0 and atoms1
subroutine assign_atoms( &
   znums0,  &
   types0, &
   coords0, &
   adjmat0, &
   weights0, &
   znums1, &
   types1, &
   coords1, &
   adjmat1, &
   weights1, &
   permlist, &
   countlist, &
   nrec, &
   error)

   use random
   use discrete
   use sorting
   use assorting
   use translation
   use adjacency
   use assignment

   integer, dimension(:), intent(in) :: znums0, types0
   integer, dimension(:), intent(in) :: znums1, types1
   logical, dimension(:, :), intent(in) :: adjmat0, adjmat1
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   real(wp), dimension(:), intent(in) :: weights0, weights1
   integer, dimension(:, :), intent(inout) :: permlist
   integer, dimension(:), intent(inout) :: countlist
   integer, intent(out) :: nrec, error
   integer :: i, h, offset
   integer :: nblk0, nblk1
   integer :: neqv0, neqv1
   integer, dimension(:), allocatable :: blkid0, blkid1
   integer, dimension(:), allocatable :: blksz0, blksz1
   integer, dimension(:), allocatable :: atomorder0, atomorder1
   integer, dimension(:), allocatable :: invatomorder0, invatomorder1
   integer, dimension(:), allocatable :: nadj0, nadj1
   integer, dimension(:, :), allocatable :: adjlist0, adjlist1
   integer, dimension(:), allocatable :: eqvid0, eqvid1
   integer, dimension(:), allocatable :: eqvsz0, eqvsz1
   real(wp) :: totalweight
   real(wp) :: center0(3), center1(3)

   ! Set error code to 0 by default

   error = 0

   ! Abort if molecules have different number of atoms

   if (natom0 /= natom1) then
      write (error_unit, '(a)') 'Error: The molecules have different number of atoms'
      error = 1
      return
   end if

   ! Allocate arrays

   allocate(atomorder0(natom0), atomorder1(natom1))
   allocate(invatomorder0(natom0), invatomorder1(natom1))
   allocate(blkid0(natom0), blkid1(natom1))
   allocate(blksz0(natom0), blksz1(natom1))
   allocate(nadj0(natom0), nadj1(natom1))
   allocate(adjlist0(maxcoord, natom0), adjlist1(maxcoord, natom1))
   allocate(eqvid0(natom0), eqvid1(natom1))
   allocate(eqvsz0(natom0), eqvsz1(natom1))

   ! Get adjacency lists

   call adjmat2list(natom0, adjmat0, nadj0, adjlist0)
   call adjmat2list(natom1, adjmat1, nadj1, adjlist1)

   ! Group atoms by label

   call grouptypes(natom0, znums0, types0, nblk0, blksz0, blkid0)
   call grouptypes(natom1, znums1, types1, nblk1, blksz1, blkid1)

   ! Group atoms by NMA at infinite level

   call groupequiv(natom0, nblk0, blkid0, nadj0, adjlist0, neqv0, eqvsz0, eqvid0)
   call groupequiv(natom1, nblk1, blkid1, nadj1, adjlist1, neqv1, eqvsz1, eqvid1)

   ! Get atom order

   atomorder0 = order(eqvid0, natom0)
   atomorder1 = order(eqvid1, natom1)

   ! Get inverse atom order

   invatomorder0 = inverseperm(atomorder0)
   invatomorder1 = inverseperm(atomorder1)

   ! Abort if molecules are not isomers

   if (any(znums0(atomorder0) /= znums1(atomorder1))) then
      write (error_unit, '(a)') 'Error: The molecules are not isomers'
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
   do h = 1, nblk0
      do i = 2, blksz0(h)
         if (weights0(atomorder0(offset+i)) /= weights0(atomorder0(offset+1))) then
            write (error_unit, '(a)') 'Error: Atoms of the same type must weight the same'
            error = 1
            return
         end if
      end do
      offset = offset + blksz0(h)
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
      natom0, nblk0, blksz0, &
      weights0(atomorder0)/totalweight, &
      centered(natom0, coords0(:, atomorder0), center0), &
      centered(natom1, coords1(:, atomorder1), center1), &
      maxrec, nrec, permlist, countlist)

   ! Reorder back to original atom ordering

   do i = 1, nrec
      permlist(:, i) = atomorder1(permlist(invatomorder0, i))
   end do

end subroutine

! Align atoms0 and atoms1
subroutine align_atoms( &
   znums0, &
   types0, &
   coords0, &
   weights0, &
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

   integer, dimension(:), intent(in) :: znums0, types0
   integer, dimension(:), intent(in) :: znums1, types1
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   real(wp), dimension(:), intent(in) :: weights0, weights1
   real(wp), intent(out) :: travec(3), rotmat(3, 3)
   integer, intent(out) :: error
   real(wp) :: totalweight
   real(wp) :: center0(3), center1(3)
   real(wp) :: aligned1(3, natom1)

   ! Set error code to 0 by default

   error = 0

   ! Abort if molecules have different number of atoms

   if (natom0 /= natom1) then
      write (error_unit, '(a)') 'Error: The molecules have different number of atoms'
      error = 1
      return
   end if

   ! Abort if molecules are not isomers

   if (any(sorted(znums0, natom0) /= sorted(znums1, natom1))) then
      write (error_unit, '(a)') 'Error: The molecules are not isomers'
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
