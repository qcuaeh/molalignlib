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
use random
use sorting
use discrete
use assorting
use rotation
use translation
use adjacency
use alignment
use assignment
use writemol
use biasing

implicit none

contains

! Assign atoms0 and atoms1
subroutine assign_atoms( &
   natom0, &
   znums0,  &
   types0, &
   weights0, &
   coords0, &
   adjmat0, &
   natom1, &
   znums1, &
   types1, &
   weights1, &
   coords1, &
   adjmat1, &
   permlist, &
   countlist, &
   nrec, &
   error)

   integer, intent(in) :: natom0, natom1
   integer, dimension(:), intent(inout) :: znums0, znums1
   integer, dimension(:), intent(inout) :: types0, types1
   logical, dimension(:, :), intent(inout) :: adjmat0, adjmat1
   real(wp), dimension(:, :), intent(inout) :: coords0, coords1
   real(wp), dimension(:), intent(inout) :: weights0, weights1
   integer, dimension(:, :), intent(inout) :: permlist
   integer, dimension(:), intent(inout) :: countlist
   integer, intent(out) :: nrec, error

   integer :: h, i, j, k, l
   integer :: nblk0, nblk1
   integer :: neqv0, neqv1
   integer, dimension(:), allocatable :: blkid0, blkid1
   integer, dimension(:), allocatable :: blksz0, blksz1
   integer, dimension(:), allocatable :: atomorder0, atomorder1
   integer, dimension(:), allocatable :: backorder0, backorder1
   integer, dimension(:), allocatable :: nadj0, nadj1
   integer, dimension(:, :), allocatable :: adjlist0, adjlist1
   integer, dimension(:), allocatable :: eqvid0, eqvid1
   integer, dimension(:), allocatable :: eqvsz0, eqvsz1
   real(wp), dimension(:), allocatable :: blkwt0, blkwt1
   real(wp) :: center0(3), center1(3)

   integer, dimension(natom0) :: offset, blkid, invatomperm
   integer, dimension(natom0, natom0) :: order01
   integer :: nbond0, bonds0(2, maxcoord*natom0)
   integer :: nbond1, bonds1(2, maxcoord*natom1)
   logical, dimension(natom0, natom0) :: adjmat00, adjmat01
   real(wp), dimension(natom0, natom0) :: d01

   ! Set error code to 0 by default

   error = 0

   !  Select bias function

   if (bias_flag) then
      if (bond_flag) then
         bias_func => mnacrossbias
      else
         bias_func => sndcrossbias
      end if
   else
      bias_func => nocrossbias
   end if

   ! Abort if molecules have different number of atoms

   if (natom0 /= natom1) then
      write (error_unit, '(a)') 'Error: The molecules are not isomers'
      error = 1
      return
   end if

   ! Allocate arrays

   allocate(atomorder0(natom0), atomorder1(natom1))
   allocate(backorder0(natom0), backorder1(natom1))
   allocate(blkid0(natom0), blkid1(natom1))
   allocate(blksz0(natom0), blksz1(natom1))
   allocate(blkwt0(natom0), blkwt1(natom1))
   allocate(nadj0(natom0), nadj1(natom1))
   allocate(adjlist0(maxcoord, natom0), adjlist1(maxcoord, natom1))
   allocate(eqvid0(natom0), eqvid1(natom1))
   allocate(eqvsz0(natom0), eqvsz1(natom1))

   ! Calculate adjacency lists

   call adjmat2list(natom0, adjmat0, nadj0, adjlist0)
   call adjmat2list(natom1, adjmat1, nadj1, adjlist1)

   ! Group atoms by label

   call groupblocks(natom0, znums0, types0, weights0, nblk0, blksz0, blkwt0, blkid0)
   call groupblocks(natom1, znums1, types1, weights1, nblk1, blksz1, blkwt1, blkid1)

   ! Group atoms by NMA at infinite level

   call groupequiv(natom0, nblk0, blkid0, nadj0, adjlist0, neqv0, eqvsz0, eqvid0)
   call groupequiv(natom1, nblk1, blkid1, nadj1, adjlist1, neqv1, eqvsz1, eqvid1)

   ! Get atom order

   atomorder0 = sortorder(eqvid0, natom0)
   atomorder1 = sortorder(eqvid1, natom1)

   ! Get inverse atom order

   backorder0 = inverseperm(atomorder0)
   backorder1 = inverseperm(atomorder1)

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

   if (any(abs(weights0(atomorder0) - weights1(atomorder1)) > 1.E-6)) then
      write (error_unit, '(a)') 'Error: There are conflicting weights'
      error = 1
      return
   end if

   ! Calculate centroids

   center0 = centroid(natom0, weights0, coords0)
   center1 = centroid(natom1, weights1, coords1)

   ! Initialize random number generator

   call initialize_random()

   ! Remap atoms to minimize distance and difference

!   call optimize_assignment( &
!      natom0, &
!      nblk0, &
!      blksz0, &
!      blkwt0, &
!      neqv0, &
!      eqvsz0, &
!      centered(natom0, coords0(:, atomorder0), center0), &
!      adjmat0(atomorder0, atomorder0), &
!      neqv1, &
!      eqvsz1, &
!      centered(natom1, coords1(:, atomorder1), center1), &
!      adjmat1(atomorder1, atomorder1), &
!      permlist, &
!      countlist, &
!      nrec)

   znums0 = znums0(atomorder0)
   znums1 = znums1(atomorder1)
   types0 = types0(atomorder0)
   types1 = types1(atomorder1)
   weights0 = weights0(atomorder0)
   weights1 = weights1(atomorder1)
   coords0 = coords0(:, atomorder0)
   coords1 = coords1(:, atomorder1)
   adjmat0 = adjmat0(atomorder0, atomorder0)
   adjmat1 = adjmat1(atomorder1, atomorder1)
   call translate(natom0, coords0, -center0)
   call translate(natom1, coords1, -center1)

   call optimize_assignment( &
      natom0, &
      nblk0, &
      blksz0, &
      blkwt0, &
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

   call rotate(natom1, coords1, &
      rotquat2rotmat(leastrotquat(natom0, weights0, coords0, coords1, permlist(:, 1))) &
   )

   offset(1) = 0
   do h = 1, nblk0 - 1
      offset(h+1) = offset(h) + blksz0(h)
   end do

   do h = 1, nblk0
      blkid(offset(h)+1:offset(h)+blksz0(h)) = h
   end do

   invatomperm = inverseperm(permlist(:, 1))

   adjmat00 = adjmat0
   adjmat01 = adjmat1

   do i = 1, natom0
      do j = i + 1, natom0
         if (adjmat00(i, j) .neqv. adjmat01(permlist(i, 1), permlist(j, 1))) then
            adjmat0(i, j) = .false.
            adjmat0(j, i) = .false.
            adjmat1(permlist(i, 1), permlist(j, 1)) = .false.
            adjmat1(permlist(j, 1), permlist(i, 1)) = .false.
            h = blkid(i)
            do k = offset(h) + 1, offset(h) + blksz0(h)
               if (sum((coords0(:, i) - coords1(:, k))**2) &
                  < sum((coords0(:, i) - coords1(:, permlist(i, 1)))**2) &
               ) then
!                  print *, i, permlist(i, 1), invatomperm(k)
                  adjmat0(invatomperm(k), j) = .false.
                  adjmat0(j, invatomperm(k)) = .false.
                  adjmat1(k, j) = .false.
                  adjmat1(j, k) = .false.
!                  adjmat0(invatomperm(k), :) = .false.
!                  adjmat0(:, invatomperm(k)) = .false.
!                  adjmat1(k, :) = .false.
!                  adjmat1(:, k) = .false.
               end if
            end do
         end if
      end do
   end do

   call optimize_assignment( &
      natom0, &
      nblk0, &
      blksz0, &
      blkwt0, &
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

   ! Print coordinates with internal order

   open(unit=99, file='ordered.mol2', action='write', status='replace')
   call adjmat2bonds(natom0, adjmat0, nbond0, bonds0)
!   call adjmat2bonds(natom1, adjmat1, nbond1, bonds1)
   call adjmat2bonds(natom1, adjmat1(permlist(:, 1), permlist(:, 1)), nbond1, bonds1)
   call writemol2(99, 'coords0', natom0, znums0, coords0, nbond0, bonds0)
!   call writemol2(99, 'coords1', natom1, znums1, coords1, nbond1, bonds1)
   call writemol2(99, 'coords1', natom1, znums1(permlist(:, 1)), coords1(:, permlist(:, 1)), nbond1, bonds1)

   ! Reorder back to original atom ordering

!   do i = 1, nrec
!      permlist(:, i) = atomorder1(permlist(backorder0, i))
!   end do

end subroutine

! Align atoms0 and atoms1
subroutine align_atoms( &
   natom0, &
   znums0, &
   types0, &
   weights0, &
   coords0, &
   natom1, &
   znums1, &
   types1, &
   weights1, &
   coords1, &
   travec, &
   rotmat, &
   error)

   integer, intent(in) :: natom0, natom1
   integer, dimension(:), intent(in) :: znums0, types0
   integer, dimension(:), intent(in) :: znums1, types1
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   real(wp), dimension(:), intent(in) :: weights0, weights1
   real(wp), intent(out) :: travec(3), rotmat(3, 3)
   integer, intent(out) :: error
   real(wp) :: center0(3), center1(3)
   real(wp) :: aligned1(3, natom1)

   ! Set error code to 0 by default

   error = 0

   ! Abort if molecules have different number of atoms

   if (natom0 /= natom1) then
      write (error_unit, '(a)') 'Error: The molecules are not isomers'
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

   ! Calculate centroids

   center0 = centroid(natom0, weights0, coords0)
   center1 = centroid(natom1, weights1, coords1)

   ! Calculate optimal rotation matrix

   rotmat = rotquat2rotmat(leastrotquat( &
      natom0, &
      weights0, &
      centered(natom0, coords0, center0), &
      centered(natom1, coords1, center1), &
      identityperm(natom0)))

   ! Calculate optimal translation vector

   travec = center0 - matmul(rotmat, center1)

   ! Calculate RMSD

   aligned1 = translated(natom1, rotated(natom1, coords1, rotmat), travec)

end subroutine

end module
