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
integer :: nreac
integer, dimension(100) :: reacatom 

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
   maplist, &
   countlist, &
   nrec, &
   error)

   integer, intent(in) :: natom0, natom1
   integer, dimension(:), intent(inout) :: znums0, znums1
   integer, dimension(:), intent(inout) :: types0, types1
   logical, dimension(:, :), intent(inout) :: adjmat0, adjmat1
   real(wp), dimension(:, :), intent(inout) :: coords0, coords1
   real(wp), dimension(:), intent(inout) :: weights0, weights1
   integer, dimension(:, :), intent(inout) :: maplist
   integer, dimension(:), intent(inout) :: countlist
   integer, intent(out) :: nrec, error

   integer :: h, i, j, k, l
   integer :: nblk0, nblk1
   integer :: neqv0, neqv1
   integer, dimension(:), allocatable :: blkidx0, blkidx1
   integer, dimension(:), allocatable :: blklen0, blklen1
   integer, dimension(:), allocatable :: atomorder0, atomorder1
   integer, dimension(:), allocatable :: backorder0, backorder1
   integer, dimension(:), allocatable :: nadj0, nadj1
   integer, dimension(:, :), allocatable :: adjlist0, adjlist1
   integer, dimension(:), allocatable :: eqvidx0, eqvidx1
   integer, dimension(:), allocatable :: eqvlen0, eqvlen1
   real(wp), dimension(:), allocatable :: blkwgt0, blkwgt1
   real(wp) :: travec0(3), travec1(3)

   integer, dimension(natom0) :: offset, blkidx, mapping
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
   allocate(blkidx0(natom0), blkidx1(natom1))
   allocate(blklen0(natom0), blklen1(natom1))
   allocate(blkwgt0(natom0), blkwgt1(natom1))
   allocate(nadj0(natom0), nadj1(natom1))
   allocate(adjlist0(maxcoord, natom0), adjlist1(maxcoord, natom1))
   allocate(eqvidx0(natom0), eqvidx1(natom1))
   allocate(eqvlen0(natom0), eqvlen1(natom1))

   ! Calculate adjacency lists

   call adjmat2list(natom0, adjmat0, nadj0, adjlist0)
   call adjmat2list(natom1, adjmat1, nadj1, adjlist1)

   ! Group atoms by label

   call groupblocks(natom0, znums0, types0, weights0, nblk0, blklen0, blkwgt0, blkidx0)
   call groupblocks(natom1, znums1, types1, weights1, nblk1, blklen1, blkwgt1, blkidx1)

   ! Group atoms by NMA at infinite level

   call groupequiv(natom0, nblk0, blkidx0, nadj0, adjlist0, neqv0, eqvlen0, eqvidx0)
   call groupequiv(natom1, nblk1, blkidx1, nadj1, adjlist1, neqv1, eqvlen1, eqvidx1)

   ! Get atom order

   atomorder0 = sortorder(eqvidx0, natom0)
   atomorder1 = sortorder(eqvidx1, natom1)

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

   travec0 = -centroid(natom0, weights0, coords0)
   travec1 = -centroid(natom1, weights1, coords1)

   ! Initialize random number generator

   call initialize_random()

   ! Remap atoms to minimize distance and difference

   call optimize_assignment( &
      natom0, &
      nblk0, &
      blklen0, &
      blkwgt0, &
      neqv0, &
      eqvlen0, &
      translated(natom0, coords0(:, atomorder0), travec0), &
      adjmat0(atomorder0, atomorder0), &
      neqv1, &
      eqvlen1, &
      translated(natom1, coords1(:, atomorder1), travec1), &
      adjmat1(atomorder1, atomorder1), &
      maplist, &
      countlist, &
      nrec)

   ! Reorder back to original atom ordering

   do i = 1, nrec
      maplist(:, i) = atomorder1(maplist(backorder0, i))
   end do

!   znums0 = znums0(atomorder0)
!   znums1 = znums1(atomorder1)
!   types0 = types0(atomorder0)
!   types1 = types1(atomorder1)
!   weights0 = weights0(atomorder0)
!   weights1 = weights1(atomorder1)
!   coords0 = coords0(:, atomorder0)
!   coords1 = coords1(:, atomorder1)
!   adjmat0 = adjmat0(atomorder0, atomorder0)
!   adjmat1 = adjmat1(atomorder1, atomorder1)
!
!   ! Mirror coordinates
!
!   if (mirror_flag) then
!      coords1(1, :) = -coords1(1, :)
!   end if
!
!   ! Translate to center of mass
!
!   call translate(natom0, coords0, travec0)
!   call translate(natom1, coords1, travec1)
!
!   call optimize_assignment( &
!      natom0, &
!      nblk0, &
!      blklen0, &
!      blkwgt0, &
!      neqv0, &
!      eqvlen0, &
!      coords0, &
!      adjmat0, &
!      neqv1, &
!      eqvlen1, &
!      coords1, &
!      adjmat1, &
!      maplist, &
!      countlist, &
!      nrec)
!
!   mapping = maplist(:, 1)
!
!   call rotate(natom1, coords1, leastrotquat(natom0, weights0, coords0, coords1, mapping))
!
!   offset(1) = 0
!   do h = 1, nblk0 - 1
!      offset(h+1) = offset(h) + blklen0(h)
!   end do
!
!   do h = 1, nblk0
!      blkidx(offset(h)+1:offset(h)+blklen0(h)) = h
!   end do
!
!   adjmat00 = adjmat0
!   adjmat01 = adjmat1
!
!   do i = 1, natom0
!      do j = 1, natom0
!         if (adjmat00(i, j) .neqv. adjmat01(mapping(i), mapping(j))) then
!            adjmat0(i, j) = .false.
!            adjmat1(mapping(i), mapping(j)) = .false.
!            do l = 1, nreac
!               adjmat0(i, reacatom(l)) = .false.
!               adjmat0(reacatom(l), i) = .false.
!               adjmat1(mapping(i), mapping(reacatom(l))) = .false.
!               adjmat1(mapping(reacatom(l)), mapping(i)) = .false.
!            end  do
!            h = blkidx(i)
!            do k = offset(h) + 1, offset(h) + blklen0(h)
!               if (sum((coords0(:, i) - coords1(:, mapping(k)))**2) < 2.0 &
!                  .or. sum((coords0(:, k) - coords1(:, mapping(i)))**2) < 2.0 &
!               ) then
!!                  print *, '<', i, mapping(i), k, mapping(k)
!                  adjmat0(k, j) = .false.
!                  adjmat0(j, k) = .false.
!                  adjmat1(mapping(k), mapping(j)) = .false.
!                  adjmat1(mapping(j), mapping(k)) = .false.
!                  do l = 1, nreac
!                     adjmat0(k, reacatom(l)) = .false.
!                     adjmat0(reacatom(l), k) = .false.
!                     adjmat1(mapping(k), mapping(reacatom(l))) = .false.
!                     adjmat1(mapping(reacatom(l)), mapping(k)) = .false.
!                  end  do
!               end if
!            end do
!         end if
!      end do
!   end do
!
!   call optimize_assignment( &
!      natom0, &
!      nblk0, &
!      blklen0, &
!      blkwgt0, &
!      neqv0, &
!      eqvlen0, &
!      coords0, &
!      adjmat0, &
!      neqv1, &
!      eqvlen1, &
!      coords1, &
!      adjmat1, &
!      maplist, &
!      countlist, &
!      nrec)
!
!   ! Print coordinates with internal order
!
!   open(unit=99, file='ordered.mol2', action='write', status='replace')
!   call adjmat2bonds(natom0, adjmat0, nbond0, bonds0)
!!   call adjmat2bonds(natom1, adjmat1, nbond1, bonds1)
!   call adjmat2bonds(natom1, adjmat1(maplist(:, 1), maplist(:, 1)), nbond1, bonds1)
!   call writemol2(99, 'coords0', natom0, znums0, coords0, nbond0, bonds0)
!!   call writemol2(99, 'coords1', natom1, znums1, coords1, nbond1, bonds1)
!   call writemol2(99, 'coords1', natom1, znums1(maplist(:, 1)), coords1(:, maplist(:, 1)), nbond1, bonds1)

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
   real(wp) :: travec0(3), travec1(3), rotquat(4)
   integer, intent(out) :: error

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

   travec0 = -centroid(natom0, weights0, coords0)
   travec1 = -centroid(natom1, weights1, coords1)

   ! Calculate optimal rotation matrix

   rotquat = leastrotquat( &
      natom0, &
      weights0, &
      translated(natom0, coords0, travec0), &
      translated(natom1, coords1, travec1), &
      identityperm(natom0) &
   )

   rotmat = rotquat2rotmat(rotquat)

   ! Calculate optimal translation vector

   travec = matmul(rotmat, travec1) - travec0

end subroutine

end module
