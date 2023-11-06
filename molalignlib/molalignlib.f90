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

module molalignlib
use stdio
use kinds
use bounds
use random
use sorting
use discrete
use moltypes
use rotation
use translation
use assorting
use adjacency
use alignment
use remapping
use assignment
use writemol
use biasing

implicit none

contains

! Assign atoms0 and atoms1
subroutine assign_atoms( &
   natom0, &
   znums0, &
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
   integer, dimension(:), intent(in) :: znums0, znums1
   integer, dimension(:), intent(in) :: types0, types1
   logical, dimension(:, :), intent(in) :: adjmat0, adjmat1
   real(wp), dimension(:), intent(in) :: weights0, weights1
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   integer, dimension(:, :), intent(inout) :: maplist
   integer, dimension(:), intent(inout) :: countlist
   integer, intent(out) :: nrec, error

   integer :: workznums0(natom0), workznums1(natom1)
   integer :: worktypes0(natom0), worktypes1(natom1)
   logical :: workadjmat0(natom0, natom0), workadjmat1(natom1, natom1)
   logical :: backadjmat0(natom0, natom0), backadjmat1(natom1, natom1)
   real(wp) :: workweights0(natom0), workweights1(natom1)
   real(wp) :: workcoords0(3, natom0), workcoords1(3, natom1)

   integer :: h, i, j, k
   integer :: nblk0, nblk1, neqv0, neqv1
   integer, dimension(:), allocatable :: nadj0, nadj1
   integer, dimension(:, :), allocatable :: adjlist0, adjlist1
   integer, dimension(:), allocatable :: blkidx0, blkidx1
   integer, dimension(:), allocatable :: blklen0, blklen1
   integer, dimension(:), allocatable :: atomorder0, atomorder1
   integer, dimension(:), allocatable :: backorder0, backorder1
   integer, dimension(:), allocatable :: eqvidx0, eqvidx1
   integer, dimension(:), allocatable :: eqvlen0, eqvlen1
   real(wp) :: travec0(3), travec1(3)

   integer, dimension(natom0) :: offset, blkidx, mapping
   integer :: nbond0, bonds0(2, maxcoord*natom0)
   integer :: nbond1, bonds1(2, maxcoord*natom1)

   ! Set error code to 0 by default

   error = 0

   !  Select assignment algorithm

   if (bond_flag) then
      mapatoms => mapatoms_bonded
      mapsetcrossbias => setcrossbias_mna
   else
      mapatoms => mapatoms_free
      mapsetcrossbias => setcrossbias_rd
   end if

   ! Select bias function

   if (bias_flag) then
      setcrossbias => mapsetcrossbias
   else
      setcrossbias => setcrossbias_none
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
   allocate(nadj0(natom0), nadj1(natom1))
   allocate(adjlist0(maxcoord, natom0), adjlist1(maxcoord, natom1))
   allocate(eqvidx0(natom0), eqvidx1(natom1))
   allocate(eqvlen0(natom0), eqvlen1(natom1))

   ! Calculate adjacency lists

   call adjmat2list(natom0, adjmat0, nadj0, adjlist0)
   call adjmat2list(natom1, adjmat1, nadj1, adjlist1)

   ! Group atoms by type

   call groupatoms(natom0, znums0, types0, weights0, nblk0, blklen0, blkidx0)
   call groupatoms(natom1, znums1, types1, weights1, nblk1, blklen1, blkidx1)

   ! Group atoms by MNA

   call groupequivatoms(natom0, nblk0, blkidx0, nadj0, adjlist0, neqv0, eqvlen0, eqvidx0)
   call groupequivatoms(natom1, nblk1, blkidx1, nadj1, adjlist1, neqv1, eqvlen1, eqvidx1)

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

   ! Reorder data arrays

   workznums0 = znums0(atomorder0)
   workznums1 = znums1(atomorder1)
   worktypes0 = types0(atomorder0)
   worktypes1 = types1(atomorder1)
   workweights0 = weights0(atomorder0)
   workweights1 = weights1(atomorder1)
   workcoords0 = coords0(:, atomorder0)
   workcoords1 = coords1(:, atomorder1)
   workadjmat0 = adjmat0(atomorder0, atomorder0)
   workadjmat1 = adjmat1(atomorder1, atomorder1)

   ! Mirror coordinates

   if (mirror_flag) then
      workcoords1(1, :) = -workcoords1(1, :)
   end if

   ! Calculate centroids

   travec0 = -centroid(natom0, workweights0, workcoords0)
   travec1 = -centroid(natom1, workweights1, workcoords1)

   ! Center coordinates at the centroids

   call translate(natom0, workcoords0, travec0)
   call translate(natom1, workcoords1, travec1)

   ! Initialize random number generator

   call initialize_random()

   ! Optimize the assignment to minimize AdjD/RMSD

   call optimize_assignment( &
      natom0, &
      nblk0, &
      blklen0, &
      neqv0, &
      eqvlen0, &
      workcoords0, &
      workadjmat0, &
      neqv1, &
      eqvlen1, &
      workcoords1, &
      workadjmat1, &
      workweights0, &
      maplist, &
      countlist, &
      nrec)

   if (react_flag) then

      offset(1) = 0
      do h = 1, nblk0 - 1
         offset(h+1) = offset(h) + blklen0(h)
      end do

      do h = 1, nblk0
         blkidx(offset(h)+1:offset(h)+blklen0(h)) = h
      end do

      mapping = maplist(:, 1)

      ! Align coordinates

      call rotate(natom1, workcoords1, leastrotquat(natom0, workweights0, workcoords0, workcoords1, mapping))

      ! Remove reactive bonds

      backadjmat0 = workadjmat0
      backadjmat1 = workadjmat1

      do i = 1, natom0
         do j = i + 1, natom0
            if (backadjmat0(i, j) .neqv. backadjmat1(mapping(i), mapping(j))) then
               call removereacbond(i, j, natom0, workznums0, workadjmat0, workadjmat1, mapping)
               h = blkidx(i)
               do k = offset(h) + 1, offset(h) + blklen0(h)
                  if (sum((workcoords0(:, i) - workcoords1(:, mapping(k)))**2) < 2.0 &
                     .or. sum((workcoords0(:, k) - workcoords1(:, mapping(i)))**2) < 2.0 &
                  ) then
!                     print *, '<', j, mapping(j), k, mapping(k)
                     call removereacbond(k, j, natom0, workznums0, workadjmat0, workadjmat1, mapping)
                  end if
               end do
               h = blkidx(j)
               do k = offset(h) + 1, offset(h) + blklen0(h)
                  if (sum((workcoords0(:, j) - workcoords1(:, mapping(k)))**2) < 2.0 &
                     .or. sum((workcoords0(:, k) - workcoords1(:, mapping(j)))**2) < 2.0 &
                  ) then
!                     print *, '>', i, mapping(i), k, mapping(k)
                     call removereacbond(i, k, natom0, workznums0, workadjmat0, workadjmat1, mapping)
                  end if
               end do
            end if
         end do
      end do

      ! Optimize the assignment to minimize AdjD/RMSD

      call optimize_assignment( &
         natom0, &
         nblk0, &
         blklen0, &
         neqv0, &
         eqvlen0, &
         workcoords0, &
         workadjmat0, &
         neqv1, &
         eqvlen1, &
         workcoords1, &
         workadjmat1, &
         workweights0, &
         maplist, &
         countlist, &
         nrec)

   end if

   ! Print coordinates with internal order

   mapping = maplist(:, 1)
   call rotate(natom1, workcoords1, leastrotquat(natom0, workweights0, workcoords0, workcoords1, mapping))
   open(unit=99, file='aligned_debug.mol2', action='write', status='replace')
   call adjmat2bonds(natom0, workadjmat0, nbond0, bonds0)
   call adjmat2bonds(natom1, workadjmat1, nbond1, bonds1)
!   call adjmat2bonds(natom1, workadjmat1(mapping, mapping), nbond1, bonds1)
   call writemol2(99, 'coords0', natom0, workznums0, workcoords0, nbond0, bonds0)
   call writemol2(99, 'coords1', natom1, workznums1, workcoords1, nbond1, bonds1)
!   call writemol2(99, 'coords1', natom1, workznums1(mapping), workcoords1(:, mapping), nbond1, bonds1)

   ! Reorder back to original atom ordering

   do i = 1, nrec
      maplist(:, i) = atomorder1(maplist(backorder0, i))
   end do

end subroutine

subroutine align_atoms( &
! Purpose: Align atoms0 and atoms1
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

subroutine removereacbond(i, j, natom, znums, adjmat0, adjmat1, mapping)
! Purpose: Remove reactive bonds
   integer, intent(in) :: i, j, natom
   integer, dimension(:), intent(in) :: znums
   logical, dimension(:, :), intent(inout) :: adjmat0, adjmat1
   integer, dimension(:), intent(in) :: mapping

   integer :: k
   integer, dimension(natom) :: nadj0, nadj1
   integer, dimension(maxcoord, natom) :: adjlist0, adjlist1

   ! Calculate adjacency lists

   call adjmat2list(natom, adjmat0, nadj0, adjlist0)
   call adjmat2list(natom, adjmat1, nadj1, adjlist1)

   adjmat0(i, j) = .false.
   adjmat0(j, i) = .false.
   adjmat1(mapping(i), mapping(j)) = .false.
   adjmat1(mapping(j), mapping(i)) = .false.
   if (znums(i) == 1) then
      do k = 1, nadj0(i)
         if (znums(adjlist0(k, i)) == 7 .or. znums(adjlist0(k, i)) == 8) then
            adjmat0(i, adjlist0(k, i)) = .false.
            adjmat0(adjlist0(k, i), i) = .false.
         end if
      end do
   end if
   if (znums(i) == 1) then
      do k = 1, nadj1(i)
         if (znums(adjlist1(k, i)) == 7 .or. znums(adjlist1(k, i)) == 8) then
            adjmat1(mapping(i), mapping(adjlist1(k, i))) = .false.
            adjmat1(mapping(adjlist1(k, i)), mapping(i)) = .false.
         end if
      end do
   end if
   if (znums(j) == 1) then
      do k = 1, nadj0(j)
         if (znums(adjlist0(k, j)) == 7 .or. znums(adjlist0(k, j)) == 8) then
            adjmat0(j, adjlist0(k, j)) = .false.
            adjmat0(adjlist0(k, j), j) = .false.
         end if
      end do
   end if
   if (znums(j) == 1) then
      do k = 1, nadj1(j)
         if (znums(adjlist1(k, j)) == 7 .or. znums(adjlist1(k, j)) == 8) then
            adjmat1(mapping(j), mapping(adjlist1(k, j))) = .false.
            adjmat1(mapping(adjlist1(k, j)), mapping(j)) = .false.
         end if
      end do
   end if

end subroutine

function get_rmsd(mol0, mol1) result(rmsd)
   class(Molecule), intent(in) :: mol0, mol1
   real(wp) :: rmsd

   rmsd = sqrt(squaredist(mol0%natom, mol0%get_weights(), mol0%get_coords(), mol1%get_coords(), identityperm(mol0%natom)) &
        / sum(mol0%get_weights()))

end function

function get_adjd(mol0, mol1) result(adjd)
   class(Molecule), intent(in) :: mol0, mol1
   integer :: adjd

   adjd = adjacencydiff(mol0%natom, mol0%adjmat, mol1%adjmat, identityperm(mol0%natom))

end function

end module
