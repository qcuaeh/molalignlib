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
use types
use bounds
use random
use sorting
use discrete
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
subroutine remap_atoms( &
   mol0, &
   mol1, &
   maplist, &
   countlist, &
   nrec, &
   error)

   type(Molecule), intent(inout) :: mol0, mol1
   integer, dimension(:, :), intent(inout) :: maplist
   integer, dimension(:), intent(inout) :: countlist
   integer, intent(out) :: nrec, error

   integer :: i
   integer, dimension(:), allocatable :: atomorder0, atomorder1
   integer, dimension(:), allocatable :: backorder0, backorder1
   real(wp) :: travec0(3), travec1(3)

   integer, dimension(mol0%natom) :: mapping
   integer :: nbond0, bonds0(2, maxcoord*mol0%natom)
   integer :: nbond1, bonds1(2, maxcoord*mol1%natom)


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

   if (mol0%natom /= mol1%natom) then
      write (stderr, '(a)') 'Error: The molecules are not isomers'
      error = 1
      return
   end if

   ! Allocate arrays

   allocate(atomorder0(mol0%natom), atomorder1(mol1%natom))
   allocate(backorder0(mol0%natom), backorder1(mol1%natom))

   ! Get atom order

   atomorder0 = sorted_order(mol0%get_eqvids(), mol0%natom)
   atomorder1 = sorted_order(mol1%get_eqvids(), mol1%natom)

   ! Get inverse atom order

   backorder0 = inverse_permut(atomorder0)
   backorder1 = inverse_permut(atomorder1)

   ! Reorder data arrays

   call mol0%permutate_atoms(atomorder0)
   call mol1%permutate_atoms(atomorder1)

   ! Abort if molecules are not isomers

   if (any(mol0%get_znums() /= mol1%get_znums())) then
      write (stderr, '(a)') 'Error: The molecules are not isomers'
      error = 1
      return
   end if

   ! Abort if there are conflicting atomic types

   if (any(mol0%get_blkids() /= mol1%get_blkids())) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      error = 1
      return
   end if

   ! Abort if there are conflicting weights

   if (any(abs(mol0%get_weights() - mol1%get_weights()) > 1.E-6)) then
      write (stderr, '(a)') 'Error: There are conflicting weights'
      error = 1
      return
   end if

   ! Mirror coordinates

   if (mirror_flag) then
      call mol1%mirror_coords()
   end if

   ! Calculate centroids

   travec0 = -mol0%get_center()
   travec1 = -mol1%get_center()

   ! Center coordinates at the centroids

   call mol0%translate_coords(travec0)
   call mol1%translate_coords(travec1)

   ! Initialize random number generator

   call initialize_random()

   ! Optimize assignment to minimize the AdjD and RMSD

    call optimize_mapping(mol0, mol1, maplist, countlist, nrec)

   ! Debond reactive sites and reoptimize assignment

   if (reac_flag) then
      call find_reactive_sites(mol0, mol1, maplist(:, 1))
      call assort_neighbors(mol0)
      call assort_neighbors(mol1)
      call optimize_mapping(mol0, mol1, maplist, countlist, nrec)
   end if

!   ! Print coordinates with internal order

!   mapping = maplist(:, 1)
!   call rotate(natom1, workcoords1, leastrotquat(natom0, workweights0, workcoords0, workcoords1, mapping))
!   open(unit=99, file='aligned_debug.mol2', action='write', status='replace')
!   call adjmat2bonds(natom0, workadjmat0, nbond0, bonds0)
!   call adjmat2bonds(natom1, workadjmat1, nbond1, bonds1)
!!   call adjmat2bonds(natom1, workadjmat1(mapping, mapping), nbond1, bonds1)
!   call writemol2(99, 'coords0', natom0, workznums0, workcoords0, nbond0, bonds0)
!   call writemol2(99, 'coords1', natom1, workznums1, workcoords1, nbond1, bonds1)
!!   call writemol2(99, 'coords1', natom1, workznums1(mapping), workcoords1(:, mapping), nbond1, bonds1)

   ! Reorder back to original atom ordering

   do i = 1, nrec
      maplist(:, i) = atomorder1(maplist(backorder0, i))
   end do

end subroutine

subroutine align_atoms( &
! Purpose: Align atoms0 and atoms1
   mol0, &
   mol1, &
   travec, &
   rotmat, &
   error)

   type(Molecule), intent(inout) :: mol0, mol1
   real(wp), intent(out) :: travec(3), rotmat(3, 3)
   real(wp) :: travec0(3), travec1(3), rotquat(4)
   integer, intent(out) :: error

   ! Set error code to 0 by default

   error = 0

   ! Abort if molecules have different number of atoms

   if (mol0%natom /= mol1%natom) then
      write (stderr, '(a)') 'Error: The molecules are not isomers'
      error = 1
      return
   end if

   ! Abort if molecules are not isomers

   if (any(sorted(mol0%get_znums(), mol0%natom) /= sorted(mol1%get_znums(), mol1%natom))) then
      write (stderr, '(a)') 'Error: The molecules are not isomers'
      error = 1
      return
   end if

   ! Abort if atoms are not ordered

   if (any(mol0%get_znums() /= mol1%get_znums())) then
      write (stderr, '(a)') 'Error: The atoms are not in the same order'
      error = 1
      return
   end if

   ! Abort if there are conflicting atomic types

   if (any(sorted(mol0%get_blkids(), mol0%natom) /= sorted(mol1%get_blkids(), mol1%natom))) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      error = 1
      return
   end if

   ! Abort if atomic types are not ordered

   if (any(mol0%get_blkids() /= mol1%get_blkids())) then
      write (stderr, '(a)') 'Error: Atomic types are not in the same order'
      error = 1
      return
   end if

   ! Calculate centroids

!CZGC: nuevo llamado (por actualizar en molalignlib/types.f90, LAZH):
    travec0 = -mol0%get_center()
    travec1 = -mol1%get_center()
!CZGC: llamado anterior
!   travec0 = -center_coords(mol0%natom, mol0%get_weights(), mol0%get_coords())
!   travec1 = -center_coords(mol1%natom, mol1%get_weights(), mol1%get_coords())

   ! Calculate optimal rotation matrix

   rotquat = leastrotquat( &
      mol0%natom, &
      mol0%get_weights(), &
      translated(mol0%natom, mol0%get_coords(), travec0), &
      translated(mol1%natom, mol1%get_coords(), travec1), &
      identity_perm(mol0%natom) &
   )

   rotmat = quat2rotmat(rotquat)

   ! Calculate optimal translation vector

   travec = matmul(rotmat, travec1) - travec0

end subroutine

function get_rmsd(mol0, mol1) result(rmsd)
   type(Molecule), intent(in) :: mol0, mol1
   real(wp) :: rmsd

   rmsd = sqrt(squaredist(mol0%natom, mol0%get_weights(), mol0%get_coords(), mol1%get_coords(), identity_perm(mol0%natom)) &
        / sum(mol0%get_weights()))

end function

function get_adjd(mol0, mol1) result(adjd)
   type(Molecule), intent(in) :: mol0, mol1
   integer :: adjd

   adjd = adjacencydiff(mol0%natom, mol0%get_adjmat(), mol1%get_adjmat(), identity_perm(mol0%natom))

end function

end module
