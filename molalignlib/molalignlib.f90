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

module molalignlib
use parameters
use settings

implicit none

contains

subroutine set_bias_flag(boolval)
   logical, intent(in) :: boolval
   bias_flag = boolval
end subroutine

subroutine set_conv_flag(boolval)
   logical, intent(in) :: boolval
   iter_flag = boolval
end subroutine

subroutine set_test_flag(boolval)
   logical, intent(in) :: boolval
   test_flag = boolval
end subroutine

subroutine set_live_flag(boolval)
   logical, intent(in) :: boolval
   live_flag = boolval
end subroutine

subroutine set_free_flag(boolval)
   logical, intent(in) :: boolval
   free_flag = boolval
end subroutine

subroutine set_debug_flag(boolval)
   logical, intent(in) :: boolval
   debug_flag = boolval
end subroutine

subroutine set_max_count(intval)
   integer, intent(in) :: intval
   max_count = intval
end subroutine

subroutine set_max_trials(intval)
   integer, intent(in) :: intval
   max_trials = intval
end subroutine

subroutine set_bias_tol(realval)
   real(wp), intent(in) :: realval
   bias_tol = realval
end subroutine

subroutine set_bias_scale(realval)
   real(wp), intent(in) :: realval
   bias_scale = realval
end subroutine

! Purpose: Check and optimize mapping
subroutine assign_atoms( &
   natom0, &
   natom1, &
   znums0,  &
   znums1, &
   types0, &
   types1, &
   coords0, &
   coords1, &
   weights0, &
   nrec, &
   nmap, &
   maplist, &
   countlist, &
   rmsdlist &
)

   use random
   use sorting
   use translation
   use assortment
   use assignment

   integer, intent(in) :: natom0, natom1, nrec
   integer, dimension(natom0), intent(in) :: znums0, types0
   integer, dimension(natom1), intent(in) :: znums1, types1
   real(wp), intent(in) :: coords0(3, natom0)
   real(wp), intent(in) :: coords1(3, natom1)
   real(wp), intent(in) :: weights0(natom0)
   integer, intent(out) :: nmap
   integer, intent(out) :: maplist(natom0, nrec)
   integer, intent(out) :: countlist(nrec)
   real(wp), intent(out) :: rmsdlist(nrec)

   integer :: i, h, offset
   integer :: nblock0, nblock1
   integer, dimension(:), allocatable :: atomorder0, atomorder1
   integer, dimension(:), allocatable :: atomunorder0, atomunorder1
   integer, dimension(:), allocatable :: blockidx0, blockidx1
   integer, dimension(:), allocatable :: blocksize0, blocksize1
   real(wp), dimension(3) :: center0, center1
   real(wp) :: weights1(natom1)

   ! Abort if clusters have different number of atoms

   if (natom0 /= natom1) then
      write (error_unit, '(a)') 'Error: The clusters have different number of atoms'
      stop
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

   atomunorder0 = inversemap(atomorder0)
   atomunorder1 = inversemap(atomorder1)

   ! Abort if clusters are not isomers

   if (any(znums0(atomorder0) /= znums1(atomorder1))) then
      write (error_unit, '(a)') 'Error: The clusters are not isomers'
      stop
   end if

   ! Abort if there are conflicting types

   if (any(types0(atomorder0) /= types1(atomorder1))) then
      write (error_unit, '(a)') 'Error: There are conflicting atom types'
      stop
   end if

   ! Assign weights for atoms1

   weights1 = weights0(atomunorder1(atomorder0))

   ! Abort if there are conflicting weights

   if (any(weights0(atomorder0) /= weights1(atomorder1))) then
      write (error_unit, '(a)') 'Error: There are conflicting weights'
      stop
   end if

   offset = 0
   do h = 1, nblock0
      do i = 2, blocksize0(h)
         if (weights0(atomorder0(offset+i)) /= weights0(atomorder0(offset+1))) then
            write (error_unit, '(a)') 'Error: Atoms of the same type must weight the same'
            stop
         end if
      end do
      offset = offset + blocksize0(h)
   end do

   ! Calculate centroids

   center0 = centroid(natom0, weights0, coords0)
   center1 = centroid(natom1, weights1, coords1)

   ! Initialize random number generator

   call initialize_random()

   ! Remap atoms to minimize distance and difference

   call optimize_assignment( &
      natom0, nblock0, blocksize0, &
      weights0(atomorder0), &
      centered(natom0, coords0(:, atomorder0), center0), &
      centered(natom1, coords1(:, atomorder1), center1), &
      nrec, nmap, maplist, countlist, rmsdlist &
   )

   ! Reorder back to original atom ordering

   do i = 1, nmap
      maplist(:, i) = atomorder1(maplist(atomunorder0, i))
   end do

end subroutine

! Purpose: Superimpose coordinates of atom sets coords0 and coords1
subroutine align_atoms( &
   natom0, &
   natom1, &
   znums0, &
   znums1, &
   types0, &
   types1, &
   coords0, &
   coords1, &
   weights0, &
   rmsd, &
   aligned1 &
)

   use sorting
   use rotation
   use translation
   use alignment

   implicit none

   integer, intent(in) :: natom0, natom1
   integer, dimension(natom0), intent(in) :: znums0, types0
   integer, dimension(natom1), intent(in) :: znums1, types1
   real(wp), intent(in) :: weights0(natom0)
   real(wp), intent(in) :: coords0(3, natom0)
   real(wp), intent(in) :: coords1(3, natom1)
   real(wp), intent(out) :: aligned1(3, natom1)
   real(wp), intent(out) :: rmsd

   real(wp) :: weights1(natom1)
   real(wp) :: center0(3), center1(3)
   real(wp) :: travec(3), rotmat(3, 3)

   ! Abort if clusters have different number of atoms

   if (natom0 /= natom1) then
      write (error_unit, '(a)') 'Error: The clusters have different number of atoms'
      stop
   end if

   ! Abort if clusters are not isomers

   if (any(sorted(znums0, natom0) /= sorted(znums1, natom1))) then
      write (error_unit, '(a)') 'Error: The clusters are not isomers'
      stop
   end if

   ! Abort if atoms are not ordered

   if (any(znums0 /= znums1)) then
      write (error_unit, '(a)') 'Error: The atoms are not in the same order'
      stop
   end if

   ! Abort if there are conflicting types

   if (any(sorted(types0, natom0) /= sorted(types1, natom1))) then
      write (error_unit, '(a)') 'Error: There are conflicting atom types'
      stop
   end if

   ! Abort if types are not ordered

   if (any(types0 /= types1)) then
      write (error_unit, '(a)') 'Error: The atom types are not in the same order'
      stop
   end if

   ! Assign weights for atoms1

   weights1 = weights0

   ! Calculate centroids

   center0 = centroid(natom0, weights0, coords0)
   center1 = centroid(natom1, weights1, coords1)

   ! Calculate optimal rotation matrix

   rotmat = rotquat2rotmat(leastrotquat(natom0, weights0, &
      centered(natom0, coords0, center0), &
      centered(natom1, coords1, center1), &
      identitymap(natom0)))

   ! Calculate optimal translation vector

   travec = center0 - matmul(rotmat, center1)

   aligned1 = translated(natom1, rotated(natom1, coords1, rotmat), travec)
   rmsd = sqrt(squaredist(natom0, weights0, coords0, aligned1, identitymap(natom0)))

end subroutine

end module
