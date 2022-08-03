!module library
!contains

subroutine remap( &
! Purpose: Check and optimize mapping
    natom0, natom1, znums0, znums1, types0, types1, weights0, weights1, &
    coords0, coords1, maxrecords, nrecord, atomaplist, countlist &
)

use iso_fortran_env, only: output_unit
use iso_fortran_env, only: error_unit

use options
use random
use sorting
use chemdata
use maputils
use rotation
use translation
use remapping
use assortment
use alignment

implicit none

integer, intent(in) :: natom0, natom1, maxrecords
integer, dimension(natom0), intent(in) :: znums0, types0
integer, dimension(natom1), intent(in) :: znums1, types1
real, intent(in) :: coords0(3, natom0)
real, intent(in) :: coords1(3, natom1)
real, intent(in) :: weights0(natom0)
real, intent(in) :: weights1(natom1)
integer, intent(out) :: nrecord
integer, intent(out) :: atomaplist(natom0, maxrecords)
integer, intent(out) :: countlist(maxrecords)

integer i
integer nblock0, nblock1
integer, dimension(:), allocatable :: seed
integer, dimension(:), allocatable :: order0, order1
integer, dimension(:), allocatable :: blockidx0, blockidx1
integer, dimension(:), allocatable :: blocksize0, blocksize1
real, dimension(3) :: center0, center1

procedure(test), pointer :: trial_test => null()
procedure(test), pointer :: match_test => null()

! Check number of atoms

if (natom0 /= natom1) then
    write (error_unit, '(a)') 'Error: The number of atoms does not match'
    stop
end if

! Allocate arrays

allocate(order0(natom0), order1(natom1))
allocate(blockidx0(natom0), blockidx1(natom1))
allocate(blocksize0(natom0), blocksize1(natom1))

! Select trial exit test

if (bounded) then
    trial_test => lower_than
else
    trial_test => dummy_test
end if

! Select match exit test

if (counting) then
    match_test => lower_than
else
    match_test => dummy_test
end if

! Group atoms by label

call getblocks(natom0, znums0, types0, weights0, nblock0, blocksize0, blockidx0, order0)
call getblocks(natom1, znums1, types1, weights1, nblock1, blocksize1, blockidx1, order1)

! Abort if there are incompatible atomic symbols

if (any(znums0(order0) /= znums1(order1))) then
    write (error_unit, '(a)') 'Error: Clusters are not isomers'
    stop
end if

! Abort if there are incompatible atomic types

if (any(types0(order0) /= types1(order1))) then
    write (error_unit, '(a)') 'Error: There are conflicting atomic types'
    stop
end if

! Abort if there are incompatible atomic weights

if (any(weights0(order0) /= weights1(order1))) then
    write (error_unit, '(a)') 'Error: There are conflicting atomic weights'
    stop
end if

! Calculate centroids

center0 = centroid(natom0, weights0, coords0)
center1 = centroid(natom1, weights1, coords1)

! Initialize random number generator

call init_random_seed(seed)
call random_seed(put=seed)

! Remap atoms to minimize distance and difference

call optimize_mapping( &
    natom0, nblock0, blocksize0, weights0(order0), &
    centered(natom0, coords0(:, order0), center0), &
    centered(natom1, coords1(:, order1), center1), &
    maxrecords, nrecord, atomaplist, countlist, &
    trial_test, match_test &
)

! Reorder back to original atom ordering

do i = 1, nrecord
    atomaplist(:, i) = order1(atomaplist(inversemap(order0), i))
end do

end subroutine

! Purpose: Superimpose coordinates of atom sets coords0 and coords1
subroutine align( &
    natom0, natom1, znums0, znums1, types0, types1, weights0, weights1, &
    coords0, coords1, travec, rotmat)

use iso_fortran_env, only: output_unit
use iso_fortran_env, only: error_unit

use options
use chemdata
use maputils
use rotation
use translation
use alignment

implicit none

integer, intent(in) :: natom0, natom1
integer, dimension(natom0), intent(in) :: znums0, types0
integer, dimension(natom1), intent(in) :: znums1, types1
real, intent(in) :: coords0(3, natom0)
real, intent(in) :: coords1(3, natom1)
real, intent(in) :: weights0(natom0)
real, intent(in) :: weights1(natom1)
real, intent(out) :: travec(3)
real, intent(out) :: rotmat(3, 3)

real, dimension(3) :: center0, center1
integer i
! Check number of atoms

if (natom0 /= natom1) then
    write (error_unit, '(a)') 'Error: The number of atoms does not match'
    stop
end if

! Abort if there are incompatible atomic symbols

if (any(znums0 /= znums1)) then
    write (error_unit, '(a)') 'Error: Clusters are not isomers or atoms are not properly ordered'
    stop
end if

! Abort if there are incompatible atomic types

if (any(types0 /= types1)) then
    write (error_unit, '(a)') 'Error: There are conflicting atomic types or atoms are not properly ordered'
    stop
end if

! Abort if there are incompatible atomic weights

if (any(weights0 /= weights1)) then
    write (error_unit, '(a)') 'Error: There are conflicting atomic weights'
    stop
end if

! Calculate centroids

center0 = centroid(natom0, weights0, coords0)
center1 = centroid(natom1, weights1, coords1)

! Calculate optimal rotation matrix

rotmat = rotquat2rotmat(leastrotquat( &
    natom0, weights0, &
    centered(natom0, coords0, center0), &
    centered(natom1, coords1, center1), &
    identitymap(natom0) &
))

! Calculate optimal translation vector

travec = center0 - matmul(rotmat, center1)

end subroutine

!end module
