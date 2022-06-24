module superposition

use iso_fortran_env, only: output_unit
use iso_fortran_env, only: error_unit

use common
use options
use random
use sorting
use maputils
use rotation
use translation
use alignment
use remapping
use assortment
use chemistry
use messages

implicit none

contains


subroutine align(natom, atoms0, atoms1, labels0, labels1, maxrecord, nrecord, atomaplist, rotmatlist)
! Purpose: Superimpose coordinates of atom sets atoms0 and atoms1

integer, intent(in) :: natom, maxrecord
real(wp), dimension(:, :), intent(inout) :: atoms0, atoms1
character(32), dimension(:), intent(inout) :: labels0, labels1

integer, intent(out) :: nrecord
real(wp), dimension(:, :, :), intent(out) :: rotmatlist
integer, dimension(:, :), intent(out) :: atomaplist

integer i
integer, dimension(:), allocatable :: znum0, znum1
real(wp), dimension(:), allocatable :: weights
real(wp), dimension(nelem) :: atomweight

procedure (generic_test), pointer :: stop_test => null()

select case (weighting)
case ('none')
    atomweight = [(1., i=1, nelem)]
case ('mass')
    atomweight = stdmatom
case default
    call error('Invalid weighting option: '//trim(weighting))
end select

! Allocate arrays

allocate(weights(natom))
allocate(znum0(natom), znum1(natom))

! Get atomic numbers

do i = 1, natom
    labels0(i) = upper(labels0(i))
    labels1(i) = upper(labels1(i))
    znum0(i) = znum(labels0(i))
    znum1(i) = znum(labels1(i))
end do

! Calculate normalized weights

weights = atomweight(znum0)/sum(atomweight(znum0))

! Abort if there are incompatible atomic symbols

if (any(labels0 /= labels1)) then
    call error('The molecules are not isomers!')
end if

! Abort if there are incompatible atomic labels

if (any(labels0 /= labels1)) then
    call error('There are incompatible atomic labels!')
end if

! Translate molecules to their centroid

call translate(natom, centroid(natom, weights, atoms0), atoms0)
call translate(natom, centroid(natom, weights, atoms1), atoms1)

! Align atoms

nrecord = 1
atomaplist(:, 1) = [(i, i=1, natom)]
call print_header()
call print_stats(0, 0, 0, 0.0_wp, 0.0_wp, 0.0_wp, leastsquaredist(natom, weights, atoms0, atoms1, atomaplist(:, 1)))
call print_footer(.false., .false., 0, 0)

! Calculate rotation quaternions

do i = 1, nrecord
    rotmatlist(:, :, i) = quat2mat(leastrotquat(natom, weights, atoms0, atoms1, atomaplist(:, i)))
end do

end subroutine


subroutine realign(natom, atoms0, atoms1, labels0, labels1, maxrecord, nrecord, atomaplist, rotmatlist)
! Purpose: Superimpose coordinates of atom sets atoms0 and atoms1

integer, intent(in) :: natom, maxrecord
real(wp), dimension(:, :), intent(inout) :: atoms0, atoms1
character(32), dimension(:), intent(inout) :: labels0, labels1

integer, intent(out) :: nrecord
real(wp), dimension(:, :, :), intent(out) :: rotmatlist
integer, dimension(:, :), intent(out) :: atomaplist

integer i
integer nblock0, nblock1
integer, dimension(:), allocatable :: seed
integer, dimension(:), allocatable :: znum0, znum1
integer, dimension(:), allocatable :: blocktype0, blocktype1
integer, dimension(:), allocatable :: blocksize0, blocksize1
integer, dimension(:), allocatable :: order0, reorder0, order1
real(wp), dimension(:), allocatable :: weights
real(wp), dimension(:, :), allocatable :: bias
real(wp), dimension(nelem) :: atomweight

procedure (generic_test), pointer :: stop_test => null()

if (maxcount < 1) then
    call error('Max count must be greater than zero')
end if

if (converge) then
    stop_test => match_stop_test
else
    stop_test => trial_stop_test
end if

select case (weighting)
case ('none')
    atomweight = [(1., i=1, nelem)]
case ('mass')
    atomweight = stdmatom
case default
    call error('Invalid weighting option: '//trim(weighting))
end select

! Allocate arrays

allocate(weights(natom))
allocate(bias(natom, natom))
allocate(order0(natom), order1(natom), reorder0(natom))
allocate(blocktype0(natom), blocktype1(natom))
allocate(blocksize0(natom), blocksize1(natom))
allocate(znum0(natom), znum1(natom))

! Get atomic numbers

do i = 1, natom
    labels0(i) = upper(labels0(i))
    labels1(i) = upper(labels1(i))
    znum0(i) = znum(labels0(i))
    znum1(i) = znum(labels1(i))
end do

! Calculate normalized weights

weights = atomweight(znum0)/sum(atomweight(znum0))

! Group atoms by label

call grouplabels(natom, labels0, nblock0, blocksize0, blocktype0)
call grouplabels(natom, labels1, nblock1, blocksize1, blocktype1)

! Sort equivalent atoms in contiguous blocks

order0 = sortorder(blocktype0, natom)
order1 = sortorder(blocktype1, natom)
reorder0 = inversemap(order0)

labels0 = labels0(order0)
labels1 = labels1(order1)

atoms0 = atoms0(:, order0)
atoms1 = atoms1(:, order1)

znum0 = znum0(order0)
znum1 = znum1(order1)

weights = weights(order0)

! Abort if there are incompatible atomic symbols

if (any(labels0 /= labels1)) then
    call error('The molecules are not isomers!')
end if

! Abort if there are incompatible atomic labels

if (any(labels0 /= labels1)) then
    call error('There are incompatible atomic labels!')
end if

! Translate molecules to their centroid

call translate(natom, centroid(natom, weights, atoms0), atoms0)
call translate(natom, centroid(natom, weights, atoms1), atoms1)

! Set bias for non equivalent atoms 

call setadjbias(natom, nblock0, blocksize0, atoms0, atoms1, bias)

! Initialize random number generator

call init_random_seed(seed)
call random_seed(put=seed)

! Remap atoms to minimize distance and difference

call remapatoms(natom, nblock0, blocksize0, atoms0, atoms1, weights, bias, maxrecord, nrecord, atomaplist, &
    stop_test)

! Calculate rotation quaternions

do i = 1, nrecord
    rotmatlist(:, :, i) = quat2mat(leastrotquat(natom, weights, atoms0, atoms1, atomaplist(:, i)))
end do

! Restore original atom order

labels0 = labels0(reorder0)
labels1 = labels1(reorder0)

atoms0 = atoms0(:, reorder0)
atoms1 = atoms1(:, reorder0)

do i = 1, nrecord
    atomaplist(:, i) = atomaplist(reorder0, i)
end do

end subroutine

end module
