subroutine superpose(natom, mol0, mol1, label0, label1, mapcount, atomaplist, rotmatlist)
! Purpose: Superimpose coordinates of atom sets mol0 and mol1

use iso_fortran_env, only: output_unit
use iso_fortran_env, only: error_unit

use options
use random
use sorting
use rotation
use translation
use alignment
use remapping
use assortment
use chemistry
use messages

implicit none

integer, intent(in) :: natom
integer, intent(out) :: mapcount
real, dimension(3, 3, maxcount), intent(out) :: rotmatlist
integer, dimension(natom, maxcount), intent(out) :: atomaplist
character(32), dimension(natom), intent(inout) :: label0, label1
real, dimension(3, natom), intent(inout) :: mol0, mol1

integer i
integer nblock0, nblock1
integer, dimension(:), allocatable :: seed
integer, dimension(:), allocatable :: znum0, znum1
integer, dimension(:), allocatable :: blocktype0, blocktype1
integer, dimension(:), allocatable :: blocksize0, blocksize1
integer, dimension(:), allocatable :: order0, reorder0, order1
real, dimension(:), allocatable :: weights
real, dimension(:, :), allocatable :: bias
real, dimension(nelem) :: atomweight

if (maxcount < 0) then
    call error('Argument is missing!')
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
    label0(i) = upper(label0(i))
    label1(i) = upper(label1(i))
    znum0(i) = znum(label0(i))
    znum1(i) = znum(label1(i))
end do

! Calculate normalized weights

weights = atomweight(znum0)/sum(atomweight(znum0))

! Group atoms by label

call grouplabels(natom, label0, nblock0, blocksize0, blocktype0)
call grouplabels(natom, label1, nblock1, blocksize1, blocktype1)

! Sort equivalent atoms in contiguous blocks

order0 = sortorder(blocktype0, natom)
order1 = sortorder(blocktype1, natom)
reorder0 = inversemap(order0)

label0 = label0(order0)
label1 = label1(order1)

mol0 = mol0(:, order0)
mol1 = mol1(:, order1)

znum0 = znum0(order0)
znum1 = znum1(order1)

weights = weights(order0)

! Abort if there are incompatible atomic symbols

if (any(label0 /= label1)) then
    call error('The molecules are not isomers!')
end if

! Abort if there are incompatible atomic labels

if (any(label0 /= label1)) then
    call error('There are incompatible atomic labels!')
end if

! Translate molecules to their centroid

call translate(natom, centroid(natom, weights, mol0), mol0)
call translate(natom, centroid(natom, weights, mol1), mol1)

! Set bias for non equivalent atoms 

call setadjbias(natom, nblock0, blocksize0, mol0, mol1, bias)

! Initialize random number generator

call init_random_seed(seed)
call random_seed(put=seed)

! Remap atoms to minimize distance and difference

if (maxcount >= 1) then
    call remapatoms(natom, nblock0, blocksize0, mol0, mol1, weights, mapcount, atomaplist, bias)
else
    mapcount = 1
    atomaplist(:, 1) = [(i, i=1, natom)]
    call print_header()
    call print_stats(0, 0, 0, 0., 0., 0., leastsquaredist(natom, weights, mol0, mol1, atomaplist(:, 1)))
    call print_footer(.false., .false., 0, 0)
end if

! Calculate rotation quaternions

do i = 1, mapcount
    rotmatlist(:, :, i) = quat2mat(leastrotquat(natom, weights, mol0, mol1, atomaplist(:, i)))
end do

! Restore original atom order

label0 = label0(reorder0)
label1 = label1(reorder0)

mol0 = mol0(:, reorder0)
mol1 = mol1(:, reorder0)

do i = 1, mapcount
    atomaplist(:, i) = atomaplist(reorder0, i)
end do

end subroutine
