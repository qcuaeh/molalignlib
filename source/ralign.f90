module ralign

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
!use utilities

implicit none

contains

logical function trial_stop_test(trials, matches) result(stop_test)
    integer, intent(in) :: trials, matches
    stop_test = trials < maxcount
end function

logical function match_stop_test(trials, matches) result(stop_test)
    integer, intent(in) :: trials, matches
    stop_test = matches < maxcount
end function

subroutine superpose(mol0, mol1, label0, label1, bonds0, bonds1, mapcount, atomaplist, rotmatlist)
! Purpose: Find best product to reactant map

    integer, intent(out) :: mapcount
    integer, dimension(:, :), intent(out) :: atomaplist
    real, dimension(:, :, :), intent(out) :: rotmatlist
    integer, dimension(:, :), intent(inout) :: bonds0, bonds1
    character(32), dimension(:), intent(inout) :: label0, label1
    real, dimension(:, :), intent(inout) :: mol0, mol1

    integer i
    integer natom0, natom1
    integer nbond0, nbond1
    integer nblock0, nblock1
    integer, dimension(:), allocatable :: seed
    integer, dimension(:), allocatable :: znum0, znum1
    integer, dimension(:), allocatable :: blocktype0, blocktype1
    integer, dimension(:), allocatable :: blocksize0, blocksize1
    integer, dimension(:), allocatable :: order0, reorder0, order1
    real, dimension(:), allocatable :: weights0
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

    ! Get number of atoms

    natom0 = size(mol0, dim=2)
    natom1 = size(mol1, dim=2)

    ! Get number of bonds

    nbond0 = size(bonds0, dim=2)
    nbond1 = size(bonds1, dim=2)

    ! Abort if number of atoms differ

    if (natom0 /= natom1) then
        call error('The molecules have different number of atoms!')
    end if

    ! Allocate arrays

    allocate(weights0(natom0))
    allocate(bias(natom0, natom1))
    allocate(order0(natom0), order1(natom1), reorder0(natom0))
    allocate(blocktype0(natom0), blocktype1(natom1))
    allocate(blocksize0(natom0), blocksize1(natom1))
    allocate(znum0(natom0), znum1(natom1))

    ! Get atomic numbers

    do i = 1, natom0
        label0(i) = upper(label0(i))
        znum0(i) = znum(label0(i))
    end do

    do i = 1, natom1
        label1(i) = upper(label1(i))
        znum1(i) = znum(label1(i))
    end do

    ! Calculate normalized weights

    weights0 = atomweight(znum0)/sum(atomweight(znum0))

    ! Group atoms by label

    call grouplabels(natom0, label0, nblock0, blocksize0, blocktype0)
    call grouplabels(natom1, label1, nblock1, blocksize1, blocktype1)

    ! Sort equivalent atoms in contiguous blocks

    order0 = sortorder(blocktype0, natom0)
    order1 = sortorder(blocktype1, natom1)
    reorder0 = inversemap(order0)

    label0 = label0(order0)
    label1 = label1(order1)

    mol0 = mol0(:, order0)
    mol1 = mol1(:, order1)

    znum0 = znum0(order0)
    znum1 = znum1(order1)

    weights0 = weights0(order0)

    ! Abort if there are incompatible atomic symbols

    if (any(label0 /= label1)) then
        call error('The molecules are not isomers!')
    end if

    ! Abort if there are incompatible atomic labels

    if (any(label0 /= label1)) then
        call error('There are incompatible atomic labels!')
    end if

    ! Translate molecules to their centroid

    call translate(natom0, centroid(natom0, weights0, mol0), mol0)
    call translate(natom1, centroid(natom1, weights0, mol1), mol1)

    ! Set bias for non equivalent atoms 

    call setadjbias(natom0, nblock0, blocksize0, mol0, mol1, bias)

    ! Initialize random number generator

    call init_random_seed(seed)
    call random_seed(put=seed)

    ! Remap atoms to minimize distance and difference

    if (maxcount >= 1) then
        call remapatoms(natom0, nblock0, blocksize0, mol0, mol1, weights0, mapcount, atomaplist, bias)
    else
        mapcount = 1
        atomaplist(:, 1) = [(i, i=1, natom0)]
        call print_header()
        call print_stats(0, 0, 0, 0., 0., 0., leastsquaredist(natom0, weights0, mol0, mol1, atomaplist(:, 1)))
        call print_footer(.false., .false., 0, 0)
    end if

    ! Calculate rotation quaternions

    do i = 1, mapcount
        rotmatlist(:, :, i) = quat2mat(leastrotquat(natom0, weights0, mol0, mol1, atomaplist(:, i)))
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

end module
