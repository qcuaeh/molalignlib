module library

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

contains

subroutine remap(natom0, natom1, znums0, znums1, types0, types1, &
    coords0, coords1, weights0, records, nrec, maplist, mapcount, mindist)
! Purpose: Check and optimize mapping

    integer, intent(in) :: natom0, natom1, records
    integer, dimension(natom0), intent(in) :: znums0, types0
    integer, dimension(natom1), intent(in) :: znums1, types1
    real, intent(in) :: coords0(3, natom0)
    real, intent(in) :: coords1(3, natom1)
    real, intent(in) :: weights0(natom0)
    integer, intent(out) :: nrec
    integer, intent(out) :: maplist(natom0, records)
    integer, intent(out) :: mapcount(records)
    real, intent(out) :: mindist(records)

    integer i, h, offset
    integer nblock0, nblock1
    integer, dimension(:), allocatable :: seed
    integer, dimension(:), allocatable :: order0, order1
    integer, dimension(:), allocatable :: blockidx0, blockidx1
    integer, dimension(:), allocatable :: blocksize0, blocksize1
    real, dimension(3) :: center0, center1

    procedure(test), pointer :: abort_test => null()
    procedure(test), pointer :: converge_test => null()

    ! Check number of atoms

    if (natom0 /= natom1) then
        write (error_unit, '(a)') 'Error: The number of atoms does not match'
        stop
    end if

    ! Allocate arrays

    allocate(order0(natom0), order1(natom1))
    allocate(blockidx0(natom0), blockidx1(natom1))
    allocate(blocksize0(natom0), blocksize1(natom1))

    ! Associate abortion test

    if (abort) then
        abort_test => lower_than
    else
        abort_test => dummy_test
    end if

    ! Associate convergence test

    if (converge) then
        converge_test => lower_than
    else
        converge_test => dummy_test
    end if

    if (associated(converge_test, dummy_test) .and. &
        associated(abort_test, dummy_test)) then
        write (error_unit, '(a)') 'Error: There is no stopping test associated'
        stop
    end if

    ! Group atoms by label

    call getblocks(natom0, znums0, types0, nblock0, blocksize0, blockidx0, order0)
    call getblocks(natom1, znums1, types1, nblock1, blocksize1, blockidx1, order1)

    ! Abort if there are incompatible atoms

    if (any(znums0(order0) /= znums1(order1))) then
        write (error_unit, '(a)') 'Error: Clusters are not isomers'
        stop
    end if

    ! Abort if there are incompatible types

    if (any(types0(order0) /= types1(order1))) then
        write (error_unit, '(a)') 'Error: There are conflicting atomic types'
        stop
    end if

    ! Abort if weights differ within a block

    offset = 0
    do h = 1, nblock0
        do i = 2, blocksize0(h)
            if (weights0(order0(offset+i)) /= weights0(order0(offset+1))) then
                write (error_unit, '(a)') 'Error: All atoms within a block must have the same weight'
                stop
            end if
        end do
        offset = offset + blocksize0(h)
    end do

    ! Calculate centroids

    center0 = centroid(natom0, weights0, coords0)
    center1 = centroid(natom1, weights0, coords1)

    ! Initialize random number generator

    call init_random_seed(seed)

    ! Remap atoms to minimize distance and difference

    call optimize_mapping(natom0, nblock0, blocksize0, weights0(order0), &
        centered(natom0, coords0(:, order0), center0), &
        centered(natom1, coords1(:, order1), center1), &
        records, nrec, maplist, mapcount, mindist, abort_test, converge_test)

    ! Reorder back to original atom ordering

    do i = 1, nrec
        maplist(:, i) = order1(maplist(inversemap(order0), i))
    end do

end subroutine

subroutine align(natom0, natom1, znums0, znums1, types0, types1, &
    coords0, coords1, weights0, travec, rotmat)
! Purpose: Superimpose coordinates of atom sets coords0 and coords1

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

    ! Calculate centroids

    center0 = centroid(natom0, weights0, coords0)
    center1 = centroid(natom1, weights0, coords1)

    ! Calculate optimal rotation matrix

    rotmat = rotquat2rotmat(leastrotquat(natom0, weights0, &
        centered(natom0, coords0, center0), &
        centered(natom1, coords1, center1), &
        identitymap(natom0)))

    ! Calculate optimal translation vector

    travec = center0 - matmul(rotmat, center1)

end subroutine

end module
