!module superposition
!contains

! Purpose: Superimpose coordinates of atom sets coords0 and coords1
subroutine superpose( &
    natom0, natom1, znums0, znums1, types0, types1, weights0, weights1, &
    coords0, coords1, maxrecord, nrecord, atomaplist, countlist &
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
use alignment
use remapping
use assortment

implicit none

integer, intent(in) :: natom0, natom1, maxrecord
integer, dimension(natom0), intent(in) :: znums0, types0
integer, dimension(natom1), intent(in) :: znums1, types1
real(wp), intent(in) :: weights0(natom0), weights1(natom1)
real(wp), intent(in) :: coords0(3, natom0), coords1(3, natom1)
integer, intent(out) :: nrecord
integer, intent(out) :: atomaplist(natom0, maxrecord)
integer, intent(out) :: countlist(maxrecord)

integer i
integer nblock0, nblock1
integer, dimension(:), allocatable :: seed
integer, dimension(:), allocatable :: order0, order1
integer, dimension(:), allocatable :: reorder0, reorder1
integer, dimension(:), allocatable :: blockidx0, blockidx1
integer, dimension(:), allocatable :: blocksize0, blocksize1
real(wp), dimension(3) :: center0, center1

procedure (generic_test), pointer :: trial_test => null()
procedure (generic_test), pointer :: match_test => null()

! Check number of atoms

if (natom0 /= natom1) then
    write (error_unit, '(a)') 'Error: The molecules have different number of atoms'
    stop
end if

! Allocate arrays

allocate(order0(natom0), order1(natom1))
allocate(reorder0(natom0), reorder1(natom1))
allocate(blockidx0(natom0), blockidx1(natom1))
allocate(blocksize0(natom0), blocksize1(natom1))

! Select trial exit test

if (trialing) then
    trial_test => lowerthan
else
    trial_test => trueconst
end if

! Select match exit test

if (matching) then
    match_test => lowerthan
else
    match_test => trueconst
end if

if (remap) then

    ! Group atoms by label

    call getblocks(natom0, znums0, types0, nblock0, blocksize0, blockidx0)
    call getblocks(natom1, znums1, types1, nblock1, blocksize1, blockidx1)

    ! Get contiguous label order

    order0 = sortorder(blockidx0, natom0)
    order1 = sortorder(blockidx1, natom1)

    ! Get inverse order

    reorder0 = inversemap(order0)
    reorder1 = inversemap(order1)

    ! Abort if there are incompatible atomic symbols

    if (any(znums0(order0) /= znums1(order1))) then
        write (error_unit, '(a)') 'Error: The molecules are not isomers'
        stop
    end if

    ! Abort if there are incompatible atomic types

    if (any(types0(order0) /= types1(order1))) then
        write (error_unit, '(a)') 'Error: There are incompatible atomic types'
        stop
    end if

    ! Abort if there are incompatible atomic weights

    if (any(weights0(order0) /= weights1(order1))) then
        write (error_unit, '(a)') 'Error: There are incompatible atomic weights'
        stop
    end if

else

    ! Abort if there are incompatible atomic symbols

    if (any(znums0 /= znums1)) then
        write (error_unit, '(a)') 'Error: The molecules are not isomers'
        stop
    end if

    ! Abort if there are incompatible atomic types

    if (any(types0 /= types1)) then
        write (error_unit, '(a)') 'Error: There are incompatible atomic types'
        stop
    end if

    ! Abort if there are incompatible atomic weights

    if (any(weights0 /= weights1)) then
        write (error_unit, '(a)') 'Error: There are incompatible atomic weights'
        stop
    end if

end if

! Calculate centroids

center0 = centroid(natom0, weights0, coords0)
center1 = centroid(natom1, weights1, coords1)

if (remap) then

    ! Initialize random number generator

    call init_random_seed(seed)
    call random_seed(put=seed)

    ! Remap atoms to minimize distance and difference

    call remapatoms( &
        natom0, nblock0, blocksize0, weights0(order0), &
        translated(natom0, center0, coords0(:, order0)), &
        translated(natom1, center1, coords1(:, order1)), &
        maxrecord, nrecord, atomaplist, countlist, &
        trial_test, match_test &
    )

    ! Restore original atom ordering

    do i = 1, nrecord
        atomaplist(:, i) = atomaplist(reorder0, i)
    end do

else

    nrecord = 1
    atomaplist(:, 1) = [(i, i=1, natom0)]

    call print_header()
    call print_stats( &
        0, 0, 0, 0.0_wp, 0.0_wp, 0.0_wp, &
        leastsquaredist( &
            natom0, weights0, &
            translated(natom0, center0, coords0), &
            translated(natom1, center1, coords1), &
            atomaplist(:, 1) &
        ) &
    )
    call print_footer(.false., .false., 0, 0)

end if

end subroutine

!end module
