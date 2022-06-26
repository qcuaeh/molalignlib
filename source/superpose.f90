subroutine superpose(natom, maxrecord, labels0, labels1, atoms0, atoms1, &
    center0, center1, nrecord, atomaplist, rotmatlist)
! Purpose: Superimpose coordinates of atom sets atoms0 and atoms1

use iso_fortran_env, only: output_unit
use iso_fortran_env, only: error_unit

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

integer, intent(in) :: natom, maxrecord
real(wp), dimension(3, natom), intent(in) :: atoms0, atoms1
character*32, dimension(natom), intent(in) :: labels0, labels1
integer, intent(out) :: nrecord
real(wp), dimension(3), intent(out) :: center0, center1
real(wp), dimension(3, 3, maxrecord), intent(out) :: rotmatlist
integer, dimension(natom, maxrecord), intent(out) :: atomaplist

integer i
integer nblock0, nblock1
integer, dimension(:), allocatable :: seed
integer, dimension(:), allocatable :: znum0, znum1
integer, dimension(:), allocatable :: order0, order1
integer, dimension(:), allocatable :: reorder0, reorder1
integer, dimension(:), allocatable :: blocktype0, blocktype1
integer, dimension(:), allocatable :: blocksize0, blocksize1
real(wp), dimension(:), allocatable :: weights
real(wp), dimension(:, :), allocatable :: bias
real(wp), dimension(nelem) :: atomweight

procedure (generic_test), pointer :: trial_test => null()
procedure (generic_test), pointer :: match_test => null()

if (trialing) then
    if (maxtrial < 1) then
        call error('maxtrial must be greater than zero')
    end if
    trial_test => lowerthan
else
    trial_test => trueconst
end if

if (matching) then
    if (maxmatch < 1) then
        call error('maxmatch must be greater than zero')
    end if
    match_test => lowerthan
else
    match_test => trueconst
end if

! Allocate arrays

allocate(weights(natom))
allocate(bias(natom, natom))
allocate(order0(natom), order1(natom))
allocate(reorder0(natom), reorder1(natom))
allocate(blocktype0(natom), blocktype1(natom))
allocate(blocksize0(natom), blocksize1(natom))
allocate(znum0(natom), znum1(natom))

! Get atomic numbers

do i = 1, natom
    znum0(i) = znum(labels0(i))
    znum1(i) = znum(labels1(i))
end do

if (remap) then

    ! Group atoms by label

    call grouplabels(natom, labels0, nblock0, blocksize0, blocktype0)
    call grouplabels(natom, labels1, nblock1, blocksize1, blocktype1)

    ! Contiguous same label order

    order0 = sortorder(blocktype0, natom)
    order1 = sortorder(blocktype1, natom)

    ! Undo contiguous same label order

    reorder0 = inversemap(order0)
    reorder1 = inversemap(order1)

    ! Abort if there are incompatible atomic symbols

    if (any(znum0(order0) /= znum1(order1))) then
        call error('The molecules are not isomers!')
    end if

    ! Abort if there are incompatible atomic labels

    if (any([(upper(labels0(order0(i))) /= upper(labels1(order1(i))), i=1, natom)])) then
        call error('There are incompatible atomic labels!')
    end if

else

    ! Abort if there are incompatible atomic symbols

    if (any(znum0 /= znum1)) then
        call error('The molecules are not isomers!')
    end if

    ! Abort if there are incompatible atomic labels

    if (any([(upper(labels0(i)) /= upper(labels1(i)), i=1, natom)])) then
        call error('There are incompatible atomic labels!')
    end if

end if

! Calculate normalized weights

select case (weighter)
case ('none')
    atomweight = [(1.0_wp, i=1, nelem)]
case ('mass')
    atomweight = stdmatom
case default
    call error('Invalid weighter option: '//trim(weighter))
end select

weights = atomweight(znum0)/sum(atomweight(znum0))

! Calculate atom sets centroids

center0 = centroid(natom, weights, atoms0)
center1 = centroid(natom, weights, atoms1)

if (remap) then

    ! Initialize random number generator

    call init_random_seed(seed)
    call random_seed(put=seed)

    ! Remap atoms to minimize distance and difference

    call remapatoms( &
        natom, nblock0, blocksize0, weights(order0), &
        translated(natom, center0, atoms0(:, order0)), &
        translated(natom, center1, atoms1(:, order1)), &
        maxrecord, nrecord, atomaplist, &
        trial_test, match_test &
    )

    ! Restore original atom ordering

    do i = 1, nrecord
        atomaplist(:, i) = atomaplist(reorder0, i)
    end do

else

    nrecord = 1
    atomaplist(:, 1) = [(i, i=1, natom)]

    call print_header()
    call print_stats( &
        0, 0, 0, 0.0_wp, 0.0_wp, 0.0_wp, &
        leastsquaredist( &
            natom, weights, &
            translated(natom, center0, atoms0), &
            translated(natom, center1, atoms1), &
            atomaplist(:, 1) &
        ) &
    )
    call print_footer(.false., .false., 0, 0)

end if

! Calculate optimal rotation matrices

do i = 1, nrecord
    rotmatlist(:, :, i) = quat2mat(leastrotquat(natom, weights, atoms0, atoms1, atomaplist(:, i)))
end do

end subroutine
