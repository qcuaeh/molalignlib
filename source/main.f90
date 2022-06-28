program ralign

use iso_fortran_env, only: input_unit
use iso_fortran_env, only: output_unit

use options
use chemdata
use strutils
use argutils
use fileutils
use chemutils
use readwrite
use translation
use rotation
use alignment
!use superposition

implicit none

integer i
integer natom0, natom1
integer nrecord, maxrecord
character(ttllen) title0, title1
real(wp) travec(3), rotmat(3, 3)
integer, allocatable :: countlist(:)
integer, allocatable :: atomaplist(:, :)
character(lbllen), dimension(:), allocatable :: labels0, labels1
integer, dimension(:), allocatable :: znums0, znums1, types0, types1
real(wp), dimension(:, :), allocatable :: coords0, coords1
real(wp), dimension(:), allocatable :: weights0, weights1
real(wp), dimension(nelem) :: property
integer first_unit, second_unit
character(optlen) arg, path

! Set default options

live = .false.
biased = .false.
iterative = .false.
trialing = .false.
counting = .false.
testing = .false.
maxrecord = 9
biasscale = 1000.0_wp
weighter = 'none'
outformat = 'xyz'

! Set defualt file units

first_unit = first_file_unit
second_unit = second_file_unit

! Get user options

call initarg()

do while (getarg(arg))

    select case (arg)
    case ('-live')
        live = .true.
    case ('-test')
        testing = .true.
    case ('-iter')
        iterative = .true.
    case ('-bias')
        biased = .true.
        call readoptarg(arg, tolerance)
    case ('-trials')
        trialing = .true.
        call readoptarg(arg, maxtrial)
    case ('-count')
        counting = .true.
        call readoptarg(arg, maxcount)
    case ('-bias-scale')
        call readoptarg(arg, biasscale)
    case ('-weight')
        call readoptarg(arg, weighter)
    case ('-maps')
        call readoptarg(arg, maxrecord)
    case ('-out')
        call readoptarg(arg, outformat)
    case ('-stdin')
        first_unit = input_unit
        second_unit = input_unit
    case default
        call readarg(arg, path)
        call open_unit(path)
    end select

end do

if (.not. opened(first_unit)) then
    write (error_unit, '(a)') 'No paths were specified'
    stop
else if (.not. opened(second_unit)) then
    second_unit = first_unit
end if

! Read coordinates

call readxyzfile(first_unit, natom0, title0, labels0, coords0)
call readxyzfile(second_unit, natom1, title1, labels1, coords1)

! Allocate arrays

allocate(znums0(natom0), znums1(natom1))
allocate(types0(natom0), types1(natom1))
allocate(weights0(natom0), weights1(natom1))
allocate(atomaplist(natom0, maxrecord), countlist(maxrecord))

! Get atomic numbers and types

do i = 1, natom0
    call getznum(labels0(i), znums0(i), types0(i))
end do

do i = 1, natom1
    call getznum(labels1(i), znums1(i), types1(i))
end do

! Select weighting property

select case (weighter)
case ('none')
    property = [(1.0_wp, i=1, nelem)]
case ('mass')
    property = stdmatom
case default
    write (error_unit, '(a,x,a)') 'Invalid weighter option:', trim(weighter)
    stop
end select

! Get normalized weights

weights0 = property(znums0)/sum(property(znums0))
weights1 = property(znums1)/sum(property(znums1))

! Superpose atoms

if (trialing .or. counting) then

    call remap(natom0, natom1, znums0, znums1, types0, types1, weights0, weights1, &
        coords0, coords1, maxrecord, nrecord, atomaplist, countlist)

    ! Write aligned coordinates

    do i = 1, nrecord

        call align( &
            natom0, natom1, &
            znums0, znums1(atomaplist(:, i)), &
            types0, types1(atomaplist(:, i)), &
            weights0, weights1(atomaplist(:, i)), &
            coords0, coords1(:, atomaplist(:, i)), &
            travec, rotmat &
        )

        open(temp_file_unit, file='aligned_'//str(i)//'.'//trim(outformat), action='write', status='replace')
        call writexyzfile(temp_file_unit, natom0, title0, znums0, coords0)
        call writexyzfile( &
            temp_file_unit, natom1, title1, znums1(atomaplist(:, i)), &
            translated(natom1, rotated(natom1, coords1(:, atomaplist(:, i)), rotmat), travec) &
        )
        close(temp_file_unit)

    end do

else

    call align(natom0, natom1, znums0, znums1, types0, types1, weights0, weights1, &
        coords0, coords1, travec, rotmat)

    call rotate(natom1, coords1, rotmat)
    call translate(natom1, coords1, travec)

    write (output_unit, '(a,x,f0.4,x,a)') 'RMSD:', &
        squaredist(natom0, weights0, coords0, coords1, identitymap(natom0)), &
        '(only alignment performed)'

    open(temp_file_unit, file='aligned.'//trim(outformat), action='write', status='replace')
    call writexyzfile(temp_file_unit, natom0, title0, znums0, coords0)
    call writexyzfile(temp_file_unit, natom1, title1, znums1, coords1) 
    close(temp_file_unit)

end if

end program
