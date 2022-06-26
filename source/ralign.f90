program ralign

use options
use messages
use strutils
use readwrite
use decoding
use rotation
use translation
!use superposition

implicit none

integer natom0, natom1
integer nrecord, maxrecord
character(ttllen) title0, title1
real(wp), dimension(3) :: center0, center1
integer, dimension(:, :), allocatable :: atomaplist
real(wp), dimension(:, :, :), allocatable :: rotmatlist
character(lbllen), dimension(:), allocatable :: labels0, labels1
integer, dimension(:), allocatable :: znums0, znums1, types0, types1
real(wp), dimension(:, :), allocatable :: coords0, coords1, auxcoords

integer i
character(optlen) arg

! Set default options

live = .false.
remap = .false.
biased = .false.
iterative = .false.
trialing = .false.
matching = .false.
testing = .false.
maxrecord = 9
biasscale = 1000.0_wp
weighter = 'none'
outformat = 'xyz'

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
        remap = .true.
        trialing = .true.
        call readoptarg(arg, maxtrial)
    case ('-matches')
        remap = .true.
        matching = .true.
        call readoptarg(arg, maxmatch)
    case ('-bias-scale')
        call readoptarg(arg, biasscale)
    case ('-weight')
        call readoptarg(arg, weighter)
    case ('-records')
        call readoptarg(arg, maxrecord)
    case ('-out')
        call readoptarg(arg, outformat)
    case default
        if (arg(:1) == '-') then
            write (error_unit, '(a, x, a)') 'Unknown option:', trim(arg)
            stop
        else
            call error('Invalid argument: '//trim(arg))
        end if
    end select

end do

! Read coordinates

call readxyzfile(natom0, title0, labels0, coords0)
call readxyzfile(natom1, title1, labels1, coords1)

! Allocate arrays

allocate(znums0(natom0), znums1(natom1))
allocate(types0(natom0), types1(natom1))
allocate(atomaplist(natom0, maxrecord))
allocate(rotmatlist(3, 3, maxrecord))

! Get atomic numbers and types

do i = 1, natom0
    call getznum(labels0(i), znums0(i), types0(i))
end do

do i = 1, natom1
    call getznum(labels1(i), znums1(i), types1(i))
end do

! Superpose atoms

call superpose(natom0, natom1, znums0, znums1, types0, types1, coords0, coords1, &
    maxrecord, nrecord, center0, center1, atomaplist, rotmatlist)

! Write aligned coordinates

allocate(auxcoords(3, natom0))

do i = 1, nrecord
    auxcoords = translated(natom1, -center0, rotated(natom1, rotmatlist(:, :, i), translated(natom1, center1, coords1)))
    open (file_unit, file='aligned_'//str(i)//'.'//trim(outformat), action='write', status='replace')
    call writexyzfile(file_unit, natom0, title0, znums0, coords0)
    call writexyzfile(file_unit, natom1, title1, znums1(atomaplist(:, i)), auxcoords(:, atomaplist(:, i)))
    close (file_unit)
end do

end program
