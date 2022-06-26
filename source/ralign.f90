program ralign

use options
use messages
use strutils
use readwrite
use decoding
use rotation
use translation
use routines

implicit none

integer natom0, natom1
integer nrecord, maxrecord
character(ttllen) title0, title1
real(wp), dimension(3) :: center0, center1
integer, dimension(:, :), allocatable :: atomaplist
real(wp), dimension(:, :, :), allocatable :: rotmatlist
character(lbllen), dimension(:), allocatable :: labels0, labels1
real(wp), dimension(:, :), allocatable :: atoms0, atoms1, atoms01

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

call readxyzfile(natom0, title0, labels0, atoms0)
call readxyzfile(natom1, title1, labels1, atoms1)

! Allocate arrays

allocate(atoms01(3, natom0))
allocate(atomaplist(natom0, maxrecord))
allocate(rotmatlist(3, 3, maxrecord))

! Superpose atoms

call superpose(natom0, natom1, labels0, labels1, atoms0, atoms1, &
    maxrecord, nrecord, center0, center1, atomaplist, rotmatlist)

! Write aligned coordinates

do i = 1, nrecord
    atoms01 = translated(natom1, -center0, rotated(natom1, rotmatlist(:, :, i), translated(natom1, center1, atoms1)))
    open (file_unit, file='aligned_'//str(i)//'.'//trim(outformat), action='write', status='replace')
    call writexyzfile(file_unit, natom0, title0, labels0, atoms0)
    call writexyzfile(file_unit, natom1, title1, labels1(atomaplist(:, i)), atoms01(:, atomaplist(:, i)))
    close (file_unit)
end do

end program
