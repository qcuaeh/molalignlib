program ralign

use common
use options
use messages
use decoding
use readwrite
use strutils
use rotation
use superposition

implicit none

integer natom0, natom1
integer nrecord, maxrecord
integer, dimension(:, :), allocatable :: atomaplist
real(wp), dimension(:, :, :), allocatable :: rotmatlist
real(wp), dimension(:, :), allocatable :: atoms0, atoms1
character(32), dimension(:), allocatable :: labels0, labels1
character(512) title0, title1

integer i
logical remap
character(32) arg

! Set default options 

remap = .false.
biased = .false.
iterative = .false.
converge = .false.
testing = .false.
live = .false.
maxrecord = 9
scaling = 1000.
weighting = 'none'
output_format = 'xyz'

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
        converge = .false.
        call readoptarg(arg, maxcount)
    case ('-matches')
        remap = .true.
        converge = .true.
        call readoptarg(arg, maxcount)
    case ('-scale')
        call readoptarg(arg, scaling)
    case ('-weight')
        call readoptarg(arg, weighting)
    case ('-record')
        call readoptarg(arg, maxrecord)
    case ('-out')
        call readoptarg(arg, output_format)
    end select

end do

! Read coordinates

call readxyzfile(natom0, title0, labels0, atoms0)
call readxyzfile(natom1, title1, labels1, atoms1)

! Check number of atoms

if (natom0 /= natom1) then
    call error('The molecules do not have the same number of atoms!')
end if

! Allocate records

allocate(rotmatlist(3, 3, maxrecord))
allocate(atomaplist(natom0, maxrecord))

! Align or realign atoms

if (remap) then
    call realign(natom0, atoms0, atoms1, labels0, labels1, maxrecord, nrecord, atomaplist, rotmatlist)
else
    call align(natom0, atoms0, atoms1, labels0, labels1, maxrecord, nrecord, atomaplist, rotmatlist)
end if

! Write aligned coordinates

do i = 1, nrecord
    open (file_unit, file='aligned_'//str(i)//'.'//trim(output_format), action='write', status='replace')
    call writexyzfile(file_unit, natom0, [(i, i=1, natom0)], title0, labels0, atoms0)
    call writexyzfile(file_unit, natom1, atomaplist(:, i), title1, labels1, rotated(natom1, rotmatlist(:, :, i), atoms1))
    close (file_unit)
end do

end program
