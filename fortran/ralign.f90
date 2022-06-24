program ralign

use options
use messages
use decoding
use readwrite
use strutils
use rotation
use superposition

implicit none

integer natom, nrecord, maxrecord
integer, dimension(:, :), allocatable :: atomaplist
real, dimension(:, :, :), allocatable :: rotmatlist
integer, dimension(:, :), allocatable :: bonds0, bonds1
real, dimension(:, :), allocatable :: atoms0, atoms1
character(32), dimension(:), allocatable :: label0, label1
character(512) title0, title1

integer i
logical argset
character(32) arg

! Set default options 

argset = .false.
biased = .false.
iterative = .false.
matching = .false.
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
    case ('-iter')
        iterative = .true.
    case ('-bias')
        biased = .true.
        call readoptarg(arg, tolerance)
    case ('-match')
        matching = .true.
    case ('-test')
        testing = .true.
    case ('-live')
        live = .true.
    case ('-scale')
        call readoptarg(arg, scaling)
    case ('-weight')
        call readoptarg(arg, weighting)
    case ('-record')
        call readoptarg(arg, maxrecord)
    case ('-out')
        call readoptarg(arg, output_format)
    case default
        argset = .true.
        call readarg(arg, maxcount)
    end select

end do

! Read coordinates

call readmol('xyz', title0, label0, atoms0, bonds0)
call readmol('xyz', title1, label1, atoms1, bonds1)

! Check number of atoms

if (size(label0) == size(label1)) then
    natom = size(label0)
else
    call error('The molecules do not have the same number of atoms!')
end if

! Allocate records

allocate(rotmatlist(3, 3, maxrecord))
allocate(atomaplist(natom, maxrecord))

! Superpose or align atoms

if (argset) then
    call superpose(natom, atoms0, atoms1, label0, label1, maxrecord, nrecord, atomaplist, rotmatlist)
else
    call align(natom, atoms0, atoms1, label0, label1, maxrecord, nrecord, atomaplist, rotmatlist)
end if

! Write aligned coordinates

do i = 1, nrecord
    open (file_unit, file='aligned_'//str(i)//'.'//trim(output_format), action='write', status='replace')
    call writemol(file_unit, [(i, i=1, natom)], title0, label0, atoms0, bonds0)
    call writemol(file_unit, atomaplist(:, i), title1, label1, rotated(natom, rotmatlist(:, :, i), atoms1), bonds1)
    close (file_unit)
end do

end program
