program ralign

use options
use messages
use decoding
use readwrite
use utilities
use rotation

implicit none

integer natom, mapcount
integer, dimension(:, :), allocatable :: atomaplist
real, dimension(:, :, :), allocatable :: rotmatlist
integer, dimension(:, :), allocatable :: bonds0, bonds1
real, dimension(:, :), allocatable :: atoms0, atoms1
character(32), dimension(:), allocatable :: label0, label1
character(512) title0, title1

integer i
character(32) arg

! Set default options 

biased = .false.
iterative = .false.
testing = .false.
live = .false.

maxcoord = 0
maxcount = -1
maxrecord = 9

scaling = 1000.
tolerance = 0.1

weighting = 'none'
output_format = 'xyz'

converge = .false.

! Get user options 

call initarg()

do while (getarg(arg))

    select case (arg)
    case ('-iter')
        iterative = .true.
    case ('-converge')
        converge = .true.
    case ('-bias')
        biased = .true.
    case ('-testing')
        testing = .true.
    case ('-live')
        live = .true.
    case ('-scale')
        call readoptarg(arg, scaling)
    case ('-tol')
        call readoptarg(arg, tolerance)
    case ('-weight')
        call readoptarg(arg, weighting)
    case ('-out')
        call readoptarg(arg, output_format)
    case ('-n')
        call readoptarg(arg, maxrecord)
    case default
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

! Align atoms

call superpose(natom, atoms0, atoms1, label0, label1, mapcount, atomaplist, rotmatlist)

! Write aligned coordinates

do i = 1, mapcount
    open (file_unit, file='aligned_'//str(i)//'.'//trim(output_format), action='write', status='replace')
    call writemol(file_unit, [(i, i=1, natom)], title0, label0, atoms0, bonds0)
    call writemol(file_unit, atomaplist(:, i), title1, label1, rotated(natom, rotmatlist(:, :, i), atoms1), bonds1)
    close (file_unit)
end do

end program
