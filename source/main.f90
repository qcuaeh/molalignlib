program main

    use options
    use decoding
    use readwrite
    use utilities
    use rotation
    use ralign

    implicit none

    integer mapcount
    integer, dimension(:, :), allocatable :: atomaplist
    real, dimension(:, :, :), allocatable :: rotmatlist
    integer, dimension(:, :), allocatable :: bonds0, bonds1
    real, dimension(:, :), allocatable :: mol0, mol1
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

    call readmol('xyz', title0, label0, mol0, bonds0)
    call readmol('xyz', title1, label1, mol1, bonds1)

    ! Allocate records

    allocate(rotmatlist(3, 3, maxrecord))
    allocate(atomaplist(size(mol0, dim=2), maxrecord))

    ! Align atoms

    call alignatoms(mol0, mol1, label0, label1, bonds0, bonds1, mapcount, atomaplist, rotmatlist)

    ! Write aligned coordinates

    do i = 1, mapcount
        open (file_unit, file='aligned_'//str(i)//'.'//trim(output_format), action='write', status='replace')
        call writemol(file_unit, [(i, i=1, size(mol0, dim=2))], title0, label0, mol0, bonds0)
        call writemol(file_unit, atomaplist(:, i), title1, label1, &
            rotated(size(mol1, dim=2), rotmatlist(:, :, i), mol1), &
            bonds1)
        close (file_unit)
    end do

end program
