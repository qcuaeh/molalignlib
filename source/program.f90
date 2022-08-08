program ralign

    use iso_fortran_env, only: input_unit
    use iso_fortran_env, only: output_unit

    use options
    use chemdata
    use strutils
    use optparse
    use chemutils
    use readwrite
    use translation
    use rotation
    use alignment
    use library

    implicit none

    integer i, u
    integer file_unit(2)
    integer nrec, records
    integer natom0, natom1
    character(arg_len) arg, path
    character(title_len) title0, title1
    integer, allocatable :: mapcount(:), maplist(:, :)
    integer, dimension(:), allocatable :: znums0, znums1, types0, types1
    character(label_len), dimension(:), allocatable :: labels0, labels1
    real, dimension(:, :), allocatable :: coords0, coords1
    real, allocatable :: weights0(:), mindist(:)
    real travec(3), rotmat(3, 3)
    real property(nelem)
    logical stdin

    procedure(writeabstractfile), pointer :: writefile => null()

    ! Set default options

    live = .false.
    stdin = .false.
    biased = .false.
    iteration = .false.
    complete = .false.
    converge = .false.
    testing = .false.
    records = 1
    lenscale = 1000.0
    property = [(1.0, i=1, nelem)]
    formatout = 'xyz'

    ! Get user options

    call initarg()

    do while (getarg(arg))

        select case (arg)
        case ('-live')
            live = .true.
        case ('-test')
            testing = .true.
        case ('-iter')
            iteration = .true.
        case ('-mass')
            property = stdmatom
        case ('-bias')
            biased = .true.
            call readoptarg(arg, tolerance)
        case ('-max')
            complete = .true.
            call readoptarg(arg, maxtrial)
        case ('-count')
            converge = .true.
            call readoptarg(arg, mincount)
        case ('-scale')
            call readoptarg(arg, lenscale)
        case ('-rec')
            call readoptarg(arg, records)
        case ('-out')
            call readoptarg(arg, formatout)
        case ('-stdin')
            stdin = .true.
        case default
            call openfile(arg, file_unit)
        end select

    end do

    if (stdin) then
        file_unit(1) = input_unit
        file_unit(2) = input_unit
    else
        if (ifile == 0) then
            write (error_unit, '(a)') 'Error: No file paths were specified'
            stop
        else if (ifile == 1) then
            file_unit(2) = file_unit(1)
        end if
    end if

    ! Select output format

    select case (formatout)
    case ('xyz')
        writefile => writexyzfile
    case ('mol2')
        writefile => writemol2file
    case default
        write (error_unit, '(a,x,a)') 'Invalid format:', trim(formatout)
        stop
    end select

    ! Read coordinates

    call readxyzfile(file_unit(1), natom0, title0, labels0, coords0)
    call readxyzfile(file_unit(2), natom1, title1, labels1, coords1)

    ! Allocate arrays

    allocate(znums0(natom0), znums1(natom1))
    allocate(types0(natom0), types1(natom1))
    allocate(maplist(natom0, records), mapcount(records), mindist(records))
    allocate(weights0(natom0))

    ! Get atomic numbers and types

    do i = 1, natom0
        call getznum(labels0(i), znums0(i), types0(i))
    end do

    do i = 1, natom1
        call getznum(labels1(i), znums1(i), types1(i))
    end do

    ! Get normalized weights

    weights0 = property(znums0)/sum(property(znums0))

    ! Superpose atoms

    if (complete .or. converge) then

        call remap(natom0, natom1, znums0, znums1, types0, types1, &
            coords0, coords1, weights0, records, nrec, maplist, mapcount, mindist)

        ! Write aligned coordinates

        do i = 1, nrec

            call align(natom0, natom1, znums0, znums1(maplist(:, i)), types0, &
                types1(maplist(:, i)), coords0, coords1(:, maplist(:, i)), &
                weights0, travec, rotmat)

            open(newunit=u, file='aligned_'//str(i)//'.'//trim(formatout), action='write', status='replace')
            call writefile(u, natom0, title0, znums0, coords0)
            call writefile(u, natom1, title1, znums1(maplist(:, i)), &
                translated(natom1, rotated(natom1, coords1(:, maplist(:, i)), rotmat), travec))
            close(u)

        end do

    else

        call align(natom0, natom1, znums0, znums1, types0, types1, &
            coords0, coords1, weights0, travec, rotmat)

        call rotate(natom1, coords1, rotmat)
        call translate(natom1, coords1, travec)

        write (output_unit, '(a,x,f0.4,x,a)') 'RMSD:', &
            squaredist(natom0, weights0, coords0, coords1, identitymap(natom0)), &
            '(only alignment performed)'

        open(newunit=u, file='aligned_1.'//trim(formatout), action='write', status='replace')
        call writefile(u, natom0, title0, znums0, coords0)
        call writefile(u, natom1, title1, znums1, coords1) 
        close(u)

    end if

end program
