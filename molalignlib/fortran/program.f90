program ralign

    use iso_fortran_env, only: input_unit
    use iso_fortran_env, only: output_unit

    use options
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
    logical sort_flag, stdin_flag
    character(arg_len) arg, format_w
    character(title_len) title0, title1
    integer, allocatable :: mapcount(:), maplist(:, :)
    integer, dimension(:), allocatable :: znums0, znums1, types0, types1
    character(label_len), dimension(:), allocatable :: labels0, labels1
    real, dimension(:, :), allocatable :: coords0, coords1
    real, allocatable, dimension(:) :: weights0
    real, allocatable :: mindist(:)
    real travec(3), rotmat(3, 3)

    procedure(f_realint), pointer :: weight_function => unity

    ! Set default options

    bias_flag = .false.
    conv_flag = .false.
    sort_flag = .false.
    abort_flag = .false.
    test_flag = .false.
    live_flag = .false.
    debug_flag = .false.
    stdin_flag = .false.

    records = 1
    maxcount = 10
    bias_tol = 0.35
    bias_scale = 1.e3
    format_w = 'xyz'

    ! Get user options

    call initarg()

    do while (getarg(arg))

        select case (arg)
        case ('-live')
            live_flag = .true.
        case ('-test')
            test_flag = .true.
        case ('-sort')
            sort_flag = .true.
        case ('-debug')
            debug_flag = .true.
        case ('-fast')
            bias_flag = .true.
            conv_flag = .true.
        case ('-mass')
            weight_function => stdmass
        case ('-count')
            call readoptarg(arg, maxcount)
        case ('-trials')
            abort_flag = .true.
            call readoptarg(arg, maxtrials)
        case ('-tol')
            call readoptarg(arg, bias_tol)
        case ('-scale')
            call readoptarg(arg, bias_scale)
        case ('-rec')
            call readoptarg(arg, records)
        case ('-out')
            call readoptarg(arg, format_w)
        case ('-stdin')
            stdin_flag = .true.
        case default
            call openfile(arg, file_unit)
        end select

    end do

    if (stdin_flag) then
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

    ! Read coordinates

    call readxyzfile(file_unit(1), natom0, title0, labels0, coords0)
    call readxyzfile(file_unit(2), natom1, title1, labels1, coords1)

    ! Allocate arrays

    allocate(znums0(natom0), znums1(natom1))
    allocate(types0(natom0), types1(natom1))
    allocate(weights0(natom0))
    allocate(maplist(natom0, records))
    allocate(mapcount(records), mindist(records))

    ! Get atomic numbers and types

    do i = 1, natom0
        call readlabel(labels0(i), znums0(i), types0(i))
    end do

    do i = 1, natom1
        call readlabel(labels1(i), znums1(i), types1(i))
    end do

    ! Set weights

    do i = 1, natom0
        weights0(i) = weight_function(znums0(i))
    end do

    ! Sort atoms to minimize MSD

    if (sort_flag) then

        call sort_atoms(natom0, natom1, znums0, znums1, types0, types1, coords0, coords1, &
            weights0, records, nrec, maplist, mapcount, mindist)

        ! Write aligned coordinates

        do i = 1, nrec

            call align_atoms(natom0, natom1, znums0, znums1(maplist(:, i)), &
                types0, types1(maplist(:, i)), coords0, coords1(:, maplist(:, i)), &
                weights0, travec, rotmat)

            open(newunit=u, file='aligned_'//str(i)//'.'//trim(format_w), action='write', status='replace')
            call writefile(u, format_w, natom0, title0, znums0, coords0)
            call writefile(u, format_w, natom1, title1, znums1(maplist(:, i)), &
                translated(natom1, rotated(natom1, coords1(:, maplist(:, i)), rotmat), travec))
            close(u)

        end do

    else

        ! Align atoms to minimize RMSD

        call align_atoms(natom0, natom1, znums0, znums1, types0, types1, coords0, coords1, &
            weights0, travec, rotmat)

        call rotate(natom1, coords1, rotmat)
        call translate(natom1, coords1, travec)

        ! Write aligned coordinates

        write (output_unit, '(a,1x,f0.4,1x,a)') 'RMSD:', &
            squaredist(natom0, weights0, coords0, coords1, identitymap(natom0)), &
            '(only alignment performed)'

        open(newunit=u, file='aligned_0.'//trim(format_w), action='write', status='replace')
        call writefile(u, format_w, natom0, title0, znums0, coords0)
        call writefile(u, format_w, natom1, title1, znums1, coords1) 
        close(u)

    end if

end program
