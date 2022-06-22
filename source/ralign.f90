program ralign

    use globals
    use decoding
    use printing
    use remapping
    use biasing

    implicit none

    character(32) arg

    iterative = .false.
    backtrack = .false.
    permutate = .false.
    topology = .false.
    debug = .false.
    live = .false.

    maxcoord = 0
    maxcount = -1
    maxrecord = 9

    scaling = 1000.
    tolerance = 0.1

    weighting = 'none'
    output_file = 'aligned.xyz'

    bias_test => bias_none
    stop_test => trial_stop_test

    call initarg()

    do while (getarg(arg))

        select case (arg)
        case ('-iter')
            iterative = .true.
        case ('-converge')
            stop_test => match_stop_test
        case ('-bias1')
            bias_test => bias_test1
        case ('-bias2')
            bias_test => bias_test2
        case ('-debug')
            debug = .true.
        case ('-live')
            live = .true.
        case ('-scale')
            call readoptarg(arg, scaling)
        case ('-tol')
            call readoptarg(arg, tolerance)
        case ('-weight')
            call readoptarg(arg, weighting)
        case ('-out')
            call readoptarg(arg, output_file)
        case ('-n')
            call readoptarg(arg, maxrecord)
        case default
            call readarg(arg, maxcount)
        end select

    end do

    call main()

end program
