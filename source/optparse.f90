module optparse

use iso_fortran_env, only: error_unit

use options

implicit none

integer iarg

private

public initarg
public getarg
public readarg
public readoptarg

interface readarg
    module procedure readstrarg
end interface

interface readoptarg
    module procedure readstroptarg
    module procedure readintoptarg
    module procedure readrealoptarg
end interface

contains

subroutine initarg()

    iarg = 1

end subroutine

logical function getarg(arg)

    character(arg_len) arg

    if (iarg <= command_argument_count()) then
        call get_command_argument(iarg, arg)
        iarg = iarg + 1
        getarg = .true.
    else
        getarg = .false.
    end if

end function

subroutine readstrarg(arg, argval)

    character(arg_len), intent(in) :: arg
    character(arg_len), intent(out) :: argval

    if (arg(1:1) == '-') then
        write (error_unit, '(a, x, a)') 'Unknown option:', trim(arg)
        stop
    end if

    argval = arg

end subroutine

subroutine getoptarg(option, optarg)

    character(arg_len), intent(in) :: option
    character(arg_len), intent(out) :: optarg

    if (iarg <= command_argument_count()) then
        call get_command_argument(iarg, optarg)
        if (optarg(1:1) /= '-') then
            iarg = iarg + 1
            return
        end if
    end if
    write (error_unit, '(a, x, a, x, a)') 'Option', trim(option), 'requires an argument'
    stop

end subroutine

subroutine readstroptarg(option, optval)

    character(arg_len), intent(in) :: option
    character(arg_len), intent(out) :: optval

    call getoptarg(option, optval)

end subroutine

subroutine readintoptarg(option, optval)

    character(arg_len), intent(in) :: option
    integer, intent(out) :: optval
    character(arg_len) optarg
    integer stat

    call getoptarg(option, optarg)
    read (optarg, *, iostat=stat) optval
    if (stat /= 0) then
        write (error_unit, '(a, x, a, x, a)') 'Option', trim(option), 'requires an integer argument'
        stop
    end if

end subroutine

subroutine readrealoptarg(option, optval)

    character(arg_len), intent(in) :: option
    real, intent(out) :: optval
    character(arg_len) optarg
    integer stat

    call getoptarg(option, optarg)
    read (optarg, *, iostat=stat) optval
    if (stat /= 0) then
        write (error_unit, '(a, x, a, x, a)') 'Option', trim(option), 'requires a numeric argument'
        stop
    end if

end subroutine

end module
