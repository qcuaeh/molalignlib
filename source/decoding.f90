module decoding
use iso_fortran_env, only: error_unit
implicit none
integer iarg

private
public getarg
public initarg
public readarg
public readoptarg

interface readarg
    module procedure readintarg
end interface

interface readoptarg
    module procedure getoptarg
    module procedure readintoptarg
    module procedure readrealoptarg
end interface

contains

subroutine initarg()
    iarg = 1
end subroutine

logical function getarg(arg)
    character(32) arg
    if (iarg <= command_argument_count()) then
        call get_command_argument(iarg, arg)
        iarg = iarg + 1
        getarg = .true.
    else
        getarg = .false.
    end if
end function

subroutine readintarg(arg, argval)
    character(32), intent(in) :: arg
    integer, intent(out) :: argval
    integer stat

    if (arg(:1) == '-') then
        write (error_unit, '(a, x, a)') 'Unknown arg:', trim(arg)
        stop
    end if
    read (arg, *, iostat=stat) argval
    if (stat /= 0) then
        write (error_unit, '(a, x, a)') 'Argument must be an integer:', trim(arg)
        stop
    end if

end subroutine

subroutine getoptarg(arg, optarg)
    character(32), intent(in) :: arg
    character(32), intent(out) :: optarg

    if (iarg <= command_argument_count()) then
        call get_command_argument(iarg, optarg)
        if (optarg(:1) /= '-') then
            iarg = iarg + 1
            return
        end if
    end if
    write (error_unit, '(a, x, a, x, a)') 'Option', trim(arg), 'requires an argument'
    stop

end subroutine

subroutine readintoptarg(arg, optval)
    character(32), intent(in) :: arg
    integer, intent(out) :: optval
    character(32) optarg
    integer stat

    call getoptarg(arg, optarg)
    read (optarg, *, iostat=stat) optval
    if (stat /= 0) then
        write (error_unit, '(a, x, a, x, a)') 'Option', trim(arg), 'requires an integer argument'
        stop
    end if

end subroutine

subroutine readrealoptarg(arg, optval)
    character(32), intent(in) :: arg
    real, intent(out) :: optval
    character(32) optarg
    integer stat

    call getoptarg(arg, optarg)
    read (optarg, *, iostat=stat) optval
    if (stat /= 0) then
        write (error_unit, '(a, x, a, x, a)') 'Option', trim(arg), 'requires a numeric argument'
        stop
    end if

end subroutine

end module
