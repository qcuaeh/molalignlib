module decoding

use iso_fortran_env, only: error_unit

use options

implicit none

integer iarg

private
public initarg
public getarg
public readoptarg

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
    character(optlen) arg

    if (iarg <= command_argument_count()) then
        call get_command_argument(iarg, arg)
        iarg = iarg + 1
        getarg = .true.
    else
        getarg = .false.
    end if

end function

subroutine getoptarg(option, optarg)
    character(optlen), intent(in) :: option
    character(optlen), intent(out) :: optarg

    if (iarg <= command_argument_count()) then
        call get_command_argument(iarg, optarg)
        if (optarg(:1) /= '-') then
            iarg = iarg + 1
            return
        end if
    end if
    write (error_unit, '(a, x, a, x, a)') 'Option', trim(option), 'requires an argument'
    stop

end subroutine

subroutine readintoptarg(option, optval)
    character(optlen), intent(in) :: option
    integer, intent(out) :: optval
    character(optlen) optarg
    integer stat

    call getoptarg(option, optarg)
    read (optarg, *, iostat=stat) optval
    if (stat /= 0) then
        write (error_unit, '(a, x, a, x, a)') 'Option', trim(option), 'requires an integer argument'
        stop
    end if

end subroutine

subroutine readrealoptarg(option, optval)
    character(optlen), intent(in) :: option
    real(wp), intent(out) :: optval
    character(optlen) optarg
    integer stat

    call getoptarg(option, optarg)
    read (optarg, *, iostat=stat) optval
    if (stat /= 0) then
        write (error_unit, '(a, x, a, x, a)') 'Option', trim(option), 'requires a numeric argument'
        stop
    end if

end subroutine

end module
