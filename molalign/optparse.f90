module optparse
use parameters

implicit none

integer :: iarg, ipos

private
public ipos
public initarg
public getarg
public readarg
public readoptarg

interface readoptarg
    module procedure readstroptarg
    module procedure readintoptarg
    module procedure readrealoptarg
end interface

contains

subroutine initarg()

    iarg = 0
    ipos = 0

end subroutine

subroutine readarg(arg, files)

    character(arg_len), intent(in) :: arg
    character(arg_len), intent(out) :: files(:)

    if (arg(1:1) == '-') then
        write (error_unit, '(a,1x,a)') 'Error: Unknown option', trim(arg)
        stop
    end if

    ipos = ipos + 1

    if (ipos > size(files)) then
        write (error_unit, '(a)') 'Error: Too many files'
        stop
    end if

    files(ipos) = arg

end subroutine

logical function getarg(arg)

    character(arg_len), intent(out) :: arg

    iarg = iarg + 1

    if (iarg <= command_argument_count()) then
        call get_command_argument(iarg, arg)
        getarg = .true.
    else
        getarg = .false.
    end if

end function

subroutine getoptarg(option, optarg)

    character(arg_len), intent(in) :: option
    character(arg_len), intent(out) :: optarg

    iarg = iarg + 1

    if (iarg <= command_argument_count()) then
        call get_command_argument(iarg, optarg)
        if (optarg(1:1) /= '-') then
            return
        end if
    end if

    write (error_unit, '(a,1x,a,1x,a)') 'Option', trim(option), 'requires an argument'
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
    character(arg_len) :: optarg
    integer :: stat

    call getoptarg(option, optarg)
    read (optarg, *, iostat=stat) optval
    if (stat /= 0) then
        write (error_unit, '(a,1x,a,1x,a)') 'Option', trim(option), 'requires an integer argument'
        stop
    end if

end subroutine

subroutine readrealoptarg(option, optval)

    character(arg_len), intent(in) :: option
    real(wp), intent(out) :: optval
    character(arg_len) :: optarg
    integer :: stat

    call getoptarg(option, optarg)
    read (optarg, *, iostat=stat) optval
    if (stat /= 0) then
        write (error_unit, '(a,1x,a,1x,a)') 'Option', trim(option), 'requires a numeric argument'
        stop
    end if

end subroutine

end module
