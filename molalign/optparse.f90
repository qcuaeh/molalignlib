module optparse
use parameters

implicit none

integer iarg, ifile

private

public ifile
public initarg
public getarg
public openfile
public readoptarg

interface readoptarg
    module procedure readstroptarg
    module procedure readintoptarg
    module procedure readrealoptarg
end interface

contains

subroutine initarg()

    iarg = 0
    ifile = 0

end subroutine

logical function getarg(arg)

    character(arg_len) arg

    iarg = iarg + 1

    if (iarg <= command_argument_count()) then
        call get_command_argument(iarg, arg)
        getarg = .true.
    else
        getarg = .false.
    end if

end function

subroutine openfile(arg, file_unit)

    character(arg_len) arg
    integer file_unit(2), stat

    if (arg(1:1) == '-') then
        write (error_unit, '(a,1x,a)') 'Error: Unknown option', trim(arg)
        stop
    end if

    ifile = ifile + 1

    if (ifile <= 2) then
        open(newunit=file_unit(ifile), file=arg, action='read', iostat=stat)
        if (stat /= 0) then
            write (error_unit, '(a,1x,a,1x,a)') 'Error:', trim(arg), 'does not exist'
            stop
        end if
    else
        write (error_unit, '(a)') 'Error: Too many files'
        stop
    end if

end subroutine

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
    character(arg_len) optarg
    integer stat

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
    character(arg_len) optarg
    integer stat

    call getoptarg(option, optarg)
    read (optarg, *, iostat=stat) optval
    if (stat /= 0) then
        write (error_unit, '(a,1x,a,1x,a)') 'Option', trim(option), 'requires a numeric argument'
        stop
    end if

end subroutine

end module
