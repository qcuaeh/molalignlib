module decoding

use iso_fortran_env, only: error_unit

use options
use chemdata
use strutils
use fileutils

implicit none

integer iarg
integer first_unit
integer second_unit

private

public getarg
public openfile
public readoptarg
public initarg
public getznum
public first_unit
public second_unit

interface readoptarg
    module procedure getoptarg
    module procedure readintoptarg
    module procedure readrealoptarg
end interface

contains

subroutine getznum(label, znum, type)
    character(*), intent(in) :: label
    integer, intent(out) :: znum, type
    integer m, n, stat

    n = len_trim(label)
    m = verify(upper(trim(label)), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')

    if (m == 0) then
        type = 0
        znum = atomic_number(label)
    else
        read (label(m:n), *, iostat=stat) type
        if (stat /= 0) then
            write (error_unit, '(2a)') 'Invalid label: ', trim(label)
            stop
        end if
        znum = atomic_number(label(1:m-1))
    end if

end subroutine

function atomic_number(symbol) result(z)
    character(*), intent(in) :: symbol
    integer z

    do z = 1, nelem
        if (upper(symbol) == upper(elsym(z))) then
            return
        end if
    end do

    write (error_unit, '(a)') 'Unknown atomic symbol: ', trim(symbol)

end function

subroutine initarg()
    iarg = 1
    first_unit = first_file_unit
    second_unit = second_file_unit
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

subroutine openfile(arg)
    character(optlen), intent(in) :: arg

    if (arg(1:1) == '-') then
        write (error_unit, '(a, x, a)') 'Unknown option:', trim(arg)
        stop
    end if

    if (opened(first_file_unit)) then
        if (opened(second_file_unit)) then
            write (error_unit, '(a)') 'Too many paths'
            stop
        else
            open(second_file_unit, file=arg, action='read')
        end if
    else
        open(first_file_unit, file=arg, action='read')
    end if

end subroutine

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
