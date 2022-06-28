module chemutils

use iso_fortran_env, only: error_unit

use chemdata
use strutils

private

public getznum

contains

subroutine getznum(label, znum, type)

    character(*), intent(in) :: label
    integer, intent(out) :: znum, type
    integer m, n, stat

    n = len_trim(label)
    m = verify(upper(trim(label)), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')

    if (m == 0) then
        type = 1
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

end module
