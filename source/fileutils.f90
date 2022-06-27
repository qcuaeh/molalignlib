module fileutils

use iso_fortran_env, only: error_unit

use options

implicit none

private

public opened
public open_unit

contains

logical function opened(unit)

    integer, intent(in) :: unit
    inquire(unit, opened=opened)

end function

subroutine open_unit(path)

    character(optlen), intent(in) :: path

    if (opened(first_file_unit)) then
        if (opened(second_file_unit)) then
            write (error_unit, '(a)') 'Too many paths'
            stop
        else
            open(second_file_unit, file=path, action='read')
        end if
    else
        open(first_file_unit, file=path, action='read')
    end if

end subroutine

end module
