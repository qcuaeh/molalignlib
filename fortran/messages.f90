module messages
use iso_fortran_env, only: output_unit
use iso_fortran_env, only: error_unit
implicit none
private
public error
public warning
public spacer

interface warning
   module procedure warning
   module procedure warning_from_caller
end interface

contains

subroutine error(message)
    character(*), intent(in) :: message
    write (error_unit, '("Error:",x,a)') message
    stop
end subroutine

subroutine warning(message)
    character(*), intent(in) :: message
    write (error_unit, '("Warning:",x,a)') message
end subroutine

subroutine warning_from_caller(message, caller)
    character(*), intent(in) :: message, caller
    write (error_unit, '("Warning from",x,a, ":",x,a)') caller, message
end subroutine

subroutine spacer(i)
    integer, intent(in) :: i
    write (error_unit, '(i0,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^",i0)') i, i
end subroutine

end module

