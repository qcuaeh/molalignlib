module strutils

use iso_fortran_env, only: error_unit

use options

implicit none

private

public str
public lower
public upper
public splitext

interface str
    module procedure int2str
    module procedure real2str
end interface

contains

function lower(str)
    character(*), intent(in) :: str
    character(26), parameter :: lowercase = 'abcdefghijklmnopqrstuvwxyz'
    character(26), parameter :: uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len(str)) :: lower

    integer :: i, j

    lower = str

    do j = 1, len(str)
        i = index(uppercase, str(j:j))
        if (i > 0) lower(j:j) = lowercase(i:i)
    end do
end function

function upper(str)
    character(*), intent(in) :: str
    character(26), parameter :: lowercase = 'abcdefghijklmnopqrstuvwxyz'
    character(26), parameter :: uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len(str)) :: upper

    integer :: i, j

    upper = str

    do j = 1, len(str)
        i = index(lowercase, str(j:j))
        if (i > 0) upper(j:j) = uppercase(i:i)
    end do
end function

function int2str(x) result(strx)
    integer, intent(in) :: x
    character(floor(log10(real(x))) + 1) strx
    write (strx, '(i0)') x
end function

function real2str(x) result(strx)
    real, intent(in) :: x
    character(floor(log10(x)) + 6) strx
    write (strx, '(f0.4)') x
end function

subroutine splitext(filename, prefix, suffix)
    character(*), intent(in) :: filename
    character(*), intent(out) :: prefix, suffix
    integer pos

    pos = scan(trim(filename), '.', BACK=.true.)
    if (pos <= 1 .or. pos >= len_trim(filename)) then
        write (error_unit, '(a,1x,a,1x,a)') trim(filename), 'is not a valid file name!'
        stop
    end if
    prefix = filename(:pos-1)
    suffix = filename(pos+1:)
end subroutine

end module
