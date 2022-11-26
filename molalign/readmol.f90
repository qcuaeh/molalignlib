module readmol

use iso_fortran_env, only: error_unit

use settings
use strutils
use maputils
use chemdata

implicit none

contains

subroutine readxyzfile(file_unit, natom, title, labels, coords)

    integer, intent(in) :: file_unit
    integer, intent(out) :: natom
    real, dimension(:, :), allocatable, intent(out) :: coords
    character(*), dimension(:), allocatable, intent(out) :: labels
    character(*), intent(out) :: title
    integer i, stat

    read (file_unit, *, iostat=stat) natom
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    allocate (labels(natom), coords(3, natom))

    read (file_unit, '(a)', iostat=stat) title
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    do i = 1, natom
        read (file_unit, *, iostat=stat) labels(i), coords(:, i)
        if (stat < 0) then
            write (error_unit, '(a)') 'Unexpected end of file!'
            stop
        end if
    end do

end subroutine

end module
