module readwrite

use iso_fortran_env, only: error_unit

use options
use strutils
use maputils
use chemdata

implicit none

abstract interface
    subroutine writeabstractfile(file_unit, natom, title, znums, coords, nbond, bonds)
        use options
        character(*), intent(in) :: title
        integer, intent(in) :: file_unit, natom
        integer, dimension(:), intent(in) :: znums
        real(wp), dimension(:, :), intent(in) :: coords
        integer, optional, intent(in) :: nbond
        integer, target, optional, intent(in) :: bonds(:, :)
    end subroutine
end interface

contains

subroutine readxyzfile(file_unit, natom, title, labels, coords)

    integer, intent(in) :: file_unit
    integer, intent(out) :: natom
    real(wp), dimension(:, :), allocatable, intent(out) :: coords
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


subroutine writexyzfile(file_unit, natom, title, znums, coords)

    character(*), intent(in) :: title
    integer, intent(in) :: file_unit, natom
    integer, dimension(:), intent(in) :: znums
    real(wp), dimension(:, :), intent(in) :: coords

    integer i

    write (file_unit, '(i0)') natom
    write (file_unit, '(a)') trim(title)
    do i = 1, natom
        write (file_unit, '(a, 3(2x, f12.6))') elsym(znums(i)), coords(:, i)
    end do

end subroutine

subroutine writemol2file(file_unit, natom, title, znums, coords, opt_nbond, opt_bonds)

    character(*), intent(in) :: title
    integer, intent(in) :: file_unit, natom
    integer, dimension(:), intent(in) :: znums
    real(wp), dimension(:, :), intent(in) :: coords
    integer, optional, intent(in) :: opt_nbond
    integer, target, optional, intent(in) :: opt_bonds(:, :)

    integer i, nbond
    integer, pointer :: bonds(:, :)
    

    if (present(opt_nbond)) then
        nbond = opt_nbond
        if (present(opt_bonds)) then
            bonds => opt_bonds
        else
            write (error_unit, '(a)') 'Bond list is missing!'
        end if
    else
        nbond = 0
    end if

    write (file_unit, '(a)') '@<TRIPOS>MOLECULE'
    write (file_unit, '(a)') trim(title)
    write (file_unit, '(5(i4, x))') natom, nbond, 0, 0, 0
    write (file_unit, '(a)') 'SMALL'
    write (file_unit, '(a)') 'NO_CHARGES'
    write (file_unit, '(a)') '@<TRIPOS>ATOM'

    do i = 1, natom
        write (file_unit, '(i4, x, a2, 3(x, f12.6), x, a8, x, i2, x, a4, x, f7.3)') &
            i, elsym(znums(i)), coords(:, i), elsym(znums(i)), 1, 'MOL1', 0.
    end do

    write (file_unit, '(a)') '@<TRIPOS>BOND'

    do i = 1, nbond
        write (file_unit, '(i4, x, 2(x, i4), x, a2)') i, bonds(:, i), '1'
    end do

end subroutine

end module
