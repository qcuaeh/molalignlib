module readwrite

use iso_fortran_env, only: input_unit
use iso_fortran_env, only: output_unit
use iso_fortran_env, only: error_unit

use options
use strutils
use maputils
use chemdata

implicit none

integer, parameter :: file_unit = 999

contains


subroutine readxyzfile(natom, title, labels, atoms)
    integer, intent(out) :: natom
    real(wp), dimension(:, :), allocatable, intent(out) :: atoms
    character(*), dimension(:), allocatable, intent(out) :: labels
    character(*), intent(out) :: title
    integer i, stat

    read (input_unit, *, iostat=stat) natom
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    allocate (labels(natom), atoms(3, natom))

    read (input_unit, '(a)', iostat=stat) title
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    do i = 1, natom
        read (input_unit, *, iostat=stat) labels(i), atoms(:, i)
        if (stat < 0) then
            write (error_unit, '(a)') 'Unexpected end of file!'
            stop
        end if
    end do

end subroutine


subroutine readmolfile(natom, title, labels, atoms, nbond, bonds)
    integer, intent(out) :: natom, nbond
    integer, dimension(:, :), allocatable, intent(out) :: bonds
    real(wp), dimension(:, :), allocatable, intent(out) :: atoms
    character(*), dimension(:), allocatable, intent(out) :: labels
    character(*), intent(out) :: title
    integer i, stat

    read (input_unit, '(a)', iostat=stat) title
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    read (input_unit, *, iostat=stat)
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    read (input_unit, *, iostat=stat)
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    read (input_unit, *, iostat=stat) natom, nbond
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    allocate (labels(natom), atoms(3, natom), bonds(2, maxcoord*natom))

    do i = 1, natom
        read (input_unit, *, iostat=stat) atoms(:, i), labels(i)
        if (stat < 0) then
            write (error_unit, '(a)') 'Unexpected end of file!'
            stop
        end if
    end do

    do i = 1, nbond
        read (input_unit, *, iostat=stat) bonds(:, i)
        if (stat < 0) then
            write (error_unit, '(a)') 'Unexpected end of file!'
            stop
        end if
    end do

    read (input_unit, *, iostat=stat)
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    read (input_unit, *, iostat=stat)
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    read (input_unit, *, iostat=stat)
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

end subroutine


subroutine writexyzfile(file_unit, natom, title, znums, atoms)
    character(*), intent(in) :: title
    integer, intent(in) :: file_unit, natom
    integer, dimension(:), intent(in) :: znums
    real(wp), dimension(:, :), intent(in) :: atoms

    integer i

    write (file_unit, '(i0)') natom
    write (file_unit, '(a)') trim(title)
    do i = 1, natom
        write (file_unit, '(a, 3(2x, f12.6))') elsym(znums(i)), atoms(:, i)
    end do

end subroutine


subroutine writemol2file(file_unit, natom, nbond, title, znums, atoms, bonds)
    character(*), intent(in) :: title
    integer, intent(in) :: file_unit, natom
    integer, dimension(:), intent(in) :: znums
    real(wp), dimension(:, :), intent(in) :: atoms
    integer, intent(in) :: nbond
    integer, dimension(:, :), intent(in) :: bonds

    integer i

    write (file_unit, '(a)') '@<TRIPOS>MOLECULE'
    write (file_unit, '(a)') trim(title)
    write (file_unit, '(5(i4, x))') natom, nbond, 0, 0, 0
    write (file_unit, '(a)') 'SMALL'
    write (file_unit, '(a)') 'NO_CHARGES'
    write (file_unit, '(a)') '@<TRIPOS>ATOM'

    do i = 1, natom
        write (file_unit, '(i4, x, a2, 3(x, f12.6), x, a8, x, i2, x, a4, x, f7.3)') &
            i, elsym(znums(i)), atoms(:, i), elsym(znums(i)), 1, 'MOL1', 0.
    end do

    write (file_unit, '(a)') '@<TRIPOS>BOND'

    do i = 1, nbond
        write (file_unit, '(i4, x, 2(x, i4), x, a2)') i, bonds(:, i), '1'
    end do

end subroutine

end module
