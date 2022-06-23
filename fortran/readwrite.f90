module readwrite
use iso_fortran_env, only: input_unit
use iso_fortran_env, only: output_unit
use iso_fortran_env, only: error_unit
use options
use strutils
use maputils
use chemistry
implicit none
private
public readmol
public writemol
public file_unit

integer, parameter :: file_unit = 999

contains

subroutine readmol(mol_format, title, label, atoms, bonds)
    character(*), intent(in) :: mol_format
    integer, dimension(:, :), allocatable, intent(out) :: bonds
    real, dimension(:, :), allocatable, intent(out) :: atoms
    character(*), intent(out) :: title
    character(*), dimension(:), allocatable, intent(out) :: label

    integer natom, nbond
    integer i

    select case (mol_format)
    case ('xyz')
        call readxyzfile(natom, title, label, atoms)
        allocate (bonds(2, maxcoord*natom))
    case ('mol')
        call readmolfile(natom, title, label, atoms, nbond, bonds)
    case default
        write (error_unit, '(a)') 'Unknown file format'
        stop
    end select

end subroutine

subroutine readxyzfile(natom, title, label, atoms)
    integer, intent(out) :: natom
    real, dimension(:, :), allocatable, intent(out) :: atoms
    character(*), dimension(:), allocatable, intent(out) :: label
    character(*), intent(out) :: title
    integer i, stat

    read (input_unit, *, iostat=stat) natom
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    allocate (label(natom), atoms(3, natom))

    read (input_unit, '(a)', iostat=stat) title
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    do i = 1, natom
        read (input_unit, *, iostat=stat) label(i), atoms(:, i)
        if (stat < 0) then
            write (error_unit, '(a)') 'Unexpected end of file!'
            stop
        end if
    end do

end subroutine

subroutine readmolfile(natom, title, label, atoms, nbond, bonds)
    integer, intent(out) :: natom, nbond
    integer, dimension(:, :), allocatable, intent(out) :: bonds
    real, dimension(:, :), allocatable, intent(out) :: atoms
    character(*), dimension(:), allocatable, intent(out) :: label
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

    allocate (label(natom), atoms(3, natom), bonds(2, maxcoord*natom))

    do i = 1, natom
        read (input_unit, *, iostat=stat) atoms(:, i), label(i)
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

subroutine writemol(file_unit, mapping, title, label, atoms, bonds)
    integer, intent(in) :: file_unit
    character(*), intent(in) :: title
    character(*), dimension(:), intent(in) :: label
    integer, dimension(:), intent(in) :: mapping
    integer, dimension(:, :), intent(in) :: bonds
    real, dimension(:, :), intent(in) :: atoms

    select case (output_format)
    case ('xyz')
        call writexyzfile(file_unit, mapping, title, label, atoms)
    case ('mol2')
        call writemol2file(file_unit, mapping, title, label, atoms, bonds)
    case default
        write (error_unit, '(a)') 'Unknown file format'
        stop
    end select

end subroutine

subroutine writexyzfile(file_unit, mapping, title, label, atoms)
    integer, intent(in) :: file_unit
    integer, dimension(:), intent(in) :: mapping
    real, dimension(:, :), intent(in) :: atoms
    character(*), dimension(:), intent(in) :: label
    character(*), intent(in) :: title

    integer i

    write (file_unit, '(i0)') size(atoms, dim=2)
    write (file_unit, '(a)') trim(title)
    do i = 1, size(atoms, dim=2)
        write (file_unit, '(a, 3(2x, f12.6))') elsym(znum(label(mapping(i)))), atoms(:, mapping(i))
    end do

end subroutine

subroutine writemol2file(file_unit, mapping, title, label, atoms, bonds)
    integer, intent(in) :: file_unit
    integer, dimension(:), intent(in) :: mapping
    integer, dimension(:, :), intent(in) :: bonds
    real, dimension(:, :), intent(in) :: atoms
    character(*), dimension(:), intent(in) :: label
    character(*), intent(in) :: title

    integer i
    integer unmapping(size(mapping))

    unmapping = inversemap(mapping)

    write (file_unit, '(a)') '@<TRIPOS>MOLECULE'
    write (file_unit, '(a)') trim(title)
    write (file_unit, '(5(i4, x))') size(atoms, dim=2), size(bonds, dim=2), 0, 0, 0
    write (file_unit, '(a)') 'SMALL'
    write (file_unit, '(a)') 'NO_CHARGES'
    write (file_unit, '(a)') '@<TRIPOS>ATOM'

    do i = 1, size(atoms, dim=2)
        write (file_unit, '(i4, x, a2, 3(x, f12.6), x, a8, x, i2, x, a4, x, f7.3)') &
            i, elsym(znum(label(mapping(i)))), atoms(:, mapping(i)), label(mapping(i)), 1, 'MOL1', 0.
    end do

    write (file_unit, '(a)') '@<TRIPOS>BOND'

    do i = 1, size(bonds, dim=2)
        write (file_unit, '(i4, x, 2(x, i4), x, a2)') i, unmapping(bonds(:, i)), '1'
    end do

end subroutine

end module
