module readwrite

use iso_fortran_env, only: error_unit

use options
use strutils
use maputils
use chemdata

implicit none

contains

subroutine writefile(file_unit, format_w, natom, title, znums, coords, opt_nbond, opt_bonds)

    integer, intent(in) :: file_unit, natom
    integer, dimension(:), intent(in) :: znums
    integer, optional, intent(in) :: opt_nbond
    integer, target, optional, intent(in) :: opt_bonds(:, :)
    real, dimension(:, :), intent(in) :: coords
    character(*), intent(in) :: title, format_w

    integer nbond
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

    select case (format_w)
    case ('xyz')
        call writexyzfile(file_unit, natom, title, znums, coords)
    case ('mol2')
        call writemol2file(file_unit, natom, title, znums, coords, nbond, bonds)
    case default
        write (error_unit, '(a,1x,a)') 'Invalid format:', trim(format_w)
        stop
    end select

end subroutine

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


subroutine writexyzfile(file_unit, natom, title, znums, coords)

    character(*), intent(in) :: title
    integer, intent(in) :: file_unit, natom
    integer, dimension(:), intent(in) :: znums
    real, dimension(:, :), intent(in) :: coords

    integer i

    write (file_unit, '(i0)') natom
    write (file_unit, '(a)') trim(title)
    do i = 1, natom
        write (file_unit, '(a,3(2x,f12.6))') elsym(znums(i)), coords(:, i)
    end do

end subroutine

subroutine writemol2file(file_unit, natom, title, znums, coords, nbond, bonds)

    character(*), intent(in) :: title
    integer, intent(in) :: file_unit, natom
    integer, dimension(:), intent(in) :: znums
    real, dimension(:, :), intent(in) :: coords
    integer, intent(in) :: nbond
    integer, target, intent(in) :: bonds(:, :)

    integer i
    
    write (file_unit, '(a)') '@<TRIPOS>MOLECULE'
    if (title /= '') then
        write (file_unit, '(a)') trim(title)
    else
        write (file_unit, '(a)') 'Untitled'
    end if
    write (file_unit, '(5(i4,1x))') natom, nbond, 0, 0, 0
    write (file_unit, '(a)') 'SMALL'
    write (file_unit, '(a)') 'NO_CHARGES'
    write (file_unit, '(a)') '@<TRIPOS>ATOM'

    do i = 1, natom
        write (file_unit, '(i4,1x,a2,3(1x,f12.6),1x,a8,1x,i2,1x,a4,1x,f7.3)') &
            i, elsym(znums(i)), coords(:, i), elsym(znums(i)), 1, 'MOL1', 0.
    end do

    write (file_unit, '(a)') '@<TRIPOS>BOND'

    do i = 1, nbond
        write (file_unit, '(i4,1x,2(1x,i4),1x,a2)') i, bonds(:, i), '1'
    end do

end subroutine

end module
