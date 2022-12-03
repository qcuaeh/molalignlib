module readmol
use parameters
use strutils
use chemdata

implicit none

contains

subroutine open2read(filepath, unit, fileext)
    character(*), intent(in) :: filepath
    integer, intent(out) :: unit
    character(*), intent(out) :: fileext
    character(len(filepath)) :: filename
    integer :: stat
    open(newunit=unit, file=filepath, action='read', status='old', iostat=stat)
    if (stat /= 0) then
        write (error_unit, '(a,1x,a,1x,a)') 'Error opening', trim(filepath), 'for reading'
        stop
    end if
    filename = basename(filepath)
    if (len_trim(filename) == 0) then
        write (error_unit, '(a,1x,a)') 'Invalid file path:', trim(filepath)
        stop
    end if
    fileext = getext(filename)
    if (len_trim(fileext) == 0) then
        write (error_unit, '(a,1x,a)') 'Missing file extension:', trim(filepath)
        stop
    end if
end subroutine

subroutine readfile(unit, fmtin, natom, title, labels, coords, opt_nbond, opt_bonds)

    integer, intent(in) :: unit
    character(*), intent(in) :: fmtin
    integer, intent(out) :: natom
    character(*), intent(out) :: title
    character(*), dimension(:), allocatable, intent(out) :: labels
    real(wp), dimension(:, :), allocatable, intent(out) :: coords
    integer, optional, intent(out) :: opt_nbond
    integer, target, optional, intent(out) :: opt_bonds(:, :)

    integer :: nbond
    integer, pointer :: bonds(:, :)
    
    if (present(opt_nbond)) then
        nbond = opt_nbond
        if (present(opt_bonds)) then
            bonds => opt_bonds
        else
            write (error_unit, '(a)') 'Bond list is missing!'
            stop
        end if
    else
        nbond = 0
    end if

    select case (fmtin)
    case ('xyz')
        call readxyzfile(unit, natom, title, labels, coords)
!    case ('mol2')
!        call readmol2file(unit, natom, title, labels, coords, nbond, bonds)
    case default
        write (error_unit, '(a,1x,a)') 'Invalid format:', trim(fmtin)
        stop
    end select

end subroutine

subroutine readxyzfile(unit, natom, title, labels, coords)

    integer, intent(in) :: unit
    integer, intent(out) :: natom
    real(wp), dimension(:, :), allocatable, intent(out) :: coords
    character(*), dimension(:), allocatable, intent(out) :: labels
    character(*), intent(out) :: title
    integer :: i, stat

    read (unit, *, iostat=stat) natom
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    allocate (labels(natom), coords(3, natom))

    read (unit, '(a)', iostat=stat) title
    if (stat < 0) then
        write (error_unit, '(a)') 'Unexpected end of file!'
        stop
    end if

    do i = 1, natom
        read (unit, *, iostat=stat) labels(i), coords(:, i)
        if (stat < 0) then
            write (error_unit, '(a)') 'Unexpected end of file!'
            stop
        end if
    end do

end subroutine

end module
