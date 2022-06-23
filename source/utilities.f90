module utilities
use iso_fortran_env, only: error_unit
use sorting
implicit none

interface str
    module procedure int2str
    module procedure real2str
end interface

contains

function inversemap(mapping)
    integer, intent(in), dimension(:) :: mapping
    integer, dimension(size(mapping)) :: inversemap

    integer i

    do i = 1, size(mapping)
        inversemap(mapping(i)) = i
    end do

end function

logical function is_badmap(mapping, n)
    integer, intent(in) :: n, mapping(:)
    integer i
    is_badmap = .false.
    if (any(sorted(mapping, n) /= [(i, i=1, n)])) then
        is_badmap = .true.
    end if
end function

subroutine union(n1, list1, n2, list2, n3, list3)
! Find the union of two sorted lists without repeated elements

    integer, intent(in) :: n1, n2, list1(:), list2(:)
    integer, intent(out) :: n3, list3(:)

    integer i1, i2

    i1 = 1
    i2 = 1
    n3 = 0

    do while (i1 <= n1 .and. i2 <= n2)
        n3 = n3 + 1
        if (list1(i1) < list2(i2)) then
            list3(n3) = list1(i1)
            i1 = i1 + 1
        else if (list1(i1) > list2(i2)) then
            list3(n3) = list2(i2)
            i2 = i2 + 1
        else
            list3(n3) = list2(i2)
            i1 = i1 + 1
            i2 = i2 + 1
        end if
    end do

    do while (i1 <= n1)
        n3 = n3 + 1
        list3(n3) = list1(i1)
        i1 = i1 + 1
    end do

    do while (i2 <= n2)
        n3 = n3 + 1
        list3(n3) = list2(i2)
        i2 = i2 + 1
    end do

end subroutine

subroutine intersection(n1, list1, n2, list2, n3, list3)
! Find the intersection of two sorted lists without repeated elements

    integer, intent(in) :: n1, n2, list1(:), list2(:)
    integer, intent(out) :: n3, list3(:)

    integer i1, i2

    i1 = 1
    i2 = 1
    n3 = 0

    do while (i1 <= n1 .and. i2 <= n2)
        if (list1(i1) < list2(i2)) then
            i1 = i1 + 1
        else if (list1(i1) > list2(i2)) then
            i2 = i2 + 1
        else
            n3 = n3 + 1
            list3(n3) = list2(i2)
            i1 = i1 + 1
            i2 = i2 + 1
        end if
    end do

end subroutine

function matmul2(a, b, m, o, n) result(ab)
    integer, intent(in) :: m, o ,n
    real, intent (in) :: a(m, o)
    real, intent (in) :: b(o, n)
    real ab(m, n)
    integer i, j

    do i = 1, n
        ab(:, i) = 0.0
        do j = 1, o
            ab(:, i) = ab(:, i) + a(:, j)*b(j, i)
        end do
    end do

end function

function trace(natom, matrix)
integer, intent(in) :: natom
real, dimension(:, :), intent(in) :: matrix
real trace

integer i

trace = 0

do i = 1, natom
    trace = trace + matrix(i, i)
end do

end function

function det33(a) result(det)
real, dimension(3,3), intent(in)  :: a

real det

det = a(1,1)*a(2,2)*a(3,3)  &
    - a(1,1)*a(2,3)*a(3,2)  &
    - a(1,2)*a(2,1)*a(3,3)  &
    + a(1,2)*a(2,3)*a(3,1)  &
    + a(1,3)*a(2,1)*a(3,2)  &
    - a(1,3)*a(2,2)*a(3,1)

end function

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

subroutine split_name(filename, prefix, suffix)
    character(*), intent(in) :: filename
    character(*), intent(out) :: prefix, suffix
    integer pos

    pos = scan(trim(filename), '.', BACK=.true.)
    if (pos <= 1 .or. pos >= len_trim(filename)) then
        write (error_unit, '(a, x, a, x, a)') trim(filename), 'is not a valid file name!'
        stop
    end if
    prefix = filename(:pos-1)
    suffix = filename(pos+1:)

end subroutine

! Returns the position in vecInd where the element ind is found or 0 otherwise
function findInd(ind, nInd, vecInd) result(pos)
    integer, intent(in) :: ind, nInd, vecInd(nInd)
    integer :: pos
    integer :: i
    pos = 0
    do i = 1, nInd
        if (ind == vecInd(i)) then
            pos = i
            return
        end if
    end do
end function findInd

! Adds a new element to a list of indices in growing order, no repeats
subroutine addInd(ind, nInd, vecInd)
    integer, intent(in) :: ind
    integer, intent(inout) :: nInd, vecInd(:)
    integer :: p, e
    p = 1
    do while (p <= nInd)   ! find position for growing order
        if (ind == vecInd(p)) return
        if (ind < vecInd(p)) then
            exit
        else
            p = p + 1
        end if
    end do
    nInd = nInd + 1
    e = nInd
    do while (e > p)
        vecInd(e) = vecInd(e - 1)   ! shifts one position from the end
        e = e - 1
    end do
    vecInd(p) = ind
end subroutine addInd

! Deletes element ind from vecInd if present.
subroutine delInd(ind,nInd,vecInd)
    integer, intent(in) :: ind
    integer, intent(inout) :: nind, vecInd(:)
    integer :: pos, i
    pos = findInd(ind,nInd,vecInd(:nInd))
    if (pos /= 0) then
        do i = pos + 1, nInd
            vecInd(i-1) = vecInd(i)
        end do
        nInd = nInd - 1
    end if
end subroutine delInd

! deletes from vec all indices already tracked (value .true. in track)
subroutine delTrack (track, nvec, vec)
    logical, dimension(:), intent(in) :: track
    integer, intent(inout) :: nvec, vec(:)
    integer i
    i = 1
    do while ( i <= nvec )
        if ( track(vec(i)) ) then
            call delInd (vec(i), nvec, vec)
        else
            i = i + 1
        end if
    end do
end subroutine delTrack

! calculates the dihedral angle between 4 atom indices
subroutine calc_dihedral (ind_add, vec4ind, atoms, calc, dihed)
    integer, intent(in) :: ind_add
    integer, intent(inout) :: vec4ind(4)
    real, dimension(:, :), intent(in) :: atoms
    logical, intent(out) :: calc
    real, intent(out) :: dihed
    real :: x, y

    real, dimension(3) :: r1, r2, r3, r4, n1, n2, m
    integer i

    ! adds last index if not equal to 0, shifts to the left the last 3 indices
    if ( ind_add /= 0 ) then
        do i = 1, 3
            vec4ind(i) = vec4ind(i+1)
        end do
        vec4ind(4) = ind_add
    end if

    ! returns if there are not 4 indices specified
    do i = 1, 4
        if ( vec4ind(i) == 0 ) then
            calc = .false.
            return
        else
            calc = .true.   ! means that the calculation can actually be performed
        end if
    end do

    ! vectors between atoms
    r1(:) =  atoms(:,vec4ind(2)) - atoms(:,vec4ind(1))
    r2(:) = -atoms(:,vec4ind(3)) + atoms(:,vec4ind(2))
    r3(:) =  atoms(:,vec4ind(4)) - atoms(:,vec4ind(3))

    ! normalize
    r1(:) = r1(:)/sqrt(dot_product(r1, r1))
    r2(:) = r2(:)/sqrt(dot_product(r2, r2))
    r3(:) = r3(:)/sqrt(dot_product(r3, r3))

    ! normal planes
    n1 = cross_prod (r1, r2)
    n2 = cross_prod (r2, r3)
    m  = cross_prod (n1, r2)

    ! coordinates of n2*
    x = dot_product (n1, n2)
    y = dot_product (m,  n2)

    ! angle between normals using arctan function
    dihed = 180.0/3.14159265*atan2(y, x)

!    print '(a,i0,3(x,i0),a,f0.2)', "Dihed(", vec4ind(:), ") = ", dihed

contains

    function cross_prod (v1, v2) result (cp)
        real, dimension(3), intent(in) :: v1, v2
        real, dimension(3) :: cp
        cp(1) =  v1(2)*v2(3) - v2(2)*v1(3)
        cp(2) = -v1(1)*v2(3) + v2(1)*v1(3)
        cp(3) =  v1(1)*v2(2) - v2(1)*v1(2)
    end function cross_prod

end subroutine calc_dihedral

end module
