module linear

implicit none
integer, parameter :: sp=4, dp=8

private
public trace
public det33
public matmul2
public syeval4
public syevec4

interface syeval4
    module procedure ssyeval4
    module procedure dsyeval4
end interface

interface syevec4
    module procedure ssyevec4
    module procedure dsyevec4
end interface

contains

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

subroutine ssyeval4(a, w)
    real(sp), dimension(:, :), intent(in) :: a
    real(sp), dimension(:), intent(out) :: w
    integer info
    real(sp) work(20)
    call ssyev('N', 'U', 4, a, size(a, dim=1), w, work, 20, info)
end subroutine

subroutine dsyeval4(a, w)
    real(dp), dimension(:, :), intent(in) :: a
    real(dp), dimension(:), intent(out) :: w
    integer info
    real(dp) work(20)
    call dsyev('N', 'U', 4, a, size(a, dim=1), w, work, 20, info)
end subroutine

subroutine ssyevec4(a, w)
    real(sp), dimension(:, :), intent(inout) :: a
    real(sp), dimension(:), intent(out) :: w
    integer info
    real(sp) work(20)
    call ssyev('V', 'U', 4, a, size(a, dim=1), w, work, 20, info)
end subroutine

subroutine dsyevec4(a, w)
    real(dp), dimension(:, :), intent(inout) :: a
    real(dp), dimension(:), intent(out) :: w
    integer info
    real(dp) work(20)
    call dsyev('V', 'U', 4, a, size(a, dim=1), w, work, 20, info)
end subroutine

end module
