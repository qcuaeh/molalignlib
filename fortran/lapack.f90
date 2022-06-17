module lapack
! Convenience interfaces to lapack subroutines
implicit none
private
public syeval3
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

interface syeval3
    module procedure ssyeval3
    module procedure dsyeval3
end interface

contains

subroutine ssyeval4(a, w)
real(4), dimension(:, :), intent(in) :: a
real(4), dimension(:), intent(out) :: w
integer info
real work(20)

call ssyev('N', 'U', 4, a, size(a, dim=1), w, work, 20, info)

end subroutine

subroutine dsyeval4(a, w)
real(8), dimension(:, :), intent(in) :: a
real(8), dimension(:), intent(out) :: w
integer info
real work(20)

call dsyev('N', 'U', 4, a, size(a, dim=1), w, work, 20, info)

end subroutine

subroutine ssyevec4(a, w)
real(4), dimension(:, :), intent(inout) :: a
real(4), dimension(:), intent(out) :: w
integer info
real work(20)

call ssyev('V', 'U', 4, a, size(a, dim=1), w, work, 20, info)

end subroutine

subroutine dsyevec4(a, w)
real(8), dimension(:, :), intent(inout) :: a
real(8), dimension(:), intent(out) :: w
integer info
real work(20)

call dsyev('V', 'U', 4, a, size(a, dim=1), w, work, 20, info)

end subroutine

subroutine ssyeval3(a, w)
real(4), dimension(:, :), intent(in) :: a
real(4), dimension(:), intent(out) :: w
integer info
real work(12)

call ssyev('N', 'U', 3, a, size(a, dim=1), w, work, 12, info)

end subroutine

subroutine dsyeval3(a, w)
real(8), dimension(:, :), intent(in) :: a
real(8), dimension(:), intent(out) :: w
integer info
real work(12)

!call dsyevc3(a, w)
call dsyev('N', 'U', 3, a, size(a, dim=1), w, work, 12, info)

end subroutine

end module
