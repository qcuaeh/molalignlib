module lapack
use parameters

implicit none

external ssyev
external dsyev

private
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

subroutine ssyeval4(a, w)
   real(sp), dimension(:, :), intent(in) :: a
   real(sp), dimension(:), intent(out) :: w
   integer :: info
   real(sp) work(20)
   call ssyev('N', 'U', 4, a, size(a, dim=1), w, work, 20, info)
end subroutine

subroutine dsyeval4(a, w)
   real(dp), dimension(:, :), intent(in) :: a
   real(dp), dimension(:), intent(out) :: w
   integer :: info
   real(dp) work(20)
   call dsyev('N', 'U', 4, a, size(a, dim=1), w, work, 20, info)
end subroutine

subroutine ssyevec4(a, w)
   real(sp), dimension(:, :), intent(inout) :: a
   real(sp), dimension(:), intent(out) :: w
   integer :: info
   real(sp) work(20)
   call ssyev('V', 'U', 4, a, size(a, dim=1), w, work, 20, info)
end subroutine

subroutine dsyevec4(a, w)
   real(dp), dimension(:, :), intent(inout) :: a
   real(dp), dimension(:), intent(out) :: w
   integer :: info
   real(dp) work(20)
   call dsyev('V', 'U', 4, a, size(a, dim=1), w, work, 20, info)
end subroutine

end module
