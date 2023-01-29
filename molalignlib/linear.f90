module linear
use parameters

implicit none

external ssyev
external dsyev

private
public det3
public syeval4
public syevec4
!public matmul

interface syeval4
   module procedure ssyeval4
   module procedure dsyeval4
end interface

interface syevec4
   module procedure ssyevec4
   module procedure dsyevec4
end interface

contains

real(wp) function det3(a) result(det)
   real(wp), dimension(3,3), intent(in) :: a

   det = a(1,1)*a(2,2)*a(3,3)  &
       - a(1,1)*a(2,3)*a(3,2)  &
       - a(1,2)*a(2,1)*a(3,3)  &
       + a(1,2)*a(2,3)*a(3,1)  &
       + a(1,3)*a(2,1)*a(3,2)  &
       - a(1,3)*a(2,2)*a(3,1)

end function

function matmul(a, b, m, o, n) result(ab)
   integer, intent(in) :: m, o ,n
   real(wp), intent (in) :: a(m, o)
   real(wp), intent (in) :: b(o, n)
   real(wp) :: ab(m, n)
   integer :: i, j

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
