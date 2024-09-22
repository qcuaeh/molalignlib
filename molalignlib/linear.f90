module linear
use stdio
use kinds

implicit none

external ssyev
external dsyev

private
public det3
public inv3
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

!function matmul(a, b, m, o, n) result(ab)
!   integer, intent(in) :: m, o ,n
!   real(wp), intent (in) :: a(m, o)
!   real(wp), intent (in) :: b(o, n)
!   real(wp) :: ab(m, n)
!   integer :: i, j
!
!   do i = 1, n
!      ab(:, i) = 0.0
!      do j = 1, o
!         ab(:, i) = ab(:, i) + a(:, j)*b(j, i)
!      end do
!   end do
!
!end function

real(wp) function det3(a) result(det)
   real(wp), dimension(3,3), intent(in) :: a

   det = a(1,1)*a(2,2)*a(3,3)  &
       - a(1,1)*a(2,3)*a(3,2)  &
       - a(1,2)*a(2,1)*a(3,3)  &
       + a(1,2)*a(2,3)*a(3,1)  &
       + a(1,3)*a(2,1)*a(3,2)  &
       - a(1,3)*a(2,2)*a(3,1)

end function

function inv3(a) result(inv)
   real(wp), intent(in) :: a(3,3)
   real(wp) :: inv(3,3)
   real(wp) :: det

   det = det3(a)

   if (abs(det) < 1.0e-6) then
      write (stderr, '(a)') 'Error: Matrix is not invertible'
      stop
   end if

   inv(1,1) = (a(2,2) * a(3,3) - a(2,3) * a(3,2)) / det
   inv(1,2) = (a(1,3) * a(3,2) - a(1,2) * a(3,3)) / det
   inv(1,3) = (a(1,2) * a(2,3) - a(1,3) * a(2,2)) / det
   inv(2,1) = (a(2,3) * a(3,1) - a(2,1) * a(3,3)) / det
   inv(2,2) = (a(1,1) * a(3,3) - a(1,3) * a(3,1)) / det
   inv(2,3) = (a(1,3) * a(2,1) - a(1,1) * a(2,3)) / det
   inv(3,1) = (a(2,1) * a(3,2) - a(2,2) * a(3,1)) / det
   inv(3,2) = (a(1,2) * a(3,1) - a(1,1) * a(3,2)) / det
   inv(3,3) = (a(1,1) * a(2,2) - a(1,2) * a(2,1)) / det

end function

subroutine ssyeval4(a, w)
   real(sp), intent(in) :: a(4, 4)
   real(sp), intent(out) :: w(4)
   integer :: info
   real(sp) work(25)
   call ssyev('N', 'U', 4, a, 4, w, work, 25, info)
end subroutine

subroutine dsyeval4(a, w)
   real(dp), intent(in) :: a(4, 4)
   real(dp), intent(out) :: w(4)
   integer :: info
   real(dp) work(25)
   call dsyev('N', 'U', 4, a, 4, w, work, 25, info)
end subroutine

subroutine ssyevec4(a, w)
   real(sp), intent(inout) :: a(4, 4)
   real(sp), intent(out) :: w(4)
   integer :: info
   real(sp) work(25)
   call ssyev('V', 'U', 4, a, 4, w, work, 25, info)
end subroutine

subroutine dsyevec4(a, w)
   real(dp), intent(inout) :: a(4, 4)
   real(dp), intent(out) :: w(4)
   integer :: info
   real(dp) work(25)
   call dsyev('V', 'U', 4, a, 4, w, work, 25, info)
end subroutine

end module
