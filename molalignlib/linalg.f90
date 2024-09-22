module linalg
use stdio
use kinds

implicit none

contains

real(rk) function matdet3(a)
   real(rk), intent(in) :: a(3, 3)

   matdet3 = a(1,1)*a(2,2)*a(3,3)  &
       - a(1,1)*a(2,3)*a(3,2)  &
       - a(1,2)*a(2,1)*a(3,3)  &
       + a(1,2)*a(2,3)*a(3,1)  &
       + a(1,3)*a(2,1)*a(3,2)  &
       - a(1,3)*a(2,2)*a(3,1)

end function

function matinv3(a)
   real(rk), intent(in) :: a(3, 3)
   real(rk) :: matinv3(3, 3)
   real(rk) :: matdet, invdet

   matdet = matdet3(a)

   if (abs(matdet) < epstol) then
      write (stderr, '(a)') 'Error: Matrix is not invertible'
      stop
   end if

   invdet = 1.0_rk / matdet

   matinv3(1,1) = (a(2,2) * a(3,3) - a(2,3) * a(3,2)) * invdet
   matinv3(1,2) = (a(1,3) * a(3,2) - a(1,2) * a(3,3)) * invdet
   matinv3(1,3) = (a(1,2) * a(2,3) - a(1,3) * a(2,2)) * invdet
   matinv3(2,1) = (a(2,3) * a(3,1) - a(2,1) * a(3,3)) * invdet
   matinv3(2,2) = (a(1,1) * a(3,3) - a(1,3) * a(3,1)) * invdet
   matinv3(2,3) = (a(1,3) * a(2,1) - a(1,1) * a(2,3)) * invdet
   matinv3(3,1) = (a(2,1) * a(3,2) - a(2,2) * a(3,1)) * invdet
   matinv3(3,2) = (a(1,2) * a(3,1) - a(1,1) * a(3,2)) * invdet
   matinv3(3,3) = (a(1,1) * a(2,2) - a(1,2) * a(2,1)) * invdet

end function

end module
