module eigen
use stdio
use kinds

implicit none

private
public leasteigval
public leasteigvec

integer, parameter :: max_iter = 100  ! Maximum number of Jacobi iterations

contains

function leasteigval(A)
   real(rk), intent(inout) :: A(4,4)
   integer :: i, j, iter
   real(rk) :: leasteigval
   real(rk) :: V(4,4)
   real(rk) :: threshold, off_diag_norm

  ! Initialize V as identity matrix
   V = 0.0_rk
   do i = 1, 4
      V(i,i) = 1.0_rk
   end do

   do iter = 1, max_iter
      ! Compute off-diagonal norm
      off_diag_norm = 0.0_rk
      do i = 1, 3
         do j = i+1, 4
            off_diag_norm = off_diag_norm + A(i,j)**2
         end do
      end do
      off_diag_norm = sqrt(2.0_rk * off_diag_norm)

      if (off_diag_norm < epstol) exit

      threshold = off_diag_norm / (4.0_rk * sqrt(real(4, rk)))

      do i = 1, 3
         do j = i+1, 4
            if (abs(A(i,j)) > threshold) then
               call jacobi_rotation(A, V, i, j)
            end if
         end do
      end do
   end do

   ! Find the smallest eigenvalue
   leasteigval = A(1,1)
   do i = 2, 4
      if (A(i,i) < leasteigval) then
         leasteigval = A(i,i)
      end if
   end do
end function

function leasteigvec(A)
   real(rk), intent(inout) :: A(4,4)
   integer :: i, j, iter, leastindex
   real(rk) :: leasteigvec(4)
   real(rk) :: V(4,4)
   real(rk) :: threshold, off_diag_norm, leasteigval

   ! Initialize V as identity matrix
   V = 0.0_rk
   do i = 1, 4
      V(i,i) = 1.0_rk
   end do

   do iter = 1, max_iter
      ! Compute off-diagonal norm
      off_diag_norm = 0.0_rk
      do i = 1, 3
         do j = i+1, 4
            off_diag_norm = off_diag_norm + A(i,j)**2
         end do
      end do
      off_diag_norm = sqrt(2.0_rk * off_diag_norm)

      if (off_diag_norm < epstol) exit

      threshold = off_diag_norm / (4.0_rk * sqrt(real(4, rk)))

      do i = 1, 3
         do j = i+1, 4
            if (abs(A(i,j)) > threshold) then
               call jacobi_rotation(A, V, i, j)
            end if
         end do
      end do
   end do

   ! Find the index of the smallest eigenvalue
   leasteigval = A(1,1)
   leastindex = 1
   do i = 2, 4
      if (A(i,i) < leasteigval) then
         leasteigval = A(i,i)
         leastindex = i
      end if
   end do

   ! Extract the eigenvector corresponding to the smallest eigenvalue
   leasteigvec = V(:, leastindex)

   ! Normalize the eigenvector
   leasteigvec = leasteigvec / sqrt(sum(leasteigvec**2))
end function

subroutine jacobi_rotation(A, V, p, q)
   real(rk), intent(inout) :: A(4,4), V(4,4)
   integer, intent(in) :: p, q
   integer :: i
   real(rk) :: c, s, t, tau, temp

   if (abs(A(p,q)) > 1.0e-15_rk) then  ! Avoid division by zero
      tau = (A(q,q) - A(p,p)) / (2.0_rk * A(p,q))
      t = sign(1.0_rk, tau) / (abs(tau) + sqrt(1.0_rk + tau**2))
      c = 1.0_rk / sqrt(1.0_rk + t**2)
      s = t * c
   else
      c = 1.0_rk
      s = 0.0_rk
   end if

   temp = A(p,p)
   A(p,p) = c**2 * temp - 2.0_rk * c * s * A(p,q) + s**2 * A(q,q)
   A(q,q) = s**2 * temp + 2.0_rk * c * s * A(p,q) + c**2 * A(q,q)
   A(p,q) = 0.0_rk
   A(q,p) = 0.0_rk

   do i = 1, 4
      if (i /= p .and. i /= q) then
         temp = A(i,p)
         A(i,p) = c * temp - s * A(i,q)
         A(p,i) = A(i,p)
         A(i,q) = s * temp + c * A(i,q)
         A(q,i) = A(i,q)
      end if

      temp = V(i,p)
      V(i,p) = c * temp - s * V(i,q)
      V(i,q) = s * temp + c * V(i,q)
   end do
end subroutine

end module
