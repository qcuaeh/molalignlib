! MolAlignLib
! Copyright (C) 2022 José M. Vásquez

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module alignment
use lapack
use rotation
use settings

implicit none

private
public squaredist
public leastsquaredist
public leastrotquat

contains

real(wp) function squaredist(natom, weights, coords0, coords1, atomap) result(dist2)
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: atomap
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords0, coords1

   dist2 = sum(weights(1:natom)*sum((coords0(:, 1:natom) - coords1(:, atomap(1:natom)))**2, dim=1))

end function

! Purpose: Calculate least square distance from eigenvalues
real(wp) function leastsquaredist(natom, weights, coords0, coords1, atomap) result(dist2)
    integer, intent(in) :: natom
    integer, dimension(:), intent(in) :: atomap
    real(wp), dimension(:), intent(in) :: weights
    real(wp), dimension(:, :), intent(in) :: coords0, coords1
    real(wp) eigmat(4, 4), eigval(4)

    call kearsleymat(natom, weights, coords0, coords1, atomap, eigmat)
    call syeval4(eigmat, eigval)
    dist2 = max(eigval(1), 0._wp)

end function

!! Purpose: Calculate least square distance from aligned coordinates
!real(wp) function leastsquaredist(natom, weights, coords0, coords1, atomap) result(dist2)
!   integer, intent(in) :: natom
!   integer, dimension(:), intent(in) :: atomap
!   real(wp), dimension(:), intent(in) :: weights
!   real(wp), dimension(:, :), intent(in) :: coords0, coords1
!   real(wp) :: eigmat(4, 4), eigval(4)
!
!   call kearsleymat(natom, weights, coords0, coords1, atomap, eigmat)
!   call syevec4(eigmat, eigval)
!   dist2 = squaredist(natom, weights, coords0, rotated(natom, coords1, eigmat(:, 1)), atomap)
!
!end function

! Purpose: Calculate rotation quaternion which minimzes the square distance
function leastrotquat(natom, weights, coords0, coords1, atomap) result(quat)
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: atomap
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   real(wp) :: quat(4), eigmat(4, 4), eigval(4)

   call kearsleymat(natom, weights, coords0, coords1, atomap, eigmat)
   call syevec4(eigmat, eigval)
   quat = eigmat(:, 1)

end function

! Purpose: Find the best orientation by least squares minimization
! Reference: Acta Cryst. (1989). A45, 208-210
subroutine kearsleymat(natom, weights, coords0, coords1, atomap, eigmat)
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: atomap
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords0, coords1

   integer :: i
   real(wp) :: eigmat(4, 4), p(3, natom), q(3, natom), auxmat(4, 4)

   eigmat = 0.0

   do i = 1, natom
      p(:, i) = coords0(:, i) + coords1(:, atomap(i))
      q(:, i) = coords0(:, i) - coords1(:, atomap(i))
   end do

! Calculate upper diagonal elements of the matrix

   do i = 1, natom
      auxmat(1, 1) = q(1, i)**2 + q(2, i)**2 + q(3, i)**2
      auxmat(1, 2) = p(2, i)*q(3, i) - q(2, i)*p(3, i)
      auxmat(1, 3) = q(1, i)*p(3, i) - p(1, i)*q(3, i)
      auxmat(1, 4) = p(1, i)*q(2, i) - q(1, i)*p(2, i)
      auxmat(2, 2) = p(2, i)**2 + p(3, i)**2 + q(1, i)**2
      auxmat(2, 3) = q(1, i)*q(2, i) - p(1, i)*p(2, i)
      auxmat(2, 4) = q(1, i)*q(3, i) - p(1, i)*p(3, i)
      auxmat(3, 3) = p(1, i)**2 + p(3, i)**2 + q(2, i)**2
      auxmat(3, 4) = q(2, i)*q(3, i) - p(2, i)*p(3, i)
      auxmat(4, 4) = p(1, i)**2 + p(2, i)**2 + q(3, i)**2
      eigmat = eigmat + weights(i)*auxmat
   end do

end subroutine

end module
