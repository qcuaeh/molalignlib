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
use kinds
use linear
use rotation

implicit none

private
public squaredist
public leastsquaredist
public leastrotquat

contains

real(wp) function squaredist(natom, weights, coords0, coords1, mapping)
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: mapping
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords0, coords1

   squaredist = sum(weights(1:natom)*sum((coords0(:, 1:natom) - coords1(:, mapping(1:natom)))**2, dim=1))

end function

real(wp) function leastsquaredist(natom, weights, coords0, coords1, mapping) result(squaredist)
! Purpose: Calculate least square distance from eigenvalues
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: mapping
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   real(wp) :: eigmat(4, 4), eigval(4)

   call kearsleymat(natom, weights, coords0, coords1, mapping, eigmat)
   call syeval4(eigmat, eigval)
   squaredist = max(eigval(1), 0._wp)

end function

!real(wp) function leastsquaredist(natom, weights, coords0, coords1, mapping) result(dist)
!! Purpose: Calculate least square distance from aligned coordinates
!   integer, intent(in) :: natom
!   integer, dimension(:), intent(in) :: mapping
!   real(wp), dimension(:), intent(in) :: weights
!   real(wp), dimension(:, :), intent(in) :: coords0, coords1
!   real(wp) :: eigmat(4, 4), eigval(4)
!
!   call kearsleymat(natom, weights, coords0, coords1, mapping, eigmat)
!   call syevec4(eigmat, eigval)
!   dist = squaredist(natom, weights, coords0, rotated(natom, coords1, eigmat(:, 1)), mapping)
!
!end function

function leastrotquat(natom, weights, coords0, coords1, mapping) result(rotquat)
! Purpose: Calculate rotation quaternion which minimzes the square distance
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: mapping
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   real(wp) :: eigmat(4, 4), eigval(4), rotquat(4)

   call kearsleymat(natom, weights, coords0, coords1, mapping, eigmat)
   call syevec4(eigmat, eigval)
   rotquat = eigmat(:, 1)

end function

subroutine kearsleymat(natom, weights, coords0, coords1, mapping, eigmat)
! Purpose: Find the best orientation by least squares minimization
! Reference: Acta Cryst. (1989). A45, 208-210
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: mapping
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords0, coords1

   integer :: i
   real(wp) :: eigmat(4, 4), p(3, natom), q(3, natom)

   eigmat = 0

   do i = 1, natom
      p(:, i) = coords0(:, i) + coords1(:, mapping(i))
      q(:, i) = coords0(:, i) - coords1(:, mapping(i))
   end do

   ! Calculate upper matrix elements

   do i = 1, natom
      eigmat(1, 1) = eigmat(1, 1) + weights(i)*(q(1, i)**2 + q(2, i)**2 + q(3, i)**2)
      eigmat(1, 2) = eigmat(1, 2) + weights(i)*(p(2, i)*q(3, i) - q(2, i)*p(3, i))
      eigmat(1, 3) = eigmat(1, 3) + weights(i)*(q(1, i)*p(3, i) - p(1, i)*q(3, i))
      eigmat(1, 4) = eigmat(1, 4) + weights(i)*(p(1, i)*q(2, i) - q(1, i)*p(2, i))
      eigmat(2, 2) = eigmat(2, 2) + weights(i)*(p(2, i)**2 + p(3, i)**2 + q(1, i)**2)
      eigmat(2, 3) = eigmat(2, 3) + weights(i)*(q(1, i)*q(2, i) - p(1, i)*p(2, i))
      eigmat(2, 4) = eigmat(2, 4) + weights(i)*(q(1, i)*q(3, i) - p(1, i)*p(3, i))
      eigmat(3, 3) = eigmat(3, 3) + weights(i)*(p(1, i)**2 + p(3, i)**2 + q(2, i)**2)
      eigmat(3, 4) = eigmat(3, 4) + weights(i)*(q(2, i)*q(3, i) - p(2, i)*p(3, i))
      eigmat(4, 4) = eigmat(4, 4) + weights(i)*(p(1, i)**2 + p(2, i)**2 + q(3, i)**2)
   end do

   ! Symmetrize matrix

   eigmat(2, 1) = eigmat(1, 2)
   eigmat(3, 1) = eigmat(1, 3)
   eigmat(4, 1) = eigmat(1, 4)
   eigmat(3, 2) = eigmat(2, 3)
   eigmat(4, 2) = eigmat(2, 4)
   eigmat(4, 3) = eigmat(3, 4)

end subroutine

end module
