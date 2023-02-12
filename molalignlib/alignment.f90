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
public biasdist
public leastsquaredist
public leastrotquat

contains

real(wp) function biasdist(natom, weights, biasmat, atomperm) result(dist)
   integer, intent(in) :: natom
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: biasmat
   integer, dimension(:), intent(in) :: atomperm
   integer :: i

   dist = 0

   do i = 1, natom
      dist = dist + weights(i)*biasmat(i, atomperm(i))
   end do

end function

real(wp) function squaredist(natom, weights, coords0, coords1, atomperm) result(dist)
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: atomperm
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords0, coords1

   dist = sum(weights(1:natom)*sum((coords0(:, 1:natom) - coords1(:, atomperm(1:natom)))**2, dim=1))

end function

real(wp) function leastsquaredist(natom, weights, coords0, coords1, atomperm) result(dist)
! Purpose: Calculate least square distance from eigenvalues
    integer, intent(in) :: natom
    integer, dimension(:), intent(in) :: atomperm
    real(wp), dimension(:), intent(in) :: weights
    real(wp), dimension(:, :), intent(in) :: coords0, coords1
    real(wp) eigmat(4, 4), eigval(4)

    call kearsleymat(natom, weights, coords0, coords1, atomperm, eigmat)
    call syeval4(eigmat, eigval)
    dist = max(eigval(1), 0._wp)

end function

!real(wp) function leastsquaredist(natom, weights, coords0, coords1, atomperm) result(dist)
!! Purpose: Calculate least square distance from aligned coordinates
!   integer, intent(in) :: natom
!   integer, dimension(:), intent(in) :: atomperm
!   real(wp), dimension(:), intent(in) :: weights
!   real(wp), dimension(:, :), intent(in) :: coords0, coords1
!   real(wp) :: eigmat(4, 4), eigval(4)
!
!   call kearsleymat(natom, weights, coords0, coords1, atomperm, eigmat)
!   call syevec4(eigmat, eigval)
!   dist = squaredist(natom, weights, coords0, rotated(natom, coords1, eigmat(:, 1)), atomperm)
!
!end function

function leastrotquat(natom, weights, coords0, coords1, atomperm) result(quat)
! Purpose: Calculate rotation quaternion which minimzes the square distance
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: atomperm
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   real(wp) :: quat(4), eigmat(4, 4), eigval(4)

   call kearsleymat(natom, weights, coords0, coords1, atomperm, eigmat)
   call syevec4(eigmat, eigval)
   quat = eigmat(:, 1)

end function

subroutine kearsleymat(natom, weights, coords0, coords1, atomperm, eigmat)
! Purpose: Find the best orientation by least squares minimization
! Reference: Acta Cryst. (1989). A45, 208-210
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: atomperm
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords0, coords1

   integer :: i
   real(wp) :: eigmat(4, 4), p(3, natom), q(3, natom), auxmat(4, 4)

   eigmat = 0

   do i = 1, natom
      p(:, i) = coords0(:, i) + coords1(:, atomperm(i))
      q(:, i) = coords0(:, i) - coords1(:, atomperm(i))
   end do

! Calculate uppercase diagonal elements of the matrix

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
