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

module spatial_transforms
use parameters
use random
use eigen

implicit none

private
public angle
public quatmul
public quatrotmat
public randrotquat
public weight_coords
public unweight_coords
public rotate_coords
public mirror_coords
public translate_coords
public centroid
public total_sqdist
public least_total_sqdist
public optimize_rotation

interface total_sqdist
   module procedure totsqdist_mol
   module procedure totsqdist_perm_mol
   module procedure totsqdist_perm_atoms
end interface

interface least_total_sqdist
   module procedure leastotsqdist_perm_atoms
end interface

interface optimize_rotation
   module procedure optrotquat_perm_mol
end interface

contains

function angle(q)
   real(rk), dimension(4), intent(in) :: q
   real(rk) :: angle

   angle = (180./asin(1.))*abs(atan(sqrt(sum(q(2:4)**2))/q(1)))
end function

function quatmul(p, q) result(pq)
   real(rk), dimension(4), intent(in) :: p, q
   real(rk) :: pq(4)

   ! Quaternion multiplication
   pq(1) = p(1)*q(1) - p(2)*q(2) - p(3)*q(3) - p(4)*q(4)
   pq(2) = p(1)*q(2) + p(2)*q(1) + p(3)*q(4) - p(4)*q(3)
   pq(3) = p(1)*q(3) - p(2)*q(4) + p(3)*q(1) + p(4)*q(2)
   pq(4) = p(1)*q(4) + p(2)*q(3) - p(3)*q(2) + p(4)*q(1)
end function

function quatrotmat(q) result(rotmat)
! Convert rotation quaternion to rotation matrix

   real(rk), intent(in) :: q(4)
   real(rk) :: rotmat(3, 3)

! Calculate the rotation matrix

   rotmat(1, 1) = 1.0_rk - 2*(q(3)**2 + q(4)**2)
   rotmat(2, 1) = 2*(q(2)*q(3) - q(1)*q(4))
   rotmat(3, 1) = 2*(q(2)*q(4) + q(1)*q(3))
   rotmat(1, 2) = 2*(q(2)*q(3) + q(1)*q(4))
   rotmat(2, 2) = 1.0_rk - 2*(q(2)**2 + q(4)**2)
   rotmat(3, 2) = 2*(q(3)*q(4) - q(1)*q(2))
   rotmat(1, 3) = 2*(q(2)*q(4) - q(1)*q(3))
   rotmat(2, 3) = 2*(q(3)*q(4) + q(1)*q(2))
   rotmat(3, 3) = 1.0_rk - 2*(q(2)**2 + q(3)**2)
end function

function randrotquat() result(rotquat)
! Description:
!    This function generates a random unit quaternion.
! References:
!    Academic Press Graphics Gems Series archive Graphics
!    Gems III archive. Pages: 129 - 132.
   real(rk) :: x(3)
   real(rk) :: rotquat(4)
   real(rk) :: pi, a1, a2, r1, r2, s1, s2, c1, c2

! Generate a random vector

   x = randvec()

! Calculate auxiliar vectors and constants

   pi = 2*asin(1.)
   a1 = 2*pi*x(1)
   a2 = 2*pi*x(2)
   r1 = sqrt(1.0 - x(3))
   r2 = sqrt(x(3))

   s1 = sin(a1)
   c1 = cos(a1)
   s2 = sin(a2)
   c2 = cos(a2)

   ! Unit quaternion (w, x, y, z)
   rotquat = [ c2*r2, s1*r1, c1*r1, s2*r2 ]
end function

subroutine weight_coords(coords, weights)
   real(rk), dimension(:,:), intent(inout) :: coords
   real(rk), dimension(:), intent(in) :: weights
   ! Local variables
   real(rk) :: total_weight
   integer :: i
   total_weight = sum(weights)
   do i = 1, size(coords, dim=2)
      coords(:, i) = coords(:, i)*sqrt(weights(i)/total_weight)
   end do
end subroutine

subroutine unweight_coords(coords, weights)
   real(rk), dimension(:,:), intent(inout) :: coords
   real(rk), dimension(:), intent(in) :: weights
   ! Local variables
   real(rk) :: total_weight
   integer :: i
   total_weight = sum(weights)
   do i = 1, size(coords, dim=2)
      coords(:, i) = coords(:, i)*sqrt(total_weight/weights(i))
   end do
end subroutine

subroutine translate_coords(coords, travec)
   real(rk), dimension(:,:), intent(inout) :: coords
   real(rk), dimension(3), intent(in) :: travec
   integer :: i

   do i = 1, size(coords, dim=2)
      coords(:, i) = coords(:, i) + travec(:)
   end do
end subroutine

subroutine rotate_coords(coords, rotquat, center)
   real(rk), dimension(:,:), intent(inout) :: coords
   real(rk), dimension(4), intent(in) :: rotquat
   real(rk), intent(in) :: center(3)
   ! Local variables
   real(rk) :: rotmat(3, 3), auxvec(3)
   integer :: i, j

   ! Covert quaternion to rotation matrix
   rotmat = quatrotmat(rotquat)

   ! Apply rotation
   do i = 1, size(coords, dim=2)
      auxvec(:) = 0.0
      do j = 1, 3
         auxvec(:) = auxvec(:) + rotmat(:, j)*(coords(j, i) - center(j))
      end do
      coords(:, i) = auxvec(:) + center(:)
   end do
end subroutine

subroutine mirror_coords(coords)
   real(rk), dimension(:,:), intent(inout) :: coords

   coords(1, :) = -coords(1, :)
end subroutine

function centroid(coords)
   real(rk), dimension(:,:), intent(in) :: coords
   ! Local variables
   integer :: i
   real(rk) :: centroid(3)

! Calculate the coordinates of the center of mass

   centroid(:) = 0

   do i = 1, size(coords, dim=2)
      centroid(:) = centroid(:) + coords(:, i)
   end do

   centroid(:) = centroid(:) / size(coords, dim=2)
end function

real(rk) function totsqdist_mol(coords1, coords2) result(total_sqdist)
   real(rk), dimension(:,:), intent(in) :: coords1, coords2

   total_sqdist = sum(sum((coords1 - coords2)**2, dim=1))
end function

real(rk) function totsqdist_perm_mol(atomperm, coords1, coords2) result(total_sqdist)
   integer, dimension(:), intent(in) :: atomperm
   real(rk), dimension(:,:), intent(in) :: coords1, coords2

   total_sqdist = sum(sum((coords1 - coords2(:, atomperm))**2, dim=1))
end function

real(rk) function totsqdist_perm_atoms(atomidcs, atomperm, coords1, coords2) result(total_sqdist)
   integer, dimension(:), intent(in) :: atomidcs, atomperm
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer :: i

   total_sqdist = 0

   do i = 1, size(atomidcs)
      total_sqdist = total_sqdist + sum((coords1(:, atomidcs(i)) - coords2(:, atomperm(atomidcs(i))))**2, dim=1)
   end do
end function

subroutine compute_residuals_matrix(coordsp, coordsm, residuals)
! Compute the 4x4 residuals matrix for optimal rotation calculation
! Reference: Acta Cryst. (1989). A45, 208-210
   real(rk), dimension(:,:), intent(in) :: coordsp, coordsm
   real(rk), dimension(4,4), intent(out) :: residuals
   integer :: i, num_atoms

   num_atoms = size(coordsp, dim=2)

   ! Initialize residuals matrix
   residuals = 0.0_rk

   ! Calculate upper matrix elements
   do i = 1, num_atoms
      residuals(1, 1) = residuals(1, 1) + (coordsm(1, i)**2 + coordsm(2, i)**2 + coordsm(3, i)**2)
      residuals(1, 2) = residuals(1, 2) + (coordsp(2, i)*coordsm(3, i) - coordsm(2, i)*coordsp(3, i))
      residuals(1, 3) = residuals(1, 3) + (coordsm(1, i)*coordsp(3, i) - coordsp(1, i)*coordsm(3, i))
      residuals(1, 4) = residuals(1, 4) + (coordsp(1, i)*coordsm(2, i) - coordsm(1, i)*coordsp(2, i))
      residuals(2, 2) = residuals(2, 2) + (coordsp(2, i)**2 + coordsp(3, i)**2 + coordsm(1, i)**2)
      residuals(2, 3) = residuals(2, 3) + (coordsm(1, i)*coordsm(2, i) - coordsp(1, i)*coordsp(2, i))
      residuals(2, 4) = residuals(2, 4) + (coordsm(1, i)*coordsm(3, i) - coordsp(1, i)*coordsp(3, i))
      residuals(3, 3) = residuals(3, 3) + (coordsp(1, i)**2 + coordsp(3, i)**2 + coordsm(2, i)**2)
      residuals(3, 4) = residuals(3, 4) + (coordsm(2, i)*coordsm(3, i) - coordsp(2, i)*coordsp(3, i))
      residuals(4, 4) = residuals(4, 4) + (coordsp(1, i)**2 + coordsp(2, i)**2 + coordsm(3, i)**2)
   end do

   ! Symmetrize matrix
   residuals(2, 1) = residuals(1, 2)
   residuals(3, 1) = residuals(1, 3)
   residuals(4, 1) = residuals(1, 4)
   residuals(3, 2) = residuals(2, 3)
   residuals(4, 2) = residuals(2, 4)
   residuals(4, 3) = residuals(3, 4)
end subroutine

subroutine optrotquat_perm_mol(atomperm, coords1, coords2, center, rotquat)
! Find the optimal rotation in quaternion representation by least squares minimization
! Reference: Acta Cryst. (1989). A45, 208-210
   integer, dimension(:), intent(in) :: atomperm
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   real(rk), intent(in) :: center(3)
   real(rk), intent(out) :: rotquat(4)
   ! Local variables
   integer :: i, num_atoms
   real(rk) :: residuals(4, 4)
   real(rk), dimension(:,:), allocatable :: coordsp, coordsm

   num_atoms = size(atomperm)

   allocate (coordsp(3, num_atoms))
   allocate (coordsm(3, num_atoms))

   do i = 1, num_atoms
      coordsp(:, i) = coords1(:, i) + coords2(:, atomperm(i)) - 2*center(:)
      coordsm(:, i) = coords1(:, i) - coords2(:, atomperm(i))
   end do

   ! Compute residuals matrix using the common procedure
   call compute_residuals_matrix(coordsp, coordsm, residuals)

   rotquat = leasteigvec(residuals)
end subroutine

function leastotsqdist_perm_atoms(atomidcs, atomperm, coords1, coords2) result(leastotsqdist)
! Find the optimal rotation in quaternion representation by least squares minimization
! Reference: Acta Cryst. (1989). A45, 208-210
   integer, dimension(:), intent(in) :: atomidcs, atomperm
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   ! Local variables
   real(rk) :: leastotsqdist
   integer :: i, num_atoms
   real(rk) :: residuals(4, 4)
   real(rk), dimension(:,:), allocatable :: coordsp, coordsm

   num_atoms = size(atomidcs)

   allocate (coordsp(3, num_atoms))
   allocate (coordsm(3, num_atoms))

   do i = 1, num_atoms
      coordsp(:, i) = coords1(:, atomidcs(i)) + coords2(:, atomperm(atomidcs(i)))
      coordsm(:, i) = coords1(:, atomidcs(i)) - coords2(:, atomperm(atomidcs(i)))
   end do

   ! Compute residuals matrix using the common procedure
   call compute_residuals_matrix(coordsp, coordsm, residuals)

   leastotsqdist = max(leasteigval(residuals), 0._rk)
end function

end module
