! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

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

module euclidean
use parameters
use types_basic
use permutation
use random
use eigen
implicit none
private

public angle
public quatmul
public quatrotmat
public randrotquat
public sqdistmean
public sqdistsum
public rotate_coords
public rotated_coords
public translate_coords
public least_rotquat
public least_sqdistsum

interface sqdistsum
   module procedure sqdistsum_base
   module procedure sqdistsum_perm
   module procedure sqdistsum_subperm
end interface

interface sqdistmean
   module procedure sqdistmean_subperm
end interface

interface least_sqdistsum
   module procedure least_sqdistsum_subperm
end interface

interface rotate_coords
   module procedure rotate_coords_all
   module procedure rotate_coords_subset
end interface

interface rotated_coords
   module procedure rotated_coords_all
   module procedure rotated_coords_all_center
end interface

interface least_rotquat
   module procedure least_rotquat_base
   module procedure least_rotquat_perm
   module procedure least_rotquat_subperm
end interface

contains

function angle(q)
   real(rk), dimension(4), intent(in) :: q
   real(rk) :: angle

!   angle = (180./asin(1.))*atan2(sqrt(sum(q(2:4)**2)), q(1))
   angle = (180./asin(1.))*atan2(sqrt(sum(q(2:4)**2)), abs(q(1)))
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

subroutine translate_coords(coords, travec)
   real(rk), dimension(:,:), intent(inout) :: coords
   real(rk), dimension(3), intent(in) :: travec
   ! Local variables
   integer :: i

   do i = 1, size(coords, dim=2)
      coords(:, i) = coords(:, i) + travec(:)
   end do
end subroutine

subroutine rotate_coords_all(coords, rotquat)
   real(rk), dimension(:,:), intent(inout) :: coords
   real(rk), dimension(4), intent(in) :: rotquat
   ! Local variables
   real(rk) :: rotmat(3, 3), auxvec(3)
   integer :: i, j

   ! Convert quaternion to rotation matrix
   rotmat = quatrotmat(rotquat)

   ! Apply rotation
   do i = 1, size(coords, dim=2)
      auxvec(:) = 0
      do j = 1, 3
         auxvec(:) = auxvec(:) + rotmat(:,j)*(coords(j,i))
      end do
      coords(:,i) = auxvec(:)
   end do
end subroutine

subroutine rotate_coords_subset(atomset1, coords, rotquat)
   integer, dimension(:), intent(in) :: atomset1
   real(rk), dimension(:,:), intent(inout) :: coords
   real(rk), dimension(4), intent(in) :: rotquat
   ! Local variables
   real(rk) :: rotmat(3, 3), auxvec(3)
   integer :: i, j

   ! Convert quaternion to rotation matrix
   rotmat = quatrotmat(rotquat)

   ! Apply rotation
   do i = 1, size(atomset1)
      auxvec(:) = 0
      do j = 1, 3
         auxvec(:) = auxvec(:) + rotmat(:,j)*(coords(j,atomset1(i)))
      end do
      coords(:,atomset1(i)) = auxvec(:)
   end do
end subroutine

function rotated_coords_all(coords, rotquat) result(rotated_coords)
   real(rk), dimension(:,:), intent(in) :: coords
   real(rk), dimension(4), intent(in) :: rotquat
   ! Local variables
   real(rk), allocatable, dimension(:,:) :: rotated_coords
   real(rk) :: rotmat(3, 3)
   integer :: i, j

   allocate (rotated_coords, mold=coords)

   ! Convert quaternion to rotation matrix
   rotmat = quatrotmat(rotquat)

   ! Apply rotation
   do i = 1, size(coords, dim=2)
      rotated_coords(:,i) = 0
      do j = 1, 3
         rotated_coords(:,i) = rotated_coords(:,i) + rotmat(:,j)*(coords(j,i))
      end do
   end do
end function

function rotated_coords_all_center(coords, rotquat, center) result(rotated_coords)
   real(rk), dimension(:,:), intent(in) :: coords
   real(rk), dimension(4), intent(in) :: rotquat
   real(rk), intent(in) :: center(3)
   ! Local variables
   real(rk), allocatable, dimension(:,:) :: rotated_coords
   real(rk) :: rotmat(3, 3)
   integer :: i, j

   allocate (rotated_coords, mold=coords)

   ! Convert quaternion to rotation matrix
   rotmat = quatrotmat(rotquat)

   ! Apply rotation
   do i = 1, size(coords, dim=2)
      rotated_coords(:,i) = center(:)
      do j = 1, 3
         rotated_coords(:,i) = rotated_coords(:,i) + rotmat(:,j)*(coords(j,i) - center(j))
      end do
   end do
end function

real(rk) function sqdistsum_base(coords1, coords2) result(sqdistsum)
   real(rk), dimension(:,:), intent(in) :: coords1, coords2

   sqdistsum = sum(sum((coords1 - coords2)**2, dim=1))
end function

real(rk) function sqdistsum_perm(atomperm1, coords1, coords2) result(sqdistsum)
   integer, dimension(:), intent(in) :: atomperm1
   real(rk), dimension(:,:), intent(in) :: coords1, coords2

   sqdistsum = sum(sum((coords1 - coords2(:, atomperm1))**2, dim=1))
end function

subroutine compute_residuals_matrix(coordsp, coordsm, residuals)
! Compute the 4x4 residuals matrix for optimal rotation calculation
! Reference: Acta Cryst. (1989). A45, 208-210
   real(rk), dimension(:,:), intent(in) :: coordsp, coordsm
   real(rk), dimension(4,4), intent(out) :: residuals
   integer :: i, n_atoms

   n_atoms = size(coordsp, dim=2)

   ! Initialize residuals matrix
   residuals = 0.0_rk

   ! Calculate upper matrix elements
   do i = 1, n_atoms
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

function least_rotquat_base(coords1, coords2) result(rotquat)
! Find the optimal rotation in quaternion representation by least squares minimization
! Reference: Acta Cryst. (1989). A45, 208-210
   real(rk), dimension(:,:), intent(in) :: coords1
   real(rk), dimension(:,:), intent(in) :: coords2
   ! Local variables
   real(rk), dimension(4) :: rotquat
   real(rk), dimension(:,:), allocatable :: coordsp, coordsm
   real(rk) :: residuals(4, 4)
   integer :: i, n_atoms

   n_atoms = size(coords1, dim=2)

   allocate (coordsp(3, n_atoms))
   allocate (coordsm(3, n_atoms))

   do i = 1, n_atoms
      coordsp(:, i) = coords1(:, i) + coords2(:, i)
      coordsm(:, i) = coords1(:, i) - coords2(:, i)
   end do

   ! Compute residuals matrix using the common procedure
   call compute_residuals_matrix(coordsp, coordsm, residuals)
   rotquat = leasteigvec(residuals)
end function

function least_rotquat_perm(atomperm1, coords1, coords2) result(rotquat)
! Find the optimal rotation in quaternion representation by least squares minimization
! Reference: Acta Cryst. (1989). A45, 208-210
   integer, dimension(:), intent(in) :: atomperm1
   real(rk), dimension(:,:), intent(in) :: coords1
   real(rk), dimension(:,:), intent(in) :: coords2
   ! Local variables
   real(rk), dimension(4) :: rotquat
   real(rk), dimension(:,:), allocatable :: coordsp, coordsm
   real(rk) :: residuals(4, 4)
   integer :: i, n_atoms

   n_atoms = size(atomperm1)

   allocate (coordsp(3, n_atoms))
   allocate (coordsm(3, n_atoms))

   do i = 1, n_atoms
      coordsp(:, i) = coords1(:, i) + coords2(:, atomperm1(i))
      coordsm(:, i) = coords1(:, i) - coords2(:, atomperm1(i))
   end do

   ! Compute residuals matrix using the common procedure
   call compute_residuals_matrix(coordsp, coordsm, residuals)
   rotquat = leasteigvec(residuals)
end function

real(rk) function sqdistsum_subperm(atomset1, atomperm1, coords1, coords2) result(sqdistsum)
   integer, dimension(:), intent(in) :: atomset1, atomperm1
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   ! Local variables
   integer :: i

   sqdistsum = 0
   do i = 1, size(atomset1)
      sqdistsum = sqdistsum + sum((coords1(:, atomset1(i)) - coords2(:, atomperm1(atomset1(i))))**2, dim=1)
   end do
end function

function least_sqdistsum_subperm(atomset1, atomperm1, coords1, coords2) result(leastotsqdist)
   integer, dimension(:), intent(in) :: atomset1, atomperm1
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   ! Local variables
   real(rk) :: leastotsqdist
   real(rk) :: residuals(4, 4)
   real(rk), dimension(:,:), allocatable :: coordsp, coordsm
   integer :: i

   allocate (coordsp(3, size(atomset1)))
   allocate (coordsm(3, size(atomset1)))

   do i = 1, size(atomset1)
      coordsp(:, i) = coords1(:, atomset1(i)) + coords2(:, atomperm1(atomset1(i)))
      coordsm(:, i) = coords1(:, atomset1(i)) - coords2(:, atomperm1(atomset1(i)))
   end do

   ! Compute residuals matrix using the common procedure
   call compute_residuals_matrix(coordsp, coordsm, residuals)

   leastotsqdist = max(leasteigval(residuals), 0._rk)
end function

real(rk) function sqdistmean_subperm(atomset1, atomperm1, weights, coords1, coords2) result(sqdistmean)
   integer, dimension(:), intent(in) :: atomset1, atomperm1
   real(rk), dimension(:), intent(in) :: weights
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   ! Local variables
   real(rk) :: total_weight, sqdistsum
   integer :: i

   total_weight = 0
   sqdistsum = 0
   do i = 1, size(atomset1)
      total_weight = total_weight + weights(atomset1(i))
      sqdistsum = sqdistsum + weights(atomset1(i))*sum((coords1(:, atomset1(i)) - coords2(:, atomperm1(atomset1(i))))**2, dim=1)
   end do
   sqdistmean = sqdistsum / total_weight
end function

function least_rotquat_subperm(atomset1, atomperm1, coords1, coords2) result(rotquat)
! Find the optimal rotation in quaternion representation by least squares minimization
! Reference: Acta Cryst. (1989). A45, 208-210
   integer, dimension(:), intent(in) :: atomset1, atomperm1
   real(rk), dimension(:,:), intent(in) :: coords1
   real(rk), dimension(:,:), intent(in) :: coords2
   ! Local variables
   real(rk), dimension(4) :: rotquat
   real(rk), dimension(:,:), allocatable :: coordsp, coordsm
   real(rk) :: residuals(4, 4)
   integer :: i

   allocate (coordsp(3, size(atomset1)))
   allocate (coordsm(3, size(atomset1)))

   do i = 1, size(atomset1)
      coordsp(:, i) = coords1(:, atomset1(i)) + coords2(:, atomperm1(atomset1(i)))
      coordsm(:, i) = coords1(:, atomset1(i)) - coords2(:, atomperm1(atomset1(i)))
   end do

   ! Compute residuals matrix using the common procedure
   call compute_residuals_matrix(coordsp, coordsm, residuals)
!   sqdistsum = max(leasteigval(residuals), 0._rk)
   rotquat = leasteigvec(residuals)
end function

end module
