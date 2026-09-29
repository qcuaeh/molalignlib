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
! Rotations (as unit quaternions) and squared-distance measures between
! mapped coordinate sets
use parameters
use common_types
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
public ORIGIN
public IDENTITY_MATRIX
public MIRROR_MATRIX
public IDENTITY_QUATERNION

! Neutral values for coordinate transformations (see get_coords)
real(rk), parameter :: ORIGIN(3) = [0.0_rk, 0.0_rk, 0.0_rk]
real(rk), parameter :: IDENTITY_MATRIX(3,3) = reshape( &
      [1.0_rk, 0.0_rk, 0.0_rk, &
       0.0_rk, 1.0_rk, 0.0_rk, &
       0.0_rk, 0.0_rk, 1.0_rk], [3, 3])
! Reflection on the YZ plane (x -> -x)
real(rk), parameter :: MIRROR_MATRIX(3,3) = reshape( &
      [-1.0_rk, 0.0_rk, 0.0_rk, &
        0.0_rk, 1.0_rk, 0.0_rk, &
        0.0_rk, 0.0_rk, 1.0_rk], [3, 3])
! Unit quaternion (w, x, y, z) of the null rotation
real(rk), parameter :: IDENTITY_QUATERNION(4) = [1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]

! Atom mappings are complete permutations of the atoms in the coordinate
! arrays: atom i of coords1 maps to atom mapping1(i) of coords2. The
! subset variants are for partial mappings under construction, where only
! the atoms of molecule 1 listed in subset1 have been assigned.

interface sqdistsum
   module procedure sqdistsum_all
   module procedure sqdistsum_subset
end interface

interface sqdistmean
   module procedure sqdistmean_all
end interface

interface least_sqdistsum
   module procedure least_sqdistsum_all
   module procedure least_sqdistsum_subset
end interface

interface least_rotquat
   module procedure least_rotquat_all
end interface

interface rotate_coords
   module procedure rotate_coords_all
end interface

interface rotated_coords
   module procedure rotated_coords_center_all
end interface

contains

function angle(q)
! Rotation angle of the unit quaternion q, in degrees, in [0, 180]
   real(rk), dimension(4), intent(in) :: q
   real(rk) :: angle

   angle = (180./asin(1.))*atan2(sqrt(sum(q(2:4)**2)), abs(q(1)))
end function

function quatmul(p, q) result(pq)
! Quaternion product p*q
   real(rk), dimension(4), intent(in) :: p, q
   real(rk) :: pq(4)

   pq(1) = p(1)*q(1) - p(2)*q(2) - p(3)*q(3) - p(4)*q(4)
   pq(2) = p(1)*q(2) + p(2)*q(1) + p(3)*q(4) - p(4)*q(3)
   pq(3) = p(1)*q(3) - p(2)*q(4) + p(3)*q(1) + p(4)*q(2)
   pq(4) = p(1)*q(4) + p(2)*q(3) - p(3)*q(2) + p(4)*q(1)
end function

function quatrotmat(q) result(rotmat)
! Rotation matrix of a unit quaternion (w, x, y, z)
   real(rk), intent(in) :: q(4)
   real(rk) :: rotmat(3, 3)

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
! Uniformly distributed random rotation, as a unit quaternion (w, x, y, z)
! (Shoemake, Graphics Gems III, pp. 129-132)
   real(rk) :: x(3)
   real(rk) :: rotquat(4)
   real(rk) :: pi, a1, a2, r1, r2, s1, s2, c1, c2

   x = randvec()

   pi = 2*asin(1.)
   a1 = 2*pi*x(1)
   a2 = 2*pi*x(2)
   r1 = sqrt(1.0 - x(3))
   r2 = sqrt(x(3))

   s1 = sin(a1)
   c1 = cos(a1)
   s2 = sin(a2)
   c2 = cos(a2)

   rotquat = [ c2*r2, s1*r1, c1*r1, s2*r2 ]
end function

subroutine translate_coords(coords, travec)
! Add travec to every point
   real(rk), dimension(:,:), intent(inout) :: coords
   real(rk), dimension(3), intent(in) :: travec
   ! Local variables
   integer(ik) :: i

   do i = 1, size(coords, dim=2)
      coords(:, i) = coords(:, i) + travec(:)
   end do
end subroutine

subroutine rotate_coords_all(coords, rotquat)
! Rotate every point about the origin
   real(rk), dimension(:,:), intent(inout) :: coords
   real(rk), dimension(4), intent(in) :: rotquat
   ! Local variables
   real(rk) :: rotmat(3, 3), auxvec(3)
   integer(ik) :: i, j

   rotmat = quatrotmat(rotquat)

   do i = 1, size(coords, dim=2)
      auxvec(:) = 0
      do j = 1, 3
         auxvec(:) = auxvec(:) + rotmat(:,j)*(coords(j,i))
      end do
      coords(:,i) = auxvec(:)
   end do
end subroutine

function rotated_coords_center_all(coords, rotquat, center) result(rotated_coords)
! Copy of coords rotated about center
   real(rk), dimension(:,:), intent(in) :: coords
   real(rk), dimension(4), intent(in) :: rotquat
   real(rk), intent(in) :: center(3)
   ! Local variables
   real(rk), allocatable, dimension(:,:) :: rotated_coords
   real(rk) :: rotmat(3, 3)
   integer(ik) :: i, j

   allocate (rotated_coords, mold=coords)

   rotmat = quatrotmat(rotquat)

   do i = 1, size(coords, dim=2)
      rotated_coords(:,i) = center(:)
      do j = 1, 3
         rotated_coords(:,i) = rotated_coords(:,i) + rotmat(:,j)*(coords(j,i) - center(j))
      end do
   end do
end function

subroutine compute_residuals_matrix(coordsp, coordsm, residuals)
! Kearsley's 4x4 residual matrix (Acta Cryst. 1989, A45, 208-210) from the
! sums coordsp and differences coordsm of the mapped points. Its smallest
! eigenvalue is the least sum of squared distances over rotations about the
! origin and the corresponding eigenvector is the optimal rotation.
   real(rk), dimension(:,:), intent(in) :: coordsp, coordsm
   real(rk), dimension(4,4), intent(out) :: residuals
   integer(ik) :: i, n_atoms

   n_atoms = size(coordsp, dim=2)

   residuals = 0.0_rk

   ! Upper triangle
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

   residuals(2, 1) = residuals(1, 2)
   residuals(3, 1) = residuals(1, 3)
   residuals(4, 1) = residuals(1, 4)
   residuals(3, 2) = residuals(2, 3)
   residuals(4, 2) = residuals(2, 4)
   residuals(4, 3) = residuals(3, 4)
end subroutine

real(rk) function sqdistsum_all(mapping1, coords1, coords2) result(sqdistsum)
! Sum of squared distances between atom i of coords1 and atom mapping1(i)
! of coords2
   integer(ik), dimension(:), intent(in) :: mapping1
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   ! Local variables
   integer(ik) :: i

   sqdistsum = 0
   do i = 1, size(mapping1)
      sqdistsum = sqdistsum + sum((coords1(:, i) - coords2(:, mapping1(i)))**2, dim=1)
   end do
end function

real(rk) function sqdistsum_subset(subset1, mapping1, coords1, coords2) result(sqdistsum)
! As sqdistsum_all, over the atoms of coords1 listed in subset1
   integer(ik), dimension(:), intent(in) :: subset1, mapping1
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   ! Local variables
   integer(ik) :: i

   sqdistsum = 0
   do i = 1, size(subset1)
      sqdistsum = sqdistsum + sum((coords1(:, subset1(i)) - coords2(:, mapping1(subset1(i))))**2, dim=1)
   end do
end function

function least_sqdistsum_all(mapping1, coords1, coords2) result(leastotsqdist)
! Least sum of squared distances over all rotations about the origin
! (coordinates are expected to be centered)
   integer(ik), dimension(:), intent(in) :: mapping1
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   ! Local variables
   real(rk) :: leastotsqdist
   real(rk) :: residuals(4, 4)
   real(rk), dimension(:,:), allocatable :: coordsp, coordsm
   integer(ik) :: i

   allocate (coordsp(3, size(mapping1)))
   allocate (coordsm(3, size(mapping1)))

   do i = 1, size(mapping1)
      coordsp(:, i) = coords1(:, i) + coords2(:, mapping1(i))
      coordsm(:, i) = coords1(:, i) - coords2(:, mapping1(i))
   end do

   call compute_residuals_matrix(coordsp, coordsm, residuals)

   leastotsqdist = max(leasteigval(residuals), 0.0_rk)
end function

function least_sqdistsum_subset(subset1, mapping1, coords1, coords2) result(leastotsqdist)
! As least_sqdistsum_all, over the atoms of coords1 listed in subset1
   integer(ik), dimension(:), intent(in) :: subset1, mapping1
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   ! Local variables
   real(rk) :: leastotsqdist
   real(rk) :: residuals(4, 4)
   real(rk), dimension(:,:), allocatable :: coordsp, coordsm
   integer(ik) :: i

   allocate (coordsp(3, size(subset1)))
   allocate (coordsm(3, size(subset1)))

   do i = 1, size(subset1)
      coordsp(:, i) = coords1(:, subset1(i)) + coords2(:, mapping1(subset1(i)))
      coordsm(:, i) = coords1(:, subset1(i)) - coords2(:, mapping1(subset1(i)))
   end do

   call compute_residuals_matrix(coordsp, coordsm, residuals)

   leastotsqdist = max(leasteigval(residuals), 0.0_rk)
end function

real(rk) function sqdistmean_all(mapping1, weights, coords1, coords2) result(sqdistmean)
! Weighted mean of the squared distances between mapped atoms
   integer(ik), dimension(:), intent(in) :: mapping1
   real(rk), dimension(:), intent(in) :: weights
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   ! Local variables
   real(rk) :: total_weight, sqdistsum
   integer(ik) :: i

   total_weight = 0
   sqdistsum = 0
   do i = 1, size(mapping1)
      total_weight = total_weight + weights(i)
      sqdistsum = sqdistsum + weights(i)*sum((coords1(:, i) - coords2(:, mapping1(i)))**2, dim=1)
   end do
   sqdistmean = sqdistsum / total_weight
end function

function least_rotquat_all(mapping1, coords1, coords2) result(rotquat)
! Rotation about the origin that minimizes the sum of squared distances
! between mapped atoms when applied to coords2 (Kearsley's method)
   integer(ik), dimension(:), intent(in) :: mapping1
   real(rk), dimension(:,:), intent(in) :: coords1
   real(rk), dimension(:,:), intent(in) :: coords2
   ! Local variables
   real(rk), dimension(4) :: rotquat
   real(rk), dimension(:,:), allocatable :: coordsp, coordsm
   real(rk) :: residuals(4, 4)
   integer(ik) :: i

   allocate (coordsp(3, size(mapping1)))
   allocate (coordsm(3, size(mapping1)))

   do i = 1, size(mapping1)
      coordsp(:, i) = coords1(:, i) + coords2(:, mapping1(i))
      coordsm(:, i) = coords1(:, i) - coords2(:, mapping1(i))
   end do

   call compute_residuals_matrix(coordsp, coordsm, residuals)
   rotquat = leasteigvec(residuals)
end function

end module
