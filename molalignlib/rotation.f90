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

module rotation
use kinds

implicit none

private
public rotate
public rotated
public matrix_rotated
public quater_rotated
public quat2rotmat
public quatmul
public rotangle
public genrotmat
public genrotquat

interface rotate
   module procedure matrix_rotate
   module procedure quater_rotate
end interface

interface rotated
   module procedure matrix_rotated
   module procedure quater_rotated
end interface

contains

function quat2rotmat(rotquat) result(rotmat)
! Purpose: Apply rotation quaternion to atomic coordinates
   real(rk), intent(in) :: rotquat(4)
   real(rk) :: rotmat(3, 3)

! Calculate the rotation matrix

   rotmat(1, 1) = 1.0_rk - 2*(rotquat(3)**2 + rotquat(4)**2)
   rotmat(1, 2) = 2*(rotquat(2)*rotquat(3) + rotquat(1)*rotquat(4))
   rotmat(1, 3) = 2*(rotquat(2)*rotquat(4) - rotquat(1)*rotquat(3))
   rotmat(2, 1) = 2*(rotquat(2)*rotquat(3) - rotquat(1)*rotquat(4))
   rotmat(2, 2) = 1.0_rk - 2*(rotquat(2)**2 + rotquat(4)**2)
   rotmat(2, 3) = 2*(rotquat(3)*rotquat(4) + rotquat(1)*rotquat(2))
   rotmat(3, 1) = 2*(rotquat(2)*rotquat(4) + rotquat(1)*rotquat(3))
   rotmat(3, 2) = 2*(rotquat(3)*rotquat(4) - rotquat(1)*rotquat(2))
   rotmat(3, 3) = 1.0_rk - 2*(rotquat(2)**2 + rotquat(3)**2)

end function

subroutine matrix_rotate(natom, coords, rotmat)
! Purpose: Apply rotation matrix to atomic coordinates
   integer, intent(in) :: natom
   real(rk), intent (in) :: rotmat(3, 3)
   real(rk), intent (inout) :: coords(3, natom)

   coords = matmul(rotmat, coords)

end subroutine

function matrix_rotated(natom, coords, rotmat) result(rotcoords)
! Purpose: Apply rotation matrix to atomic coordinates
   integer, intent(in) :: natom
   real(rk), intent (in) :: rotmat(3, 3)
   real(rk), intent (in) :: coords(3, natom)
   real(rk) :: rotcoords(3, natom)

   rotcoords = matmul(rotmat, coords)

end function

subroutine quater_rotate(natom, coords, rotquat)
! Purpose: Apply rotation quaternion to atomic coordinates
   integer, intent(in) :: natom
   real(rk), intent(in) :: rotquat(4)
   real(rk), intent(inout) :: coords(3, natom)

! Rotate atomic coordinates

!    call matrix_rotate(natom, quat2rotmat(rotquat), coords)
    coords = matmul(quat2rotmat(rotquat), coords)

end subroutine

function quater_rotated(natom, coords, rotquat) result(rotcoords)
! Purpose: Apply rotation quaternion to atomic coordinates
   integer, intent(in) :: natom
   real(rk), intent(in) :: rotquat(4)
   real(rk), intent(in) :: coords(3, natom)
   real(rk) :: rotcoords(3, natom)

   rotcoords = matmul(quat2rotmat(rotquat), coords)

end function

function quatmul(p, q) result(prod)

   real(rk), dimension(4), intent(in) :: p, q
   real(rk) :: prod(4)
   real(rk) :: left(4, 4)

   left(:, 1) = [ p(1), -p(2), -p(3), -p(4) ]
   left(:, 2) = [ p(2),  p(1), -p(4),  p(3) ]
   left(:, 3) = [ p(3),  p(4),  p(1), -p(2) ]
   left(:, 4) = [ p(4), -p(3),  p(2),  p(1) ]

   prod = matmul(left, q)

end function

function rotangle(q) result(angle)
   real(rk), dimension(4), intent(in) :: q
   real(rk) :: angle

!    angle = 2*acos(q(1))
!    angle = 2*asin(sqrt(sum(q(2:4)**2)))
   angle = 2*abs(atan(sqrt(sum(q(2:4)**2))/q(1)))
!    angle = 2*atan2(sqrt(sum(q(2:4)**2)), q(1))

!    print *, sum(q**2), &
!        90./asin(1.)*2*acos(q(1)), &
!        90./asin(1.)*2*asin(sqrt(sum(q(2:4)**2))), &
!        90./asin(1.)*2*atan(sqrt(sum(q(2:4)**2))/q(1)), &
!        90./asin(1.)*2*atan2(sqrt(sum(q(2:4)**2)), q(1))

end function

function genrotquat(x) result(rotquat)
! Purpose: Generate a random rotation quaternion.
! Reference: Academic Press Graphics Gems Series archive Graphics
!            Gems III archive. Pages: 129 - 132.
   real(rk), intent(in) :: x(3)
   real(rk) :: rotquat(4)
   real(rk) :: pi, a1, a2, r1, r2, s1, s2, c1, c2

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

!   Quaternion (w, x, y, z), w must be first as required by quat2rotmat
   rotquat = [ c2*r2, s1*r1, c1*r1, s2*r2 ]

!    print *, sum(rotquat**2)

end function

function genrotmat(x) result(rotmat)
! Purpose: Generate a random rotation matrix.
! Reference: Academic Press Graphics Gems Series archive Graphics
!            Gems III archive. Pages: 117 - 120.
   real(rk), intent(in) :: x(3)
   real(rk) :: rotmat(3, 3)

   integer :: i
   real(rk) :: pi, a1, a2, r1, r2, s1, c1, s2, c2
   real(rk) :: v(3), right(3, 3), left(3, 3)

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

   right(1, :) = [ c1, s1, 0.0_rk ]
   right(2, :) = [ -s1, c1, 0.0_rk ]
   right(3, :) = [ 0.0_rk, 0.0_rk, 1.0_rk ]

   v = [ c2*r2, s2*r2, r1 ]

   do i = 1, 3
      left(:, i) = 2*v(i)*v(:)
   end do

   do i = 1, 3
      left(i, i) = left(i, i) - 1.0_rk
   end do

! Calculate the rotation matrix

   rotmat = matmul(left, right)

!    print *, matdet3(rotmat)

end function

end module
