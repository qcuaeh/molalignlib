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
use random

implicit none

private
public rotate
public rotated
public quatmul
public rotangle
public randrotquat

contains

function quatrotmat(q) result(rotmat)
! Purpose: Apply rotation quaternion to atomic coordinates
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

subroutine rotate(natom, coords, q)
! Purpose: Apply rotation quaternion to atomic coordinates
   integer, intent(in) :: natom
   real(rk), intent(in) :: q(4)
   real(rk), intent(inout) :: coords(3, natom)
   real(rk) :: rotmat(3, 3)

! Rotate atomic coordinates

   rotmat = quatrotmat(q)
   coords = matmul(rotmat, coords)

end subroutine

function rotated(natom, coords, q) result(rotated_coords)
! Purpose: Apply rotation quaternion to atomic coordinates
   integer, intent(in) :: natom
   real(rk), intent(in) :: q(4)
   real(rk), intent(in) :: coords(3, natom)

   integer :: i
   real(rk) :: rotmat(3, 3)
   real(rk) :: rotated_coords(3, natom)

   rotmat = quatrotmat(q)
   rotated_coords = matmul(rotmat, coords)

end function

function quatmul(p, q) result(pq)

   real(rk), dimension(4), intent(in) :: p, q
   real(rk) :: pq(4)
!   real(rk) :: matp(4, 4)

! This matrix was wrongly defined in previous versions
! (it was erroneously transposed), now it is correct
!   matp(1, :) = [ p(1), -p(2), -p(3), -p(4) ]
!   matp(2, :) = [ p(2),  p(1), -p(4),  p(3) ]
!   matp(3, :) = [ p(3),  p(4),  p(1), -p(2) ]
!   matp(4, :) = [ p(4), -p(3),  p(2),  p(1) ]
!
!   pq = matmul(matp, q)

! Better to use direct quaternion multiplication
   pq(1) = p(1)*q(1) - p(2)*q(2) - p(3)*q(3) - p(4)*q(4)
   pq(2) = p(1)*q(2) + p(2)*q(1) + p(3)*q(4) - p(4)*q(3)
   pq(3) = p(1)*q(3) - p(2)*q(4) + p(3)*q(1) + p(4)*q(2)
   pq(4) = p(1)*q(4) + p(2)*q(3) - p(3)*q(2) + p(4)*q(1)

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

function randrotquat() result(rotquat)
! Purpose: Generate a random unit quaternion.
! Reference: Academic Press Graphics Gems Series archive Graphics
!            Gems III archive. Pages: 129 - 132.
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

!   Quaternion (w, x, y, z), w must be first as required by quatrotmat
   rotquat = [ c2*r2, s1*r1, c1*r1, s2*r2 ]

!    print *, sum(rotquat**2)

end function

end module
