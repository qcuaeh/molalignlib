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

module rigid_body
use parameters
use rotation

implicit none

private
public rotate_coords
public rotated_coords
public mirror_coords
public reflected_coords
public translate_coords
public translated_coords
public centroid

contains

subroutine translate_coords(coords, travec)
   real(rk), dimension(:, :), intent(inout) :: coords
   real(rk), dimension(3), intent(in) :: travec
   integer :: i

   do i = 1, size(coords, dim=2)
      coords(:, i) = coords(:, i) + travec(:)
   end do

end subroutine

function translated_coords(coords, travec)
   real(rk), dimension(:, :), intent(in) :: coords
   real(rk), dimension(3), intent(in) :: travec
   ! Local variables
   real(rk), allocatable :: translated_coords(:, :)
   integer :: i

   allocate (translated_coords(3, size(coords, dim=2)))

   do i = 1, size(coords, dim=2)
      translated_coords(:, i) = coords(:, i) + travec(:)
   end do

end function

subroutine rotate_coords(coords, q)
   real(rk), dimension(:, :), intent(inout) :: coords
   real(rk), dimension(4), intent(in) :: q
   ! Local variables
   real(rk) :: rotmat(3, 3)

! Rotate atomic coordinates

   rotmat = quatrotmat(q)
   coords = matmul(rotmat, coords)

end subroutine

function rotated_coords(coords, q)
   real(rk), dimension(:, :), intent(in) :: coords
   real(rk), dimension(4), intent(in) :: q
   ! Local variables
   real(rk) :: rotmat(3, 3)
   real(rk), allocatable :: rotated_coords(:, :)

   rotmat = quatrotmat(q)
   rotated_coords = matmul(rotmat, coords)

end function

subroutine mirror_coords(coords)
   real(rk), dimension(:, :), intent(inout) :: coords

   coords(1, :) = -coords(1, :)

end subroutine

function reflected_coords(coords)
   real(rk), dimension(:, :), intent(in) :: coords
   ! Local variables
   real(rk), allocatable :: reflected_coords(:, :)

   reflected_coords(1, :) = -coords(1, :)
   reflected_coords(2:3, :) = coords(2:3, :)

end function

function centroid(coords)
   real(rk), dimension(:, :), intent(in) :: coords
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

end module
