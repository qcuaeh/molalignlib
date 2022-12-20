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

module translation
use settings

implicit none

private
public centroid
public centered
public translate
public translated

contains

function centroid(natom, weights, coords)
! Purpose: Get coordinates of centroid
   integer, intent(in) :: natom
   real(wp), intent(in) :: weights(natom), coords(3, natom)
   real(wp) :: centroid(3)
   integer :: i

! Calculate the coordinates of the center of mass

   centroid(:) = 0

   do i = 1, natom
      centroid(:) = centroid(:) + weights(i)*coords(:, i)
   end do

   centroid(:) = centroid(:)/sum(weights)

end function

function translated(natom, coords, vector)
! Purpose: Translate atomic coordinates to its vector of geometry
   integer, intent(in) :: natom
   real(wp), intent(in) :: coords(3, natom), vector(3)
   real(wp) :: translated(3, natom)
   integer :: i

   do i = 1, natom
      translated(:, i) = coords(:, i) + vector(:)
   end do

end function

subroutine translate(natom, coords, vector)
! Purpose: Translate atomic coordinates to its vector of geometry
   integer, intent(in) :: natom
   real(wp), intent(in) :: vector(3)
   real(wp), intent(inout) :: coords(3, natom)
   integer :: i

   do i = 1, natom
      coords(:, i) = coords(:, i) + vector(:)
   end do

end subroutine

function centered(natom, coords, vector)
! Purpose: Translate atomic coordinates to its vector of geometry
   integer, intent(in) :: natom
   real(wp), intent(in) :: coords(3, natom), vector(3)
   real(wp) :: centered(3, natom)
   integer :: i

   do i = 1, natom
      centered(:, i) = coords(:, i) - vector(:)
   end do

end function

end module
