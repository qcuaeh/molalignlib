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

module biasing

use settings
use sorting

implicit none

contains

subroutine setadjbias(natom, nblock, bsize, coords0, coords1, bias)
   integer, intent(in) :: natom, nblock
   integer, dimension(:), intent(in) :: bsize
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   real(wp), dimension(:, :), intent(out) :: bias

   integer :: h, i, j, offset
!    real(wp) d0(natom, natom), d1(natom, natom)
   real(wp), allocatable :: d0(:, :), d1(:, :)
   allocate(d0(natom, natom), d1(natom, natom))

   do i = 1, natom
      offset = 0
      do h = 1, nblock
         do j = offset + 1, offset + bsize(h)
            d0(j, i) = sqrt(sum((coords0(:, j) - coords0(:, i))**2))
         end do
         call sort(d0(:, i), offset + 1, offset + bsize(h))
         offset = offset + bsize(h)
      end do
   end do

   do i = 1, natom
      offset = 0
      do h = 1, nblock
         do j = offset + 1, offset + bsize(h)
            d1(j, i) = sqrt(sum((coords1(:, j) - coords1(:, i))**2))
         end do
         call sort(d1(:, i), offset + 1, offset + bsize(h))
         offset = offset + bsize(h)
      end do
   end do

   offset = 0
   bias(:, :) = 0.

   if (bias_flag) then
      do h = 1, nblock
         do i = offset + 1, offset + bsize(h)
            do j = offset + 1, offset + bsize(h)
               if (any(abs(d1(:, j) - d0(:, i)) > bias_tol)) then
                  bias(i, j) = bias_scale**2
               end if
            end do
         end do
         offset = offset + bsize(h)
      end do
   end if

end subroutine

end module
