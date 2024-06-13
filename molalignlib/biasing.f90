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
use kinds
use flags
use bounds
use sorting

implicit none

abstract interface
   subroutine proc_setcrossbias(natom, ntype, typelenlist, coords0, coords1, equivmat, biasmat)
      use kinds
      integer, intent(in) :: natom, ntype
      integer, dimension(:), intent(in) :: typelenlist
      real(wp), dimension(:, :), intent(in) :: coords0, coords1
      integer, dimension(:, :), intent(in) :: equivmat
      real(wp), dimension(:, :), intent(out) :: biasmat
   end subroutine
end interface

real(wp) :: bias_tol
real(wp) :: bias_scale
real(wp) :: bias_ratio

procedure(proc_setcrossbias), pointer :: setcrossbias
procedure(proc_setcrossbias), pointer :: mapsetcrossbias

contains

subroutine setcrossbias_none(natom, ntype, typelenlist, coords0, coords1, equivmat, biasmat)
! Purpose: Set biases from sorted neighbors' distances equivalence

   integer, intent(in) :: natom, ntype
   integer, dimension(:), intent(in) :: typelenlist
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   integer, dimension(:, :), intent(in) :: equivmat
   real(wp), dimension(:, :), intent(out) :: biasmat

   integer :: h, i, j, offset

   offset = 0

   do h = 1, ntype
      do i = offset + 1, offset + typelenlist(h)
         do j = offset + 1, offset + typelenlist(h)
            biasmat(j, i) = 0
         end do
      end do
      offset = offset + typelenlist(h)
   end do

end subroutine

subroutine setcrossbias_rd(natom, ntype, typelenlist, coords0, coords1, equivmat, biasmat)
! Purpose: Set biases from sorted neighbors' distances equivalence

   integer, intent(in) :: natom, ntype
   integer, dimension(:), intent(in) :: typelenlist
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   integer, dimension(:, :), intent(in) :: equivmat
   real(wp), dimension(:, :), intent(out) :: biasmat

   integer :: h, i, j, offset
   real(wp), allocatable :: d0(:, :), d1(:, :)

   allocate(d0(natom, natom), d1(natom, natom))

   do i = 1, natom
      offset = 0
      do h = 1, ntype
         do j = offset + 1, offset + typelenlist(h)
            d0(j, i) = sqrt(sum((coords0(:, j) - coords0(:, i))**2))
         end do
         call sort(d0(:, i), offset + 1, offset + typelenlist(h))
         offset = offset + typelenlist(h)
      end do
   end do

   do i = 1, natom
      offset = 0
      do h = 1, ntype
         do j = offset + 1, offset + typelenlist(h)
            d1(j, i) = sqrt(sum((coords1(:, j) - coords1(:, i))**2))
         end do
         call sort(d1(:, i), offset + 1, offset + typelenlist(h))
         offset = offset + typelenlist(h)
      end do
   end do

   offset = 0

   do h = 1, ntype
      do i = offset + 1, offset + typelenlist(h)
         do j = offset + 1, offset + typelenlist(h)
            if (all(abs(d1(:, j) - d0(:, i)) < bias_tol)) then
               biasmat(j, i) = 0
            else
               biasmat(j, i) = bias_scale**2
            end if
         end do
      end do
      offset = offset + typelenlist(h)
   end do

end subroutine

subroutine setcrossbias_mna(natom, ntype, typelenlist, coords0, coords1, equivmat, biasmat)
! Purpose: Set biases from sorted distances to neighbors equivalence

   integer, intent(in) :: natom, ntype
   integer, dimension(:), intent(in) :: typelenlist
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   integer, dimension(:, :), intent(in) :: equivmat
   real(wp), dimension(:, :), intent(out) :: biasmat
   integer :: h, i, j, offset, maxequiv

   maxequiv = 0

   offset = 0
   do h = 1, ntype
      do i = offset + 1, offset + typelenlist(h)
         do j = offset + 1, offset + typelenlist(h)
            maxequiv = max(maxequiv, equivmat(j, i))
         end do
      end do
      offset = offset + typelenlist(h)
   end do

   offset = 0
   do h = 1, ntype
      do i = offset + 1, offset + typelenlist(h)
         do j = offset + 1, offset + typelenlist(h)
!            biasmat(j, i) = bias_scale**2*bias_ratio**(equivmat(j, i))
            biasmat(j, i) = bias_scale**2*(maxequiv - equivmat(j, i))
         end do
      end do
      offset = offset + typelenlist(h)
   end do

end subroutine

end module
