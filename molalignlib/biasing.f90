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
   subroutine f_cross(natom, neltype, eltypepartlens, coords0, coords1, equivmat, prunemat, biasmat)
      use kinds
      integer, intent(in) :: natom, neltype
      integer, dimension(:), intent(in) :: eltypepartlens
      real(rk), dimension(:, :), intent(in) :: coords0, coords1
      integer, dimension(:, :), intent(in) :: equivmat
      logical, dimension(:, :), intent(out) :: prunemat
      real(rk), dimension(:, :), intent(out) :: biasmat
   end subroutine
end interface

real(rk) :: bias_tol
real(rk) :: bias_scale
real(rk) :: bias_ratio

procedure(f_cross), pointer :: cross_function

contains

subroutine cross_none(natom, neltype, eltypepartlens, coords0, coords1, equivmat, prunemat, biasmat)
! Purpose: Set biases from sorted neighbors' distances equivalence

   integer, intent(in) :: natom, neltype
   integer, dimension(:), intent(in) :: eltypepartlens
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   integer, dimension(:, :), intent(in) :: equivmat
   logical, dimension(:, :), intent(out) :: prunemat
   real(rk), dimension(:, :), intent(out) :: biasmat

   integer :: h, i, j, offset

   offset = 0

   do h = 1, neltype
      do i = offset + 1, offset + eltypepartlens(h)
         do j = offset + 1, offset + eltypepartlens(h)
            biasmat(j, i) = 0
            prunemat(j, i) = .true.
         end do
      end do
      offset = offset + eltypepartlens(h)
   end do

end subroutine

subroutine cross_prune_rd(natom, neltype, eltypepartlens, coords0, coords1, equivmat, prunemat, biasmat)
! Purpose: Set biases from sorted neighbors' distances equivalence

   integer, intent(in) :: natom, neltype
   integer, dimension(:), intent(in) :: eltypepartlens
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   integer, dimension(:, :), intent(in) :: equivmat
   logical, dimension(:, :), intent(out) :: prunemat
   real(rk), dimension(:, :), intent(out) :: biasmat

   integer :: h, i, j, offset
   real(rk), allocatable :: d0(:, :), d1(:, :)

   allocate (d0(natom, natom), d1(natom, natom))

   do i = 1, natom
      offset = 0
      do h = 1, neltype
         do j = offset + 1, offset + eltypepartlens(h)
            d0(j, i) = sqrt(sum((coords0(:, j) - coords0(:, i))**2))
         end do
         call sort(d0(:, i), offset + 1, offset + eltypepartlens(h))
         offset = offset + eltypepartlens(h)
      end do
   end do

   do i = 1, natom
      offset = 0
      do h = 1, neltype
         do j = offset + 1, offset + eltypepartlens(h)
            d1(j, i) = sqrt(sum((coords1(:, j) - coords1(:, i))**2))
         end do
         call sort(d1(:, i), offset + 1, offset + eltypepartlens(h))
         offset = offset + eltypepartlens(h)
      end do
   end do

   offset = 0

   do h = 1, neltype
      do i = offset + 1, offset + eltypepartlens(h)
         do j = offset + 1, offset + eltypepartlens(h)
            if (all(abs(d1(:, j) - d0(:, i)) < bias_tol)) then
               prunemat(j, i) = .true.
            else
               prunemat(j, i) = .false.
            end if
         end do
      end do
      offset = offset + eltypepartlens(h)
   end do

end subroutine

subroutine cross_bias_mna(natom, neltype, eltypepartlens, coords0, coords1, equivmat, prunemat, biasmat)
! Purpose: Set biases from sorted distances to neighbors equivalence

   integer, intent(in) :: natom, neltype
   integer, dimension(:), intent(in) :: eltypepartlens
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   integer, dimension(:, :), intent(in) :: equivmat
   logical, dimension(:, :), intent(out) :: prunemat
   real(rk), dimension(:, :), intent(out) :: biasmat
   integer :: h, i, j, offset, maxequiv

   maxequiv = 0

   offset = 0
   do h = 1, neltype
      do i = offset + 1, offset + eltypepartlens(h)
         do j = offset + 1, offset + eltypepartlens(h)
            maxequiv = max(maxequiv, equivmat(j, i))
         end do
      end do
      offset = offset + eltypepartlens(h)
   end do

   offset = 0
   do h = 1, neltype
      do i = offset + 1, offset + eltypepartlens(h)
         do j = offset + 1, offset + eltypepartlens(h)
!            biasmat(j, i) = bias_scale**2*bias_ratio**(equivmat(j, i))
            biasmat(j, i) = bias_scale**2*(maxequiv - equivmat(j, i))
         end do
      end do
      offset = offset + eltypepartlens(h)
   end do

end subroutine

end module
