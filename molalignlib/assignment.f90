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

module assignment
use stdio
use kinds
use bounds
use lapkm
use lapjv

implicit none

private
public assign_atoms_function
public assign_atoms_fast
public assign_atoms_pruned
public assign_atoms_biased
public f_assign

abstract interface
   subroutine f_assign(natom, neltype, eltypepartlens, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
      nadjmna1, adjmnalen1, adjmnalist1, coords1, weights, equivmat, prunemat, biasmat, mapping)
      use kinds
      integer, intent(in) :: natom, neltype
      integer, dimension(:), intent(in) :: eltypepartlens
      integer, dimension(:, :), intent(in) :: nadjmna0, nadjmna1
      integer, dimension(:, :, :), intent(in) :: adjmnalen0, adjmnalen1
      integer, dimension(:, :, :), intent(in) :: adjmnalist0, adjmnalist1
      real(rk), dimension(:, :), intent(in) :: coords0, coords1
      real(rk), dimension(:), intent(in) :: weights
      integer, dimension(:, :), intent(in) :: equivmat
      logical, dimension(:, :), intent(in) :: prunemat
      real(rk), dimension(:, :), intent(in) :: biasmat
      integer, dimension(:), intent(out) :: mapping
   end subroutine
end interface

procedure(f_assign), pointer :: assign_atoms_function

contains

subroutine assign_atoms_fast(natom, neltype, eltypepartlens, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
   nadjmna1, adjmnalen1, adjmnalist1, coords1, weights, equivmat, prunemat, biasmat, mapping)
! Find best correspondence between points sets with fixed orientation
   integer, intent(in) :: natom, neltype
   integer, dimension(:), intent(in) :: eltypepartlens
   integer, dimension(:, :), intent(in) :: nadjmna0, nadjmna1
   integer, dimension(:, :, :), intent(in) :: adjmnalen0, adjmnalen1
   integer, dimension(:, :, :), intent(in) :: adjmnalist0, adjmnalist1
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   real(rk), dimension(:), intent(in) :: weights
   integer, dimension(:, :), intent(in) :: equivmat
   logical, dimension(:, :), intent(in) :: prunemat
   real(rk), dimension(:, :), intent(in) :: biasmat
   integer, dimension(:), intent(out) :: mapping

   integer :: h, offset
   real(rk) :: dummy

   ! Fill distance matrix for each block

   offset = 0

   do h = 1, neltype
      call minperm_fast(eltypepartlens(h), coords0(:, offset+1:offset+eltypepartlens(h)), &
         coords1(:, offset+1:offset+eltypepartlens(h)), prunemat(offset+1:offset+eltypepartlens(h), &
         offset+1:offset+eltypepartlens(h)), mapping(offset+1:), dummy)
      mapping(offset+1:offset+eltypepartlens(h)) = mapping(offset+1:offset+eltypepartlens(h)) + offset
      offset = offset + eltypepartlens(h)
   end do

end subroutine

subroutine assign_atoms_pruned(natom, neltype, eltypepartlens, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
   nadjmna1, adjmnalen1, adjmnalist1, coords1, weights, equivmat, prunemat, biasmat, mapping)
! Find best correspondence between points sets with fixed orientation
   integer, intent(in) :: natom, neltype
   integer, dimension(:), intent(in) :: eltypepartlens
   integer, dimension(:, :), intent(in) :: nadjmna0, nadjmna1
   integer, dimension(:, :, :), intent(in) :: adjmnalen0, adjmnalen1
   integer, dimension(:, :, :), intent(in) :: adjmnalist0, adjmnalist1
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   real(rk), dimension(:), intent(in) :: weights
   integer, dimension(:, :), intent(in) :: equivmat
   logical, dimension(:, :), intent(in) :: prunemat
   real(rk), dimension(:, :), intent(in) :: biasmat
   integer, dimension(:), intent(out) :: mapping

   integer :: h, offset
   real(rk) :: dummy

   ! Fill distance matrix for each block

   offset = 0

   do h = 1, neltype
      call minperm_pruned(eltypepartlens(h), coords0(:, offset+1:offset+eltypepartlens(h)), &
         coords1(:, offset+1:offset+eltypepartlens(h)), prunemat(offset+1:offset+eltypepartlens(h), &
         offset+1:offset+eltypepartlens(h)), mapping(offset+1:), dummy)
      mapping(offset+1:offset+eltypepartlens(h)) = mapping(offset+1:offset+eltypepartlens(h)) + offset
      offset = offset + eltypepartlens(h)
   end do

end subroutine

subroutine assign_atoms_biased(natom, neltype, eltypepartlens, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
   nadjmna1, adjmnalen1, adjmnalist1, coords1, weights, equivmat, prunemat, biasmat, mapping)
! Find best correspondence between points sets with fixed orientation
   integer, intent(in) :: natom, neltype
   integer, dimension(:), intent(in) :: eltypepartlens
   integer, dimension(:, :), intent(in) :: nadjmna0, nadjmna1
   integer, dimension(:, :, :), intent(in) :: adjmnalen0, adjmnalen1
   integer, dimension(:, :, :), intent(in) :: adjmnalist0, adjmnalist1
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   real(rk), dimension(:), intent(in) :: weights
   integer, dimension(:, :), intent(in) :: equivmat
   logical, dimension(:, :), intent(in) :: prunemat
   real(rk), dimension(:, :), intent(in) :: biasmat
   integer, dimension(:), intent(out) :: mapping

   integer :: h, i, j, offset
!   real(rk) :: costs(natom, natom) ! Causes allocation errors
   real(rk), allocatable :: costs(:, :)
   real(rk) :: dummy

   allocate(costs(natom, natom))

   offset = 0
   do h = 1, neltype
      do i = offset + 1, offset + eltypepartlens(h)
         do j = offset + 1, offset + eltypepartlens(h)
!            costs(i, j) = sum((coords0(:, i) - coords1(:, j))**2)
            costs(i, j) = biasmat(j, i) + sum((coords0(:, i) - coords1(:, j))**2)
!            costs(i, j) = biasmat(j, i) &
!                        + equivdist(natom, i, j, equivmat(j, i), nadjmna0, adjmnalen0, adjmnalist0, &
!                             coords0, nadjmna1, adjmnalen1, adjmnalist1, coords1, weights)
         end do
      end do
      call assndx(1, costs(offset+1:, offset+1:), eltypepartlens(h), eltypepartlens(h), mapping(offset+1:), dummy)
      mapping(offset+1:offset+eltypepartlens(h)) = mapping(offset+1:offset+eltypepartlens(h)) + offset
      offset = offset + eltypepartlens(h)
   end do

end subroutine

real(rk) function equivdist(natom, i, j, maxlevel, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
      nadjmna1, adjmnalen1, adjmnalist1, coords1, weights)
   integer, intent(in) :: natom, i, j, maxlevel
   integer, dimension(:, :), intent(in) :: nadjmna0, nadjmna1
   integer, dimension(:, :, :), intent(in) :: adjmnalen0, adjmnalen1
   integer, dimension(:, :, :), intent(in) :: adjmnalist0, adjmnalist1
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   real(rk), dimension(:), intent(in) :: weights

   real(rk) :: totdist, totweight
   logical, dimension(natom) :: mapped0, mapped1

   totdist = 0
   totweight = 0
   mapped0(:) = .false.
   mapped1(:) = .false.

   call recursivemap(i, j, 1, maxlevel, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
      nadjmna1, adjmnalen1, adjmnalist1, coords1, weights, totdist, totweight, mapped0, mapped1)
   equivdist = totdist/totweight

end function

recursive subroutine recursivemap(i, j, level, maxlevel, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
      nadjmna1, adjmnalen1, adjmnalist1, coords1, weights, totdist, totweight, mapped0, mapped1)
   integer, intent(in) :: i, j, level, maxlevel
   integer, dimension(:, :), intent(in) :: nadjmna0, nadjmna1
   integer, dimension(:, :, :), intent(in) :: adjmnalen0, adjmnalen1
   integer, dimension(:, :, :), intent(in) :: adjmnalist0, adjmnalist1
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   real(rk), dimension(:), intent(in) :: weights
   real(rk), intent(out) :: totdist, totweight
   logical, dimension(:), intent(inout) :: mapped0, mapped1

   integer :: h, k, l, m, n, offset
   integer, dimension(maxcoord) :: mapping, idx0, idx1
   real(rk) :: distmat(maxcoord, maxcoord)
   real(rk) :: dummy

!   print *, '>>>', i, j, ':', level, maxlevel - level
!!   print *, i, ':', adjmnalist0(:4, i, maxlevel - level)
!!   print *, j, ':', adjmnalist1(:4, j, maxlevel - level)
   mapped0(i) = .true.
   mapped1(j) = .true.
   totdist = totdist + weights(i)*sum((coords0(:, i) - coords1(:, j))**2)
   totweight = totweight + weights(i)
   if (level < maxlevel) then
      offset = 0
!      print *, 'num:', nadjmna0(i, maxlevel - level), &
!         nadjmna0(i, maxlevel - level) == nadjmna1(j, maxlevel - level)
      do h = 1, nadjmna0(i, maxlevel - level)
!         print *, 'len:', adjmnalen0(h, i, maxlevel - level), &
!            adjmnalen0(h, i, maxlevel - level) == adjmnalen1(h, j, maxlevel - level)
         m = 0
         do k = offset + 1, offset + adjmnalen0(h, i, maxlevel - level)
            if (.not. mapped0(adjmnalist0(k, i, maxlevel - level))) then
!!               print *, 'not mapped0:', adjmnalist0(k, i, maxlevel - level)
               m = m + 1
               idx0(m) = adjmnalist0(k, i, maxlevel - level)
               n = 0
               do l = offset + 1, offset + adjmnalen1(h, j, maxlevel - level)
                  if (.not. mapped1(adjmnalist1(l, j, maxlevel - level))) then
!!                     print *, 'not mapped1:', adjmnalist1(l, j, maxlevel - level)
                     n = n + 1
                     idx1(n) = adjmnalist1(l, j, maxlevel - level)
                     distmat(m, n) = sum((coords0(:, idx0(m)) - coords1(:, idx1(n)))**2)
!!                  else
!!                     print *, 'mapped1:', adjmnalist1(l, j, maxlevel - level)
                  end if
               end do
!!            else
!!               print *, 'mapped0:', adjmnalist0(k, i, maxlevel - level)
            end if
         end do
         if (m > 0) then
            call assndx(1, distmat, m, m, mapping, dummy)
            do n = 1, m
               call recursivemap(idx0(n), idx1(mapping(n)), level + 1, maxlevel, &
                  nadjmna0, adjmnalen0, adjmnalist0, coords0, nadjmna1, adjmnalen1, adjmnalist1, &
                  coords1, weights, totdist, totweight, mapped0, mapped1)
            end do
         end if
         offset = offset + adjmnalen0(h, i, maxlevel - level)
      end do
   end if
!   print *, '<<<'

end subroutine

end module
