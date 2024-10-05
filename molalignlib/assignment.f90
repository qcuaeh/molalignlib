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
use partition
use lap_dense
use lap_sparse

implicit none

private
public assign_atoms
public assign_atoms_nearest
public assign_atoms_pruned
public assign_atoms_biased
public f_assign

abstract interface
   subroutine f_assign( eltypes0, coords0, coords1, prunemat, biasmat, mapping)
      use kinds
      use partition
      type(atompartition_type), intent(in) :: eltypes0
      real(rk), dimension(:, :), intent(in) :: coords0, coords1
      logical, dimension(:, :), intent(in) :: prunemat
      real(rk), dimension(:, :), intent(in) :: biasmat
      integer, dimension(:), intent(out) :: mapping
   end subroutine
end interface

procedure(f_assign), pointer :: assign_atoms

contains

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms_nearest( eltypes0, coords0, coords1, prunemat, biasmat, mapping)
   type(atompartition_type), intent(in) :: eltypes0
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   logical, dimension(:, :), intent(in) :: prunemat
   real(rk), dimension(:, :), intent(in) :: biasmat
   integer, dimension(:), intent(out) :: mapping
   ! Local variables
   integer :: h, offset
   integer, allocatable :: eltypesubsetlens(:)
   real(rk) :: dummy

   allocate (eltypesubsetlens(size(eltypes0%subsets)))
   do h = 1, size(eltypes0%subsets)
      eltypesubsetlens(h) = size(eltypes0%subsets(h)%atomidcs)
   end do

   ! Fill distance matrix for each block

   offset = 0
   do h = 1, size(eltypes0%subsets)
      call minperm_nearest(eltypesubsetlens(h), coords0(:, offset+1:offset+eltypesubsetlens(h)), &
         coords1(:, offset+1:offset+eltypesubsetlens(h)), prunemat(offset+1:offset+eltypesubsetlens(h), &
         offset+1:offset+eltypesubsetlens(h)), mapping(offset+1:), dummy)
      mapping(offset+1:offset+eltypesubsetlens(h)) = mapping(offset+1:offset+eltypesubsetlens(h)) + offset
      offset = offset + eltypesubsetlens(h)
   end do

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms_pruned( eltypes0, coords0, coords1, prunemat, biasmat, mapping)
   type(atompartition_type), intent(in) :: eltypes0
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   logical, dimension(:, :), intent(in) :: prunemat
   real(rk), dimension(:, :), intent(in) :: biasmat
   integer, dimension(:), intent(out) :: mapping
   ! Local variables
   integer :: h, offset
   integer, allocatable :: eltypesubsetlens(:)
   real(rk) :: dummy

   allocate (eltypesubsetlens(size(eltypes0%subsets)))
   do h = 1, size(eltypes0%subsets)
      eltypesubsetlens(h) = size(eltypes0%subsets(h)%atomidcs)
   end do

   ! Fill distance matrix for each block

   offset = 0

   do h = 1, size(eltypes0%subsets)
      call minperm_pruned(eltypesubsetlens(h), coords0(:, offset+1:offset+eltypesubsetlens(h)), &
         coords1(:, offset+1:offset+eltypesubsetlens(h)), prunemat(offset+1:offset+eltypesubsetlens(h), &
         offset+1:offset+eltypesubsetlens(h)), mapping(offset+1:), dummy)
      mapping(offset+1:offset+eltypesubsetlens(h)) = mapping(offset+1:offset+eltypesubsetlens(h)) + offset
      offset = offset + eltypesubsetlens(h)
   end do

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms_biased( eltypes0, coords0, coords1, prunemat, biasmat, mapping)
   type(atompartition_type), intent(in) :: eltypes0
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   logical, dimension(:, :), intent(in) :: prunemat
   real(rk), dimension(:, :), intent(in) :: biasmat
   integer, dimension(:), intent(out) :: mapping
   ! Local variables
   integer :: h, i, j, offset
   real(rk), allocatable :: costs(:, :)
   integer, allocatable :: eltypesubsetlens(:)
   real(rk) :: dummy

   allocate (eltypesubsetlens(size(eltypes0%subsets)))
   do h = 1, size(eltypes0%subsets)
      eltypesubsetlens(h) = size(eltypes0%subsets(h)%atomidcs)
   end do

   allocate (costs(size(coords0, dim=2), size(coords0, dim=2)))

   offset = 0
   do h = 1, size(eltypes0%subsets)
      do i = offset + 1, offset + eltypesubsetlens(h)
         do j = offset + 1, offset + eltypesubsetlens(h)
            costs(i, j) = sum((coords0(:, i) - coords1(:, j))**2) + biasmat(j, i)
         end do
      end do
      call assndx(1, costs(offset+1:, offset+1:), eltypesubsetlens(h), eltypesubsetlens(h), mapping(offset+1:), dummy)
      mapping(offset+1:offset+eltypesubsetlens(h)) = mapping(offset+1:offset+eltypesubsetlens(h)) + offset
      offset = offset + eltypesubsetlens(h)
   end do

end subroutine

end module
