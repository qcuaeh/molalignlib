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
use permutation
use lapdense
use lapsparse

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
   integer :: h
   integer, allocatable :: atomidcs(:)
   integer, allocatable :: auxmap(:)
   real(rk) :: dummy

   allocate (auxmap(eltypes0%largest))

   ! Fill distance matrix for each block

   do h = 1, size(eltypes0%subsets)
      atomidcs = eltypes0%subsets(h)%atomidcs 
      call minperm_nearest(size(atomidcs), coords0(:, atomidcs), coords1(:, atomidcs), auxmap, dummy)
      mapping(atomidcs) = atomidcs(auxmap(:size(atomidcs)))
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
   integer :: h
   integer, allocatable :: atomidcs(:)
   integer, allocatable :: auxmap(:)
   real(rk) :: dummy

   allocate (auxmap(eltypes0%largest))

   ! Fill distance matrix for each block

   do h = 1, size(eltypes0%subsets)
      atomidcs = eltypes0%subsets(h)%atomidcs 
      call minperm_pruned(size(atomidcs), coords0(:, atomidcs), coords1(:, atomidcs), &
           prunemat(atomidcs, atomidcs), auxmap, dummy)
      mapping(atomidcs) = atomidcs(auxmap(:size(atomidcs)))
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
   integer :: h
   integer, allocatable :: atomidcs(:)
   integer, allocatable :: auxmap(:)

   allocate (auxmap(eltypes0%largest))

   ! Fill distance matrix for each block

   do h = 1, size(eltypes0%subsets)
      atomidcs = eltypes0%subsets(h)%atomidcs 
      call minperm_biased(coords0(:, atomidcs), coords1(:, atomidcs), biasmat(atomidcs, atomidcs), auxmap)
      mapping(atomidcs) = atomidcs(auxmap(:size(atomidcs)))
   end do

end subroutine

end module
