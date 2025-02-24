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

module lap_driver
use parameters
use globals
use common_types
use lcrs_tree
use permutation
use lap_solvers

implicit none

private
public assign_atoms
public assign_atoms_biased
public assign_atoms_nearest
public assign_atoms_pruned

contains

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms_nearest( eltypes, coords1, coords2, pruned, mnadiffs, atomperm)
   type(bipartition_container), target, intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(boolmatrix_type), dimension(:), intent(in) :: pruned
   type(realmatrix_type), dimension(:), intent(in) :: mnadiffs
   integer, dimension(:), intent(out) :: atomperm
   ! Local variables
   integer :: h, num_items1
   integer, allocatable :: auxperm(:)
   integer, dimension(:), pointer :: indices1, indices2
   real(rk) :: dist

   allocate (auxperm(maxval(eltypes%parts%num_items1)))

   ! Fill distance matrix for each block

   do h = 1, eltypes%num_parts
      num_items1 = eltypes%parts(h)%num_items1
      indices1 => eltypes%parts(h)%indices1
      indices2 => eltypes%parts(h)%indices2
      call minperm_nearest(num_items1, indices1, indices2, coords1, coords2, auxperm, dist)
      atomperm(indices1) = indices2(auxperm(:num_items1))
   end do

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms_pruned( eltypes, coords1, coords2, pruned, atomperm)
   type(bipartition_container), target, intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(boolmatrix_type), dimension(:), intent(in) :: pruned
   integer, dimension(:), intent(out) :: atomperm
   ! Local variables
   integer :: h, num_items1
   integer, allocatable :: auxperm(:)
   integer, dimension(:), pointer :: indices1, indices2
   real(rk) :: dist

   allocate (auxperm(maxval(eltypes%parts%num_items1)))

   ! Optimize atomperm for each block
   do h = 1, eltypes%num_parts
      num_items1 = eltypes%parts(h)%num_items1
      indices1 => eltypes%parts(h)%indices1
      indices2 => eltypes%parts(h)%indices2
      call minperm_pruned(num_items1, indices1, indices2, coords1, coords2, pruned(h)%b, auxperm, dist)
      atomperm(indices1) = indices2(auxperm(:num_items1))
   end do

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms( eltypes, coords1, coords2, atomperm, dist)
   type(bipartition_container), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   integer, dimension(:), intent(out) :: atomperm
   real(rk), intent(out) :: dist
   ! Local variables
   integer :: h
   integer, allocatable :: auxperm(:)

   allocate (auxperm(maxval(eltypes%parts%num_items1)))

   ! Optimize atomperm for each block
   do h = 1, eltypes%num_parts
      call minperm(eltypes%parts(h), coords1, coords2, auxperm, dist)
      atomperm(eltypes%parts(h)%indices1) = eltypes%parts(h)%indices2(auxperm(:eltypes%parts(h)%num_items1))
   end do

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms_biased( eltypes, coords1, coords2, mnadiffs, atomperm)
   type(bipartition_container), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(intmatrix_type), dimension(:), intent(in) :: mnadiffs
   integer, dimension(:), intent(out) :: atomperm
   ! Local variables
   integer :: h
   integer, allocatable :: auxperm(:)
   real(rk) :: dist

   allocate (auxperm(maxval(eltypes%parts%num_items1)))

   ! Optimize atomperm for each block
   do h = 1, eltypes%num_parts
      call minperm_biased(eltypes%parts(h), coords1, coords2, mnadiffs(h)%n, auxperm, dist)
      atomperm(eltypes%parts(h)%indices1) = eltypes%parts(h)%indices2(auxperm(:eltypes%parts(h)%num_items1))
   end do

end subroutine

end module
