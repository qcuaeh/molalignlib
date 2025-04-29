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

module pruning
use parameters
use basetypes
use molecule
use globals
use sorting
use lcrs_tree

implicit none

real(rk) :: prune_tol
procedure(prune_proc), pointer :: prune_procedure

abstract interface
   subroutine prune_proc( eltypes, mol1, mol2, prunes)
      use parameters
      use basetypes
      use molecule
      use lcrs_tree
      type(item_partition), intent(in) :: eltypes
      type(mol_type), intent(in) :: mol1, mol2
      type(bool_matrix), dimension(:), allocatable, intent(out) :: prunes
   end subroutine
end interface

contains

subroutine prune_none( eltypes, mol1, mol2, prunes)
   type(item_partition), intent(in) :: eltypes
   type(mol_type), intent(in) :: mol1, mol2
   type(bool_matrix), dimension(:), allocatable, intent(out) :: prunes
   ! Local variables
   integer :: h, i, j
   integer :: num_items1, num_items2

   allocate (prunes(eltypes%num_parts))
   do h = 1, eltypes%num_parts
      num_items1 = eltypes%parts(h)%num_items1
      num_items2 = eltypes%parts(h)%num_items2
      allocate (prunes(h)%ee(num_items1, num_items2))
      do i = 1, num_items1
         do j = 1, num_items2
            prunes(h)%ee(j, i) = .false.
         end do
      end do
   end do

end subroutine

subroutine prune_rd( eltypes, mol1, mol2, prunes)
   type(item_partition), intent(in) :: eltypes
   type(mol_type), intent(in) :: mol1, mol2
   type(bool_matrix), dimension(:), allocatable, intent(out) :: prunes
   ! Local variables
   type(real_listlist), allocatable, dimension(:) :: dists1, dists2
   real(rk), dimension(:,:), allocatable :: coords1, coords2
   integer :: num_items1, num_items2
   integer :: h, i, j, k, iatom, jatom

   coords1 = get_coords(mol1)
   coords2 = get_coords(mol2)
   allocate (dists1(size(coords1, dim=2)))
   allocate (dists2(size(coords2, dim=2)))
   allocate (prunes(eltypes%num_parts))

   do i = 1, size(coords1, dim=2)
      allocate (dists1(i)%e(eltypes%num_parts))
      allocate (dists2(i)%e(eltypes%num_parts))
      do h = 1, eltypes%num_parts
         allocate (dists1(i)%e(h)%e(eltypes%parts(h)%num_items1))
         allocate (dists2(i)%e(h)%e(eltypes%parts(h)%num_items2))
      end do
   end do

   do i = 1, size(coords1, dim=2)
      do h = 1, eltypes%num_parts
         do j = 1, eltypes%parts(h)%num_items1
            jatom = eltypes%parts(h)%items1(j)
            dists1(i)%e(h)%e(j) = sqrt(sum((coords1(:, jatom) - coords1(:, i))**2))
         end do
         call sort(dists1(i)%e(h)%e)
      end do
   end do

   do i = 1, size(coords2, dim=2)
      do h = 1, eltypes%num_parts
         do j = 1, eltypes%parts(h)%num_items2
            jatom = eltypes%parts(h)%items2(j)
            dists2(i)%e(h)%e(j) = sqrt(sum((coords2(:, jatom) - coords2(:, i))**2))
         end do
         call sort(dists2(i)%e(h)%e)
      end do
   end do

   do h = 1, eltypes%num_parts
      num_items1 = eltypes%parts(h)%num_items1
      num_items2 = eltypes%parts(h)%num_items2
      allocate (prunes(h)%ee(num_items1, num_items2))
      prunes(h)%ee = .false.
      do i = 1, num_items1
         iatom = eltypes%parts(h)%items1(i)
         do j = 1, num_items2
            jatom = eltypes%parts(h)%items2(j)
            do k = 1, eltypes%num_parts
               if (any(abs(dists2(jatom)%e(k)%e - dists1(iatom)%e(k)%e) > prune_tol)) then
                  prunes(h)%ee(j, i) = .true.
                  exit
               end if
            end do
         end do
      end do
   end do

end subroutine

end module
