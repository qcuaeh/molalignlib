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
use derived_types
use molecule
use options
use sorting
use lcrs_tree

implicit none

real(rk) :: prune_tol
procedure(prune_proc), pointer :: prune_procedure

abstract interface
   subroutine prune_proc( atomtypes, atoms1, atoms2, prunes)
      use parameters
      use derived_types
      use molecule
      use lcrs_tree
      type(partition_t), intent(in) :: atomtypes
      type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
      type(bool_matrix), dimension(:), allocatable, intent(out) :: prunes
   end subroutine
end interface

contains

subroutine prune_none( atomtypes, atoms1, atoms2, prunes)
   type(partition_t), intent(in) :: atomtypes
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(bool_matrix), dimension(:), allocatable, intent(out) :: prunes
   ! Local variables
   integer :: h, i, j
   integer :: num_items1, num_items2

   allocate (prunes(atomtypes%num_parts))
   do h = 1, atomtypes%num_parts
      num_items1 = atomtypes%parts(h)%num_items1
      num_items2 = atomtypes%parts(h)%num_items2
      allocate (prunes(h)%ee(num_items1, num_items2))
      do i = 1, num_items1
         do j = 1, num_items2
            prunes(h)%ee(j, i) = .false.
         end do
      end do
   end do

end subroutine

subroutine prune_rd( atomtypes, atoms1, atoms2, prunes)
   type(partition_t), intent(in) :: atomtypes
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(bool_matrix), dimension(:), allocatable, intent(out) :: prunes
   ! Local variables
   type(real_listlist), allocatable, dimension(:) :: dists1, dists2
   real(rk), dimension(:,:), allocatable :: coords1, coords2
   integer :: num_items1, num_items2
   integer :: h, i, j, k, iatom, jatom

   coords1 = get_coords(atoms1)
   coords2 = get_coords(atoms2)
   allocate (dists1(size(coords1, dim=2)))
   allocate (dists2(size(coords2, dim=2)))
   allocate (prunes(atomtypes%num_parts))

   do i = 1, size(coords1, dim=2)
      allocate (dists1(i)%e(atomtypes%num_parts))
      allocate (dists2(i)%e(atomtypes%num_parts))
      do h = 1, atomtypes%num_parts
         allocate (dists1(i)%e(h)%e(atomtypes%parts(h)%num_items1))
         allocate (dists2(i)%e(h)%e(atomtypes%parts(h)%num_items2))
      end do
   end do

   do i = 1, size(coords1, dim=2)
      do h = 1, atomtypes%num_parts
         do j = 1, atomtypes%parts(h)%num_items1
            jatom = atomtypes%parts(h)%items1(j)
            dists1(i)%e(h)%e(j) = sqrt(sum((coords1(:, jatom) - coords1(:, i))**2))
         end do
         call sort(dists1(i)%e(h)%e)
      end do
   end do

   do i = 1, size(coords2, dim=2)
      do h = 1, atomtypes%num_parts
         do j = 1, atomtypes%parts(h)%num_items2
            jatom = atomtypes%parts(h)%items2(j)
            dists2(i)%e(h)%e(j) = sqrt(sum((coords2(:, jatom) - coords2(:, i))**2))
         end do
         call sort(dists2(i)%e(h)%e)
      end do
   end do

   do h = 1, atomtypes%num_parts
      num_items1 = atomtypes%parts(h)%num_items1
      num_items2 = atomtypes%parts(h)%num_items2
      allocate (prunes(h)%ee(num_items1, num_items2))
      prunes(h)%ee = .false.
      do i = 1, num_items1
         iatom = atomtypes%parts(h)%items1(i)
         do j = 1, num_items2
            jatom = atomtypes%parts(h)%items2(j)
            do k = 1, atomtypes%num_parts
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
