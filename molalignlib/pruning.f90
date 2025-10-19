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
use sorting
use molecule
use lcrs_trees
use options
implicit none

real(rk) :: prune_tol
procedure(prune_proc), pointer :: prune_procedure

abstract interface
   subroutine prune_proc( atomtypes, coords1, coords2, prunes)
      use parameters
      use derived_types
      use molecule
      use lcrs_trees
      type(partition_t), intent(in) :: atomtypes
      real(rk), dimension(:,:), intent(in) :: coords1, coords2
      type(bool_matrix), dimension(:), allocatable, intent(out) :: prunes
   end subroutine
end interface

contains

subroutine prune_none( atomtypes, coords1, coords2, prunes)
   type(partition_t), intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(bool_matrix), dimension(:), allocatable, intent(out) :: prunes
   ! Local variables
   integer :: h, i, j
   integer :: num_items1, num_items2

   allocate (prunes(atomtypes%num_parts))
   do h = 1, atomtypes%num_parts
      num_items1 = atomtypes%parts(h)%num_items1
      num_items2 = atomtypes%parts(h)%num_items2
      allocate (prunes(h)%a(num_items1, num_items2))
      do i = 1, num_items1
         do j = 1, num_items2
            prunes(h)%a(j, i) = .false.
         end do
      end do
   end do

end subroutine

subroutine prune_rd( atomtypes, coords1, coords2, prunes)
   type(partition_t), intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(bool_matrix), dimension(:), allocatable, intent(out) :: prunes
   ! Local variables
   type(real_listlist), allocatable, dimension(:) :: dists1, dists2
   integer :: num_items1, num_items2
   integer :: h, i, j, k, iatom, jatom

   allocate (dists1(size(coords1, dim=2)))
   allocate (dists2(size(coords2, dim=2)))
   allocate (prunes(atomtypes%num_parts))

   do i = 1, size(coords1, dim=2)
      allocate (dists1(i)%u(atomtypes%num_parts))
      allocate (dists2(i)%u(atomtypes%num_parts))
      do h = 1, atomtypes%num_parts
         allocate (dists1(i)%u(h)%u(atomtypes%parts(h)%num_items1))
         allocate (dists2(i)%u(h)%u(atomtypes%parts(h)%num_items2))
      end do
   end do

   do i = 1, size(coords1, dim=2)
      do h = 1, atomtypes%num_parts
         do j = 1, atomtypes%parts(h)%num_items1
            jatom = atomtypes%parts(h)%items1(j)
            dists1(i)%u(h)%u(j) = sqrt(sum((coords1(:, jatom) - coords1(:, i))**2))
         end do
         call quicksort(dists1(i)%u(h)%u)
      end do
   end do

   do i = 1, size(coords2, dim=2)
      do h = 1, atomtypes%num_parts
         do j = 1, atomtypes%parts(h)%num_items2
            jatom = atomtypes%parts(h)%items2(j)
            dists2(i)%u(h)%u(j) = sqrt(sum((coords2(:, jatom) - coords2(:, i))**2))
         end do
         call quicksort(dists2(i)%u(h)%u)
      end do
   end do

   do h = 1, atomtypes%num_parts
      num_items1 = atomtypes%parts(h)%num_items1
      num_items2 = atomtypes%parts(h)%num_items2
      allocate (prunes(h)%a(num_items1, num_items2))
      prunes(h)%a = .false.
      do i = 1, num_items1
         iatom = atomtypes%parts(h)%items1(i)
         do j = 1, num_items2
            jatom = atomtypes%parts(h)%items2(j)
            do k = 1, atomtypes%num_parts
               if (any(abs(dists2(jatom)%u(k)%u - dists1(iatom)%u(k)%u) > 3.4641*prune_tol)) then
                  prunes(h)%a(j, i) = .true.
                  exit
               end if
            end do
         end do
      end do
   end do

end subroutine

end module
