! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

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

module pruning_atoms
! Distance-based pruning of atom pairs for the cluster assignment
use parameters
use common_types
use sorting
use molecule
use linked_list_types
use flags
implicit none

! Pruning tolerance (Angstrom)
real(rk) :: prunetol
! Expected value factor = 2*sqrt(3)
real(rk), parameter :: EVALFAC = 3.4641
! Selected pruning method (prune_none or prune_rd)
procedure(prune_proc), pointer :: prune_procedure

abstract interface
   subroutine prune_proc( atomtypes, coords1, coords2, prunes)
      use parameters
      use common_types
      use molecule
      use linked_list_types
      type(partition_t), intent(in) :: atomtypes
      real(rk), dimension(:,:), intent(in) :: coords1, coords2
      type(bool_matrix), dimension(:), allocatable, intent(out) :: prunes
   end subroutine
end interface

contains

subroutine prune_none( atomtypes, coords1, coords2, prunes)
! No pair is pruned
   type(partition_t), intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(bool_matrix), dimension(:), allocatable, intent(out) :: prunes
   ! Local variables
   integer(ik) :: h, i, j
   integer(ik) :: n_items1, n_items2

   allocate (prunes(atomtypes%n_parts))
   do h = 1, atomtypes%n_parts
      n_items1 = atomtypes%parts(h)%n_items1
      n_items2 = atomtypes%parts(h)%n_items2
      allocate (prunes(h)%a(n_items1, n_items2))
      do i = 1, n_items1
         do j = 1, n_items2
            prunes(h)%a(j, i) = .FALSE.
         end do
      end do
   end do

end subroutine

subroutine prune_rd( atomtypes, coords1, coords2, prunes)
! Prune the pairs whose atoms have incompatible environments: atoms i and j
! are not paired if, for some atom type, their sorted distances to the
! atoms of that type differ by more than EVALFAC*prunetol. The test is
! orientation independent, so it holds for all trials.
   type(partition_t), intent(in) :: atomtypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   type(bool_matrix), dimension(:), allocatable, intent(out) :: prunes
   ! Local variables
   type(real_listlist), allocatable, dimension(:) :: dists1, dists2
   integer(ik) :: n_items1, n_items2
   integer(ik) :: h, i, j, k, iatom, jatom

   allocate (dists1(size(coords1, dim=2)))
   allocate (dists2(size(coords2, dim=2)))
   allocate (prunes(atomtypes%n_parts))

   do i = 1, size(coords1, dim=2)
      allocate (dists1(i)%u(atomtypes%n_parts))
      allocate (dists2(i)%u(atomtypes%n_parts))
      do h = 1, atomtypes%n_parts
         allocate (dists1(i)%u(h)%u(atomtypes%parts(h)%n_items1))
         allocate (dists2(i)%u(h)%u(atomtypes%parts(h)%n_items2))
      end do
   end do

   do i = 1, size(coords1, dim=2)
      do h = 1, atomtypes%n_parts
         do j = 1, atomtypes%parts(h)%n_items1
            jatom = atomtypes%parts(h)%items1(j)
            dists1(i)%u(h)%u(j) = sqrt(sum((coords1(:, jatom) - coords1(:, i))**2))
         end do
         call quicksort(dists1(i)%u(h)%u)
      end do
   end do

   do i = 1, size(coords2, dim=2)
      do h = 1, atomtypes%n_parts
         do j = 1, atomtypes%parts(h)%n_items2
            jatom = atomtypes%parts(h)%items2(j)
            dists2(i)%u(h)%u(j) = sqrt(sum((coords2(:, jatom) - coords2(:, i))**2))
         end do
         call quicksort(dists2(i)%u(h)%u)
      end do
   end do

   do h = 1, atomtypes%n_parts
      n_items1 = atomtypes%parts(h)%n_items1
      n_items2 = atomtypes%parts(h)%n_items2
      allocate (prunes(h)%a(n_items1, n_items2))
      prunes(h)%a = .FALSE.
      do i = 1, n_items1
         iatom = atomtypes%parts(h)%items1(i)
         do j = 1, n_items2
            jatom = atomtypes%parts(h)%items2(j)
            do k = 1, atomtypes%n_parts
               if (any(abs(dists2(jatom)%u(k)%u - dists1(iatom)%u(k)%u) > EVALFAC*prunetol)) then
                  prunes(h)%a(j, i) = .TRUE.
                  exit
               end if
            end do
         end do
      end do
   end do

end subroutine

end module
