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
use kinds
use flags
use bounds
use sorting
use partition

implicit none

type dists_type
   real(rk), allocatable :: dists(:)
end type

type subsets_type
   type(dists_type), allocatable :: subsets(:)
end type

abstract interface
   subroutine f_prune( eltypes0, eltypes1, coords0, coords1, prunemat)
      use kinds
      use partition
      type(atompartition_type), intent(in) :: eltypes0, eltypes1
      real(rk), dimension(:, :), intent(in) :: coords0, coords1
      logical, dimension(:, :), intent(out) :: prunemat
   end subroutine
end interface

real(rk) :: prune_tol
procedure(f_prune), pointer :: prune_procedure

contains

subroutine prune_none( eltypes0, eltypes1, coords0, coords1, prunemat)
   type(atompartition_type), intent(in) :: eltypes0, eltypes1
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   logical, dimension(:, :), intent(out) :: prunemat
   ! Local variables
   integer :: h, i, j
   integer :: atomidx_i, atomidx_j

   do h = 1, size(eltypes0%subsets)
      do i = 1, size(eltypes0%subsets(h)%atomidcs)
         atomidx_i = eltypes0%subsets(h)%atomidcs(i)
         do j = 1, size(eltypes1%subsets(h)%atomidcs)
            atomidx_j = eltypes1%subsets(h)%atomidcs(j)
            prunemat(atomidx_j, atomidx_i) = .true.
         end do
      end do
   end do

end subroutine

subroutine prune_rd( eltypes0, eltypes1, coords0, coords1, prunemat)
   type(atompartition_type), intent(in) :: eltypes0, eltypes1
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   logical, dimension(:, :), intent(out) :: prunemat
   ! Local variables
   integer :: h, i, j, k
   integer :: atomidx_i, atomidx_j
   type(subsets_type), allocatable, dimension(:) :: d0, d1
   logical :: pruned

   allocate (d0(size(coords0, dim=2)))
   allocate (d1(size(coords1, dim=2)))
   do i = 1, size(coords0, dim=2)
      allocate (d0(i)%subsets(size(eltypes0%subsets)))
      allocate (d1(i)%subsets(size(eltypes1%subsets)))
      do h = 1, size(eltypes0%subsets)
         allocate (d0(i)%subsets(h)%dists(size(eltypes0%subsets(h)%atomidcs)))
         allocate (d1(i)%subsets(h)%dists(size(eltypes1%subsets(h)%atomidcs)))
      end do
   end do

   do i = 1, size(coords0, dim=2)
      do h = 1, size(eltypes0%subsets)
         do j = 1, size(eltypes0%subsets(h)%atomidcs)
            atomidx_j = eltypes0%subsets(h)%atomidcs(j)
            d0(i)%subsets(h)%dists(j) = sqrt(sum((coords0(:, atomidx_j) - coords0(:, i))**2))
         end do
         call sort(d0(i)%subsets(h)%dists)
      end do
   end do

   do i = 1, size(coords1, dim=2)
      do h = 1, size(eltypes1%subsets)
         do j = 1, size(eltypes1%subsets(h)%atomidcs)
            atomidx_j = eltypes1%subsets(h)%atomidcs(j)
            d1(i)%subsets(h)%dists(j) = sqrt(sum((coords1(:, atomidx_j) - coords1(:, i))**2))
         end do
         call sort(d1(i)%subsets(h)%dists)
      end do
   end do

   do h = 1, size(eltypes0%subsets)
      do i = 1, size(eltypes0%subsets(h)%atomidcs)
         atomidx_i = eltypes0%subsets(h)%atomidcs(i)
         do j = 1, size(eltypes1%subsets(h)%atomidcs)
            atomidx_j = eltypes1%subsets(h)%atomidcs(j)
            pruned = .false.
            do k = 1, size(eltypes0%subsets)
               if (any(abs(d1(atomidx_j)%subsets(k)%dists - d0(atomidx_i)%subsets(k)%dists) > prune_tol)) then
                  pruned = .true.
                  exit
               end if
            end do
            prunemat(atomidx_j, atomidx_i) = .not. pruned
         end do
      end do
   end do

end subroutine

end module
