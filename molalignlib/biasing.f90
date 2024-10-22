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
use partition
use assorting

implicit none

abstract interface
   subroutine f_bias( eltypes0, eltypes1, adjlists0, adjlists1, biasmat)
      use kinds
      use partition
      type(partition_type), intent(in) :: eltypes0, eltypes1
      type(indexlist_type), allocatable, dimension(:), intent(in) :: adjlists0, adjlists1
      real(rk), intent(out) :: biasmat(:, :)
   end subroutine
end interface

real(rk) :: bias_scale
procedure(f_bias), pointer :: bias_procedure

contains

subroutine bias_none( eltypes0, eltypes1, adjlists0, adjlists1, biasmat)
   type(partition_type), intent(in) :: eltypes0, eltypes1
   type(indexlist_type), allocatable, dimension(:), intent(in) :: adjlists0, adjlists1
   real(rk), intent(out) :: biasmat(:, :)
   ! Local variables
   integer :: h, i, j
   integer :: atomidx_i, atomidx_j

   do h = 1, eltypes0%size
      do i = 1, eltypes0%parts(h)%size
         atomidx_i = eltypes0%parts(h)%indices(i)
         do j = 1, eltypes1%parts(h)%size
            atomidx_j = eltypes1%parts(h)%indices(j)
            biasmat(atomidx_j, atomidx_i) = 0
         end do
      end do
   end do

end subroutine

subroutine bias_mna( eltypes0, eltypes1, adjlists0, adjlists1, biasmat)
   type(partition_type), intent(in) :: eltypes0, eltypes1
   type(indexlist_type), allocatable, dimension(:), intent(in) :: adjlists0, adjlists1
   real(rk), intent(out) :: biasmat(:, :)
   ! Local variables
   integer :: h, i, j, maxequiv
   integer :: atomidx_i, atomidx_j
   integer, allocatable :: equivmat(:, :)

   ! Calculate MNA equivalence matrix
   call compute_equivmat(eltypes0, eltypes1, adjlists0, adjlists1, equivmat)

   maxequiv = 0
   do h = 1, eltypes0%size
      do i = 1, eltypes0%parts(h)%size
         atomidx_i = eltypes0%parts(h)%indices(i)
         do j = 1, eltypes1%parts(h)%size
            atomidx_j = eltypes1%parts(h)%indices(j)
            maxequiv = max(maxequiv, equivmat(atomidx_j, atomidx_i))
         end do
      end do
   end do

   do h = 1, eltypes0%size
      do i = 1, eltypes0%parts(h)%size
         atomidx_i = eltypes0%parts(h)%indices(i)
         do j = 1, eltypes1%parts(h)%size
            atomidx_j = eltypes1%parts(h)%indices(j)
            biasmat(atomidx_j, atomidx_i) = bias_scale**2*(maxequiv - equivmat(atomidx_j, atomidx_i))
         end do
      end do
   end do

end subroutine

! Calculate the maximum common MNA level for all atom cross assignments
subroutine compute_equivmat( eltypes0, eltypes1, adjlists0, adjlists1, equivmat)
   type(partition_type), intent(in) :: eltypes0, eltypes1
   type(indexlist_type), allocatable, dimension(:), intent(in) :: adjlists0, adjlists1
   integer, allocatable, dimension(:, :), intent(out) :: equivmat
   ! Local variables
   integer :: h, i, j
   integer :: atomidx_i, atomidx_j
   integer :: ntype, nintype, level
   integer, allocatable, dimension(:) :: types0, types1
   integer, allocatable, dimension(:) :: intypes0, intypes1

   allocate (equivmat(eltypes1%total_size, eltypes0%total_size))
   allocate (types0(eltypes0%total_size), types1(eltypes1%total_size))
   allocate (intypes0(eltypes0%total_size), intypes1(eltypes1%total_size))

   level = 1
   nintype = eltypes0%size
   intypes0 = eltypes0%index_part_map
   intypes1 = eltypes1%index_part_map

   do

      do h = 1, eltypes0%size
         do i = 1, eltypes0%parts(h)%size
            atomidx_i = eltypes0%parts(h)%indices(i)
            do j = 1, eltypes1%parts(h)%size
               atomidx_j = eltypes1%parts(h)%indices(j)
               if (intypes0(atomidx_i) == intypes1(atomidx_j)) then
                  equivmat(atomidx_j, atomidx_i) = level
               end if
            end do
         end do
      end do

      call compute_crossmnatypes(adjlists0, adjlists1, nintype, intypes0, intypes1, ntype, types0, types1)

      if (all(types0 == intypes0) .and. all(types1 == intypes1)) exit

      nintype = ntype
      intypes0 = types0
      intypes1 = types1
      level = level + 1

   end do

end subroutine

end module
