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

module writemol
use stdio
use types
use bounds
use strutils
use chemdata

implicit none

contains

subroutine writexyz(unit, mol)
   ! Arguments
   integer, intent(in) :: unit
   type(t_mol) :: mol
   ! Local varibles
   integer :: i, natom
   character(:), allocatable :: title
   integer, allocatable :: elnums(:)
   real(rk), allocatable :: coords(:, :)

   title = mol%get_title()
   natom = mol%get_natom()
   elnums = mol%get_elnums()
   coords = mol%get_coords()

   write (unit, '(i0)') natom
   write (unit, '(a)') title

   do i = 1, natom
      write (unit, '(a,3(2x,f12.6))') elsym(elnums(i)), coords(:, i)
   end do

end subroutine

subroutine writemol2(unit, mol)
   ! Arguments
   integer, intent(in) :: unit
   type(t_mol) :: mol
   ! Local varibles
   integer :: i, natom
   character(:), allocatable :: title
   type(t_bond), allocatable :: bonds(:)
   integer, allocatable :: elnums(:), nadjs(:)
   real(rk), allocatable :: coords(:, :)

   title = mol%get_title()
   nadjs = mol%get_nadjs()
   bonds = mol%get_bonds()
   natom = mol%get_natom()
   elnums = mol%get_elnums()
   coords = mol%get_coords()

   write (unit, '(a)') '@<TRIPOS>MOLECULE'
   if (title /= '') then
      write (unit, '(a)') title
   else
      write (unit, '(a)') 'Untitled'
   end if
   write (unit, '(5(i4,1x))') natom, size(bonds), 0, 0, 0
   write (unit, '(a)') 'SMALL'
   write (unit, '(a)') 'NO_CHARGES'
   write (unit, '(a)') '@<TRIPOS>ATOM'

   do i = 1, natom
      write (unit, '(i4,2x,a2,3(1x,f12.6),2x,a4,1x,i2,1x,a4,1x,f7.3)') &
         i, elsym(elnums(i)), coords(:, i), atomtype(elnums(i), nadjs(i)), 1, 'MOL1', 0.
   end do

   write (unit, '(a)') '@<TRIPOS>BOND'

   do i = 1, size(bonds)
      write (unit, '(i4,1x,2(1x,i4),1x,a2)') i, bonds(i)%atomidx1, bonds(i)%atomidx2, '1'
   end do

end subroutine

function atomtype(elnum, nadj)
   integer, intent(in) :: elnum, nadj
   integer :: hybnum
   character(4) :: atomtype

   select case (elnum)
   case (6)
      hybnum = nadj - 1
      atomtype = trim(elsym(elnum))//'.'//intstr(hybnum)
   case (7, 15)
      hybnum = min(nadj, 3)
      atomtype = trim(elsym(elnum))//'.'//intstr(hybnum)
   case (8, 16)
      hybnum = min(nadj + 1, 3)
      atomtype = trim(elsym(elnum))//'.'//intstr(hybnum)
   case default
      atomtype = elsym(elnum)
   end select

end function

end module
