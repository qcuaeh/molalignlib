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
use io
use bounds
use strutils
use chemdata

implicit none

contains

subroutine open2write(filename, unit)
   character(*), intent(in) :: filename
   integer, intent(out) :: unit
   integer :: stat

   open(newunit=unit, file=filename, action='write', status='replace', iostat=stat)
   if (stat /= 0) then
      write (error_unit, '(a,1x,a,1x,a)') 'Error opening', filename, 'for writing'
      stop
   end if

end subroutine


subroutine writefile(unit, fmtout, natom, title, znums, coords, opt_nbond, opt_bonds)
   integer, intent(in) :: unit, natom
   integer, dimension(:), intent(in) :: znums
   integer, optional, intent(in) :: opt_nbond
   integer, target, optional, intent(in) :: opt_bonds(:, :)
   real(wp), dimension(:, :), intent(in) :: coords
   character(*), intent(in) :: title, fmtout

   integer :: nbond
   integer, pointer :: bonds(:, :)
   
   if (present(opt_nbond)) then
      nbond = opt_nbond
      if (present(opt_bonds)) then
         bonds => opt_bonds
      else
         write (error_unit, '(a)') 'Bond list is missing!'
         stop
      end if
   else
      nbond = 0
   end if

   select case (fmtout)
   case ('xyz')
      call writexyzfile(unit, natom, title, znums, coords)
   case ('mol2')
      call writemol2file(unit, natom, title, znums, coords, nbond, bonds)
   case default
      write (error_unit, '(a,1x,a)') 'Invalid format:', fmtout
      stop
   end select

end subroutine

subroutine writexyzfile(unit, natom, title, znums, coords)
   character(*), intent(in) :: title
   integer, intent(in) :: unit, natom
   integer, dimension(:), intent(in) :: znums
   real(wp), dimension(:, :), intent(in) :: coords
   integer :: i

   write (unit, '(i0)') natom
   write (unit, '(a)') title

   do i = 1, natom
      write (unit, '(a,3(2x,f12.6))') elsym(znums(i)), coords(:, i)
   end do

end subroutine

subroutine writemol2file(unit, natom, title, znums, coords, nbond, bonds)
   character(*), intent(in) :: title
   integer, intent(in) :: unit, natom
   integer, dimension(:), intent(in) :: znums
   real(wp), dimension(:, :), intent(in) :: coords
   integer, intent(in) :: nbond
   integer, target, intent(in) :: bonds(:, :)
   integer :: i
   
   write (unit, '(a)') '@<TRIPOS>MOLECULE'
   if (title /= '') then
      write (unit, '(a)') title
   else
      write (unit, '(a)') 'Untitled'
   end if
   write (unit, '(5(i4,1x))') natom, nbond, 0, 0, 0
   write (unit, '(a)') 'SMALL'
   write (unit, '(a)') 'NO_CHARGES'
   write (unit, '(a)') '@<TRIPOS>ATOM'

   do i = 1, natom
      write (unit, '(i4,1x,a2,3(1x,f12.6),1x,a8,1x,i2,1x,a4,1x,f7.3)') &
         i, elsym(znums(i)), coords(:, i), elsym(znums(i)), 1, 'MOL1', 0.
   end do

   write (unit, '(a)') '@<TRIPOS>BOND'

   do i = 1, nbond
      write (unit, '(i4,1x,2(1x,i4),1x,a2)') i, bonds(:, i), '1'
   end do

end subroutine

end module
