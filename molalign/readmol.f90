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

module readmol
use parameters
use strutils
use chemdata

implicit none

contains

subroutine open2read(filepath, unit, fileext)
   character(*), intent(in) :: filepath
   integer, intent(out) :: unit
   character(:), allocatable, intent(out) :: fileext
   character(:), allocatable :: filename
   integer :: stat

   if (len(filepath) == 0) then
      write (error_unit, '(a)') 'Error: File path is empty'
      stop
   end if

   fileext = baseext(filepath)

   if (len(fileext) == 0) then
      write (error_unit, '(a,1x,a)') 'Missing file extension:', filepath
      stop
   end if

   open(newunit=unit, file=filepath, action='read', status='old', iostat=stat)
   if (stat /= 0) then
      write (error_unit, '(a,1x,a,1x,a)') 'Error opening', filepath, 'for reading'
      stop
   end if

end subroutine

subroutine readfile(unit, fmtin, natom, title, labels, coords, opt_nbond, opt_bonds)
   integer, intent(in) :: unit
   character(*), intent(in) :: fmtin
   integer, intent(out) :: natom
   character(:), allocatable, intent(out) :: title
   character(*), dimension(:), allocatable, intent(out) :: labels
   real(wp), dimension(:, :), allocatable, intent(out) :: coords
   integer, optional, intent(out) :: opt_nbond
   integer, target, optional, intent(out) :: opt_bonds(:, :)

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

   select case (fmtin)
   case ('xyz')
      call readxyzfile(unit, natom, title, labels, coords)
!    case ('mol2')
!        call readmol2file(unit, natom, title, labels, coords, nbond, bonds)
   case default
      write (error_unit, '(a,1x,a)') 'Invalid format:', fmtin
      stop
   end select

end subroutine

subroutine readxyzfile(unit, natom, title, labels, coords)
   integer, intent(in) :: unit
   integer, intent(out) :: natom
   real(wp), dimension(:, :), allocatable, intent(out) :: coords
   character(*), dimension(:), allocatable, intent(out) :: labels
   character(:), allocatable, intent(out) :: title
   character(maxstrlen) :: buffer
   integer :: i, stat

   read (unit, *, iostat=stat) natom
   if (stat < 0) then
      write (error_unit, '(a)') 'Unexpected end of file!'
      stop
   end if

   allocate (labels(natom), coords(3, natom))

   read (unit, '(a)', iostat=stat) buffer
   title = trim(buffer)
   if (stat < 0) then
      write (error_unit, '(a)') 'Unexpected end of file!'
      stop
   end if

   do i = 1, natom
      read (unit, *, iostat=stat) labels(i), coords(:, i)
      if (stat < 0) then
         write (error_unit, '(a)') 'Unexpected end of file!'
         stop
      end if
   end do

end subroutine

end module
