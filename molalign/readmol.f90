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
use stdio
use bounds
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

subroutine readfile(unit, fmtin, title, natom, labels, coords)
   integer, intent(in) :: unit
   character(*), intent(in) :: fmtin
   integer, intent(out) :: natom
   character(:), allocatable, intent(out) :: title
   character(*), dimension(:), allocatable, intent(out) :: labels
   real(wp), dimension(:, :), allocatable, intent(out) :: coords

   integer i
   integer, allocatable :: nbond
   integer, allocatable :: bonds(:, :)
   
   nbond = 0

   select case (fmtin)
   case ('xyz')
      call readxyzfile(unit, title, natom, labels, coords)
   case ('mol2')
      call readmol2file(unit, title, natom, labels, coords, nbond, bonds)
   case default
      write (error_unit, '(a,1x,a)') 'Invalid format:', fmtin
      stop
   end select

end subroutine

subroutine readxyzfile(unit, title, natom, labels, coords)
   integer, intent(in) :: unit
   integer, intent(out) :: natom
   real(wp), dimension(:, :), allocatable, intent(out) :: coords
   character(*), dimension(:), allocatable, intent(out) :: labels
   character(:), allocatable, intent(out) :: title
   character(ll) :: buffer

   integer :: i, stat

   read (unit, *, end=99) natom

   allocate(labels(natom), coords(3, natom))

   read (unit, '(a)', end=99) buffer
   title = trim(buffer)

   do i = 1, natom
      read (unit, *, end=99) labels(i), coords(:, i)
   end do

   return

   99 continue
   write (error_unit, '(a)') 'Unexpected end of file!'
   stop

end subroutine

subroutine readmol2file(unit, title, natom, labels, coords, nbond, bonds)
   integer, intent(in) :: unit
   integer, intent(out) :: natom, nbond
   integer, dimension(:, :), allocatable, intent(out) :: bonds
   real(wp), dimension(:, :), allocatable, intent(out) :: coords
   character(*), dimension(:), allocatable, intent(out) :: labels
   character(:), allocatable, intent(out) :: title
   character(ll) :: buffer
   integer :: i, id

   do
      read (unit, '(a)', end=99) buffer
      if (buffer == '@<TRIPOS>MOLECULE') exit
   end do

   read (unit, '(a)', end=99) buffer
   title = trim(buffer)
   read (unit, *, end=99) natom, nbond

   allocate(labels(natom), coords(3, natom))
   allocate(bonds(2, nbond))

   do
      read (unit, '(a)', end=99) buffer
      if (buffer == '@<TRIPOS>ATOM') exit
   end do

   do i = 1, natom
      read (unit, *, end=99) id, labels(i), coords(:,i)
   end do

   do
      read (unit, '(a)', end=99) buffer
      if (buffer == '@<TRIPOS>BOND') exit
   end do

   do i = 1, nbond
      read (unit, *, end=99) id, bonds(:,i)
   end do

   return

   99 continue
   write (error_unit, '(a)') 'Unexpected end of file!'
   stop

end subroutine

end module
