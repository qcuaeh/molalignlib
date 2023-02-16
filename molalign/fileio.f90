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

module fileio
use stdio
use readmol
use writemol

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

subroutine readfile(unit, fmtin, title, natom, labels, coords, adjmat)
   integer, intent(in) :: unit
   character(*), intent(in) :: fmtin
   integer, intent(out) :: natom
   character(:), allocatable, intent(out) :: title
   character(*), dimension(:), allocatable, intent(out) :: labels
   real(wp), dimension(:, :), allocatable, intent(out) :: coords
   logical, dimension(:, :), allocatable, intent(out) :: adjmat

   integer i
   integer, allocatable :: nbond
   integer, allocatable :: bonds(:, :)
   
   nbond = 0

   select case (fmtin)
   case ('xyz')
      call readxyz(unit, title, natom, labels, coords)
   case ('mol2')
      call readmol2(unit, title, natom, labels, coords, nbond, bonds)
   case default
      write (error_unit, '(a,1x,a)') 'Invalid format:', fmtin
      stop
   end select

   allocate(adjmat(natom, natom))

   adjmat(:, :) = .false.

   do i = 1, nbond
      adjmat(bonds(1, i), bonds(2, i)) = .true.
      adjmat(bonds(2, i), bonds(1, i)) = .true.
   end do

end subroutine

subroutine writefile(unit, fmtout, title, natom, znums, coords, adjmat)
   integer, intent(in) :: unit, natom
   integer, dimension(:), intent(in) :: znums
   real(wp), dimension(:, :), intent(in) :: coords
   logical, dimension(:, :), intent(in) :: adjmat
   character(*), intent(in) :: title, fmtout
   integer :: i, j, nbond, bonds(2, natom*maxcoord)

   nbond = 0
   do i = 1, natom
      do j = i + 1, natom
         if (adjmat(i, j)) then
            nbond = nbond + 1
            bonds(1, nbond) = i
            bonds(2, nbond) = j
         end if
      end do
   end do

   select case (fmtout)
   case ('xyz')
      call writexyz(unit, title, natom, znums, coords)
   case ('mol2')
      call writemol2(unit, title, natom, znums, coords, nbond, bonds)
   case default
      write (error_unit, '(a,1x,a)') 'Invalid format:', fmtout
      stop
   end select

end subroutine

end module
