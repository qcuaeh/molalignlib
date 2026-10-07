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

module file_utils
use parameters
use str_utils
use chemdata
implicit none
private
public open2read
public open2write
public bondtype_code
public bondtype_str
public SDF_BONDTYPES

! Labels given by the readers to the MOL/SDF bond type numbers 1..10 (9
! and 10 are V3000 only), and written back by the writers: 1 single (1),
! 2 double (2), 3 triple (3), 4 aromatic (ar), 5 single or double (sd),
! 6 single or aromatic (sa), 7 double or aromatic (da), 9 coordination (co)
! and 10 hydrogen bond (hb). Number 8 (any) has the blank label, which
! stands for the undefined type (see bondtype_code), so it is read as
! UNDEFINED_BOND_TYPE.
character(2), dimension(*), parameter :: SDF_BONDTYPES = &
      [character(2) :: '1', '2', '3', 'ar', 'sd', 'sa', 'da', '', 'co', 'hb']

contains

subroutine split_path(path, dirname, filename)
   character(*), intent(in) :: path
   character(:), allocatable, intent(out) :: dirname, filename
   ! Local variables
   integer(ik) :: pos
   pos = index(path, '/', back=.TRUE.)
   if (pos /= 0) then
      dirname = path(:pos-1)
      filename = path(pos+1:)
      if (len(filename) == 0) then
         write (stderr, '(A,1X,A)') 'Invalid path'
         stop
      end if
   else
      dirname = ''
      filename = path
   end if
end subroutine

subroutine split_filename(filename, basename, extension)
   character(*), intent(in) :: filename
   character(:), allocatable, intent(out) :: basename, extension
   ! Local variables
   integer(ik) :: pos
   pos = index(filename, '.', back=.TRUE.)
   if (pos /= 0) then
      basename = filename(:pos-1)
      extension = filename(pos+1:)
      if (len(basename) == 0 .or. len(extension) == 0) then
         write (stderr, '(A,1X,A)') 'Invalid file name', filename
         stop
      end if
   else
      basename = ''
      extension = filename
   end if
end subroutine

subroutine open2read(path, extension, unit)
! Open path for reading and return its extension. A path without a base
! name (e.g. 'xyz') returns the stdin unit.
   character(*), intent(in) :: path
   character(:), allocatable, intent(out) :: extension
   integer(ik), intent(out) :: unit
   ! Local variables
   character(:), allocatable :: dirname, filename, basename
   integer(ik) :: stat

   call split_path(path, dirname, filename)
   call split_filename(filename, basename, extension)

   if (basename == '') then
      unit = stdin
      return
   end if

   open(newunit=unit, file=path, action='read', status='old', iostat=stat)
   if (stat /= 0) then
      write (stderr, '(A,1X,A,1X,A)') 'opening', path, 'for reading'
      stop
   end if
end subroutine

subroutine open2write(path, extension, unit)
! Open path for writing and return its extension. A bare extension (e.g.
! 'xyz') selects standard output.
   character(*), intent(in) :: path
   character(:), allocatable, intent(out) :: extension
   integer(ik), intent(out) :: unit
   ! Local variables
   character(:), allocatable :: dirname, filename, basename
   integer(ik) :: stat

   call split_path(path, dirname, filename)
   call split_filename(filename, basename, extension)

   if (basename == '') then
      unit = stdout
      return
   end if

   open(newunit=unit, file=path, action='write', status='replace', iostat=stat)
   if (stat /= 0) then
      write (stderr, '(A,1X,A,1X,A)') 'Can''t open', path, 'for writing'
      stop
   end if
end subroutine

function bondtype_code(typestr) result(code)
! Integer code of a bond type string, case insensitive, surrounding blanks
! ignored (see parameters):
!   d     a digit 1-9                         -> d
!   ab    a letter, then a letter or a digit  -> BOND_BLOCK2 + 36*a + b
!   nsd   a digit, a separator, a digit       -> BOND_BLOCK3 + 100*s + 10*n + d
! with a = 0..25 for a..z, b = 0..35 for 0..9, a..z, n, d = 0..9 and s the
! position (from 0) of the separator in BONDTYPE_SEPARATORS. The code
! identifies the string and nothing else: equal codes mean equal strings.
! The blank string stands for the undefined type and gives
! UNDEFINED_BOND_TYPE. -1 for a string of any other form (including 0),
! which is not a bond type.
   character(*), intent(in) :: typestr
   integer(ik) :: code
   ! Local variables
   character(:), allocatable :: trimmed
   integer(ik) :: char1, char2, separator

   code = -1
   trimmed = lowercase(trim(adjustl(typestr)))

   select case (len(trimmed))
   case (0)
      code = UNDEFINED_BOND_TYPE
   case (1)
      ! The digit's value; 0 is not a bond type
      char1 = index(BONDTYPE_DIGITS, trimmed) - 1
      if (char1 >= 1) code = char1
   case (2)
      char1 = index(BONDTYPE_LETTERS, trimmed(1:1)) - 1
      char2 = index(BONDTYPE_ALNUM, trimmed(2:2)) - 1
      if (char1 >= 0 .and. char2 >= 0) then
         code = BOND_BLOCK2 + len(BONDTYPE_ALNUM)*char1 + char2
      end if
   case (3)
      char1 = index(BONDTYPE_DIGITS, trimmed(1:1)) - 1
      char2 = index(BONDTYPE_DIGITS, trimmed(3:3)) - 1
      separator = index(BONDTYPE_SEPARATORS, trimmed(2:2)) - 1
      if (separator >= 0 .and. char1 >= 0 .and. char2 >= 0) then
         code = BOND_BLOCK3 + 100*separator + 10*char1 + char2
      end if
   end select
end function

function bondtype_str(code) result(typestr)
! Bond type string of a valid code, 1..MAX_BOND_TYPE, or of
! UNDEFINED_BOND_TYPE (the inverse of bondtype_code): one to three
! characters, lowercase, or the blank string for the undefined type
   integer(ik), intent(in) :: code
   character(:), allocatable :: typestr
   ! Local variables
   integer(ik) :: offset, i, j, k

   if (code >= 1 .and. code < BOND_BLOCK2) then
      typestr = BONDTYPE_DIGITS(code+1:code+1)
   else if (code >= BOND_BLOCK2 .and. code < BOND_BLOCK3) then
      offset = code - BOND_BLOCK2
      i = offset/len(BONDTYPE_ALNUM) + 1
      j = mod(offset, len(BONDTYPE_ALNUM)) + 1
      ! GNU Fortran 7.3.1 crashes if automatically allocated
      allocate (character(2)::typestr)
      typestr(:) = BONDTYPE_LETTERS(i:i)//BONDTYPE_ALNUM(j:j)
   else if (code >= BOND_BLOCK3 .and. code <= MAX_BOND_TYPE) then
      offset = code - BOND_BLOCK3
      k = offset/100 + 1
      i = mod(offset, 100)/10 + 1
      j = mod(offset, 10) + 1
      ! GNU Fortran 7.3.1 crashes if automatically allocated
      allocate (character(3)::typestr)
      typestr(:) = BONDTYPE_DIGITS(i:i)//BONDTYPE_SEPARATORS(k:k)//BONDTYPE_DIGITS(j:j)
   else if (code == UNDEFINED_BOND_TYPE) then
      typestr = ''
   else
      error stop 'bondtype_str: bond type code out of range'
   end if
end function

end module
