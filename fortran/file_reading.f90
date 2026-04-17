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

module file_reading
use parameters
use str_utils
use chemdata
use molecule
use flags

implicit none

private
public read_file

contains

subroutine parse_label( label, elnum, group)
   character(*), intent(in) :: label
   integer(ik), intent(out) :: elnum, group
   ! Local variables
   character(:), allocatable :: normalized_label, elsym
   character(symlen), dimension(:), allocatable :: normalized_atomic_symbols
   integer(ik) :: pos, z

   normalized_label = lowercase(trim(adjustl(label)))
   normalized_atomic_symbols = lowercase(atomic_symbols)
   pos = verify(normalized_label, LOWERCHAR)

   if (pos == 0) then
      elsym = normalized_label
      group = 0
   else
      if (verify(normalized_label(pos:), NUMCHAR) == 0) then
         elsym = normalized_label(:pos-1)
         group = int(normalized_label(pos:))
      else
         write (stderr, '(A,1X,A)') 'Invalid atomic label:', normalized_label
         stop
      end if
   end if

   elnum = 0
   do z = 1, num_elems
      if (elsym == normalized_atomic_symbols(z)) then
         elnum = z
         return
      end if
   end do
end subroutine

subroutine read_file(unit, in_format, title, atoms, bonds)
   integer(ik), intent(in) :: unit
   character(*), intent(in) :: in_format
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds

   select case (in_format)
   case ('xyz')
      call read_file_xyz(unit, title, atoms, bonds)
   case ('mol')
      call read_file_mol(unit, title, atoms, bonds)
   case ('sdf')
      call read_file_sdf(unit, title, atoms, bonds)
   case ('mol2')
      call read_file_mol2(unit, title, atoms, bonds)
   case default
      write (stderr, '(A,A,A)') 'File format "', in_format, '" is not supported'
      stop
   end select
end subroutine

subroutine read_file_xyz(unit, title, atoms, bonds)
   integer(ik), intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   character(ll) :: label
   character(ll) :: buffer
   integer(ik) :: elnum, group
   integer(ik) :: i, n_atoms, stat
   real(rk) :: coords(3)

   ! Read number of atoms
   read (unit, *, iostat=stat) n_atoms
   if (stat /= 0) then
      stop 'Invalid XYZ format'
   end if

   ! Check for empty file
   if (n_atoms <= 0) then
      stop 'File contains no atoms'
   end if

   allocate (atoms(n_atoms))
   allocate (bonds(0))

   ! Read title line
   read (unit,'(A)',iostat=stat) buffer
   if (stat /= 0) then
      stop 'Invalid XYZ format'
   end if
   title = trim(buffer)

   do i = 1, n_atoms
      read (unit,*,iostat=stat) label, coords
      if (stat /= 0) then
         stop 'Invalid XYZ format'
      end if

      call parse_label(label, elnum, group)
      atoms(i)%elnum = elnum
      atoms(i)%coords = coords
      if (label_flag) then
         atoms(i)%group = group
      end if
   end do
end subroutine

subroutine read_file_mol(unit, title, atoms, bonds)
   integer(ik), intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, target, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   character(ll) :: buffer
   logical(lk) :: is_v3000
   integer(ik) :: stat

   ! Read header block (3 lines)
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      stop 'Invalid MOL format'
   end if
   title = trim(buffer)

   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      stop 'Invalid MOL format'
   end if

   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      stop 'Invalid MOL format'
   end if

   ! Read counts line and determine format version
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      stop 'Invalid MOL format'
   end if

   ! Check if this is V3000 format
   is_v3000 = index(buffer, 'V3000') > 0

   if (is_v3000) then
      call read_v3000_format(unit, buffer, atoms, bonds)
   else
      call read_v2000_format(unit, buffer, atoms, bonds)
   end if

   ! Check for empty molecule
   if (size(atoms) <= 0) then
      stop 'File contains no atoms'
   end if
end subroutine

subroutine read_file_sdf(unit, title, atoms, bonds)
   integer(ik), intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   character(ll) :: buffer
   integer(ik) :: stat

   ! Read the molecule using MOL format reader
   call read_file_mol(unit, title, atoms, bonds)

   ! Continue reading until we reach $$$$ or end of file
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat < 0) then
         ! End of file reached, this is acceptable here
         exit
      else if (stat > 0) then
         stop 'Reading SDF property data'
      end if
      if (trim(buffer) == '$$$$') then
         ! End of molecule reached
         exit
      end if
   end do
end subroutine

subroutine read_v2000_format(unit, counts_line, atoms, bonds)
   integer(ik), intent(in) :: unit
   character(*), intent(in) :: counts_line
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   real(rk) :: coords(3)
   character(ll) :: label
   character(ll) :: buffer
   integer(ik) :: elnum, group
   integer(ik) :: atomidx1, atomidx2, bondtype
   integer(ik) :: n_atoms, n_bonds, stat, i

   ! Parse counts from V2000 format
   read (counts_line(1:3), '(I3)', iostat=stat) n_atoms
   if (stat /= 0) then
      stop 'Invalid MOL/SDF format'
   end if
   read (counts_line(4:6), '(I3)', iostat=stat) n_bonds
   if (stat /= 0) then
      stop 'Invalid MOL/SDF format'
   end if

   allocate (atoms(n_atoms))
   ! Read atom block
   do i = 1, n_atoms
      read (unit,'(A)',iostat=stat) buffer
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if
      read (buffer(1:10), '(F10.4)', iostat=stat) coords(1)
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if
      read (buffer(11:20), '(F10.4)', iostat=stat) coords(2)
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if
      read (buffer(21:30), '(F10.4)', iostat=stat) coords(3)
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if

      label = adjustl(buffer(32:34))
      call parse_label(label, elnum, group)
      atoms(i)%elnum = elnum
      atoms(i)%coords = coords
      if (label_flag) then
         atoms(i)%group = group
      end if
   end do

   ! Read bond block
   allocate (bonds(n_bonds))
   do i = 1, n_bonds
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if
      read (buffer(1:3), '(I3)', iostat=stat) atomidx1
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if
      read (buffer(4:6), '(I3)', iostat=stat) atomidx2
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if
      read (buffer(7:9), '(I3)', iostat=stat) bondtype
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if

      bonds(i)%atomidx1 = atomidx1
      bonds(i)%atomidx2 = atomidx2
      bonds(i)%bondtype = bondtype
   end do
end subroutine

subroutine read_v3000_format(unit, counts_line, atoms, bonds)
   integer(ik), intent(in) :: unit
   character(*), intent(in) :: counts_line
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   character(ll) :: buffer
   integer(ik) :: elnum, group
   real(rk) :: coords(3)
   integer(ik) :: n_atoms, n_bonds, stat, i, pos
   integer(ik) :: bondtype, atomidx1, atomidx2, dummy_int
   character(3) :: label

   ! Find BEGIN CTAB
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if
      if (index(buffer, 'BEGIN CTAB') > 0) exit
   end do

   ! Find COUNTS line in V3000 format
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if
      if (index(buffer, 'COUNTS') > 0) then
         ! Parse: M  V30 COUNTS n_atoms n_bonds ...
         pos = index(buffer, 'COUNTS')
         read (buffer(pos+6:), *, iostat=stat) n_atoms, n_bonds
         if (stat /= 0) then
            stop 'Invalid MOL/SDF format'
         end if
         exit
      end if
   end do

   allocate (atoms(n_atoms))

   ! Find BEGIN ATOM
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if
      if (index(buffer, 'BEGIN ATOM') > 0) exit
   end do

   ! Read atoms in V3000 format
   ! Format: M  V30 atom_id element_symbol x y z charge ...
   do i = 1, n_atoms
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if
      pos = index(buffer, 'V30')
      read (buffer(pos+3:), *, iostat=stat) dummy_int, label, coords
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if

      call parse_label(label, elnum, group)
      atoms(i)%elnum = elnum
      atoms(i)%coords = coords
      if (label_flag) then
         atoms(i)%group = group
      end if
   end do

   ! Find END ATOM
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if
      if (index(buffer, 'END ATOM') > 0) exit
   end do

   allocate (bonds(n_bonds))
   if (n_bonds > 0) then
      ! Find BEGIN BOND
      do
         read (unit, '(A)', iostat=stat) buffer
         if (stat /= 0) then
            stop 'Invalid MOL/SDF format'
         end if
         if (index(buffer, 'BEGIN BOND') > 0) exit
      end do

      ! Read bonds in V3000 format
      ! Format: M  V30 bond_id type atom1 atom2 ...
      do i = 1, n_bonds
         read (unit, '(A)', iostat=stat) buffer
         if (stat /= 0) then
            stop 'Invalid MOL/SDF format'
         end if
         pos = index(buffer, 'V30')
         read (buffer(pos+3:), *, iostat=stat) dummy_int, bondtype, atomidx1, atomidx2
         if (stat /= 0) then
            stop 'Invalid MOL/SDF format'
         end if

         bonds(i)%atomidx1 = atomidx1
         bonds(i)%atomidx2 = atomidx2
         bonds(i)%bondtype = bondtype
      end do

      ! Find END BOND
      do
         read (unit, '(A)', iostat=stat) buffer
         if (stat /= 0) then
            stop 'Invalid MOL/SDF format'
         end if
         if (index(buffer, 'END BOND') > 0) exit
      end do
   end if

   ! Find END CTAB
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         stop 'Invalid MOL/SDF format'
      end if
      if (index(buffer, 'END CTAB') > 0) exit
   end do
end subroutine

subroutine read_file_mol2(unit, title, atoms, bonds)
   integer(ik), intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, target, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   real(rk) :: coords(3)
   character(ll) :: label, dummy, typestr
   character(ll) :: buffer
   integer(ik) :: elnum, group
   integer(ik) :: n_atoms, n_bonds, stat
   integer(ik) :: i, atomidx1, atomidx2

   ! Find @<TRIPOS>MOLECULE section
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         stop 'Invalid MOL2 format'
      end if
      if (trim(buffer) == '@<TRIPOS>MOLECULE') exit
   end do

   ! Read molecule name
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      stop 'Invalid MOL2 format'
   end if
   title = trim(buffer)

   ! Read counts line
   read (unit, *, iostat=stat) n_atoms, n_bonds
   if (stat /= 0) then
      stop 'Invalid MOL2 format'
   end if

   ! Check for empty file
   if (n_atoms <= 0) then
      stop 'File contains no atoms'
   end if

   allocate (atoms(n_atoms))

   ! Find @<TRIPOS>ATOM section
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         stop 'Invalid MOL2 format'
      end if
      if (trim(buffer) == '@<TRIPOS>ATOM') exit
   end do

   ! Read atom section
   ! Format: atom_id atom_name x y z typestr [subst_id subst_name charge]
   do i = 1, n_atoms
      read (unit,*,iostat=stat) dummy, label, coords, typestr
      if (stat /= 0) then
         stop 'Invalid MOL2 format'
      end if

      call parse_label(label, elnum, group)
      atoms(i)%elnum = elnum
      atoms(i)%coords = coords
      if (label_flag) then
         atoms(i)%group = group
      end if
   end do

   ! Read bonds
   allocate (bonds(n_bonds))
   if (n_bonds > 0) then
      ! Find @<TRIPOS>BOND section
      do
         read (unit, '(A)', iostat=stat) buffer
         if (stat /= 0) then
            stop 'Invalid MOL2 format'
         end if
         if (trim(buffer) == '@<TRIPOS>BOND') exit
      end do

      ! Read bond section
      ! Format: bond_id origin_atom_id target_atom_id typestr
      do i = 1, n_bonds
         read (unit, *, iostat=stat) dummy, atomidx1, atomidx2, typestr
         if (stat /= 0) then
            stop 'Invalid MOL2 format'
         end if

         bonds(i)%atomidx1 = atomidx1
         bonds(i)%atomidx2 = atomidx2
!         bonds(i)%bondtype = type_index(typestr)
      end do
   end if
end subroutine

end module
