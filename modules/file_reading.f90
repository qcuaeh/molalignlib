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

module file_reading
use parameters
use permutation
use chemistry
use molecule
use options
implicit none

contains

subroutine open2read(filepath, filetype, fileunit)
   character(*), intent(in) :: filepath
   character(:), allocatable, intent(out) :: filetype
   integer, intent(out) :: fileunit
   integer :: stat

   call parse_path(filepath, filetype)
   open(newunit=fileunit, file=filepath, action='read', status='old', iostat=stat)
   if (stat /= 0) then
      write (stderr, '(A,1X,A,1X,A)') 'Error: opening', filepath, 'for reading'
      stop 1
   end if
end subroutine

subroutine read_file(unit, typein, title, atoms, bonds)
   integer, intent(in) :: unit
   character(*), intent(in) :: typein
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds

   select case (typein)
   case ('xyz')
      call readfile_xyz(unit, title, atoms, bonds)
   case ('mol')
      call readfile_mol(unit, title, atoms, bonds)
   case ('sdf')
      call readfile_sdf(unit, title, atoms, bonds)
   case ('mol2')
      call readfile_mol2(unit, title, atoms, bonds)
   case default
      write (stderr, '(A,A,A)') 'Error: File format "', typein, '" is not supported'
      stop 1
   end select
end subroutine

subroutine readfile_xyz(unit, title, atoms, bonds)
   integer, intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   character(ll) :: buffer
   character(wl) :: elsym
   character(:), allocatable :: label
   integer :: i, atoms_size, elnum, stat
   real(rk) :: coords(3)

   ! Read number of atoms
   read (unit, *, iostat=stat) atoms_size
   if (stat /= 0) then
      write (stderr, '(A)') 'Error: Invalid XYZ format'
      stop 1
   end if

   ! Check for empty file
   if (atoms_size <= 0) then
      write (stderr, '(A)') 'Error: File contains no atoms'
      stop 1
   end if

   allocate (atoms(atoms_size))
   allocate (bonds(0))

   ! Read title line
   read (unit,'(A)',iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Error: Invalid XYZ format'
      stop 1
   end if
   title = trim(buffer)

   do i = 1, atoms_size
      read (unit,*,iostat=stat) elsym, coords
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid XYZ format'
         stop 1
      end if

      call split_symbol(elsym, elnum, label)
      atoms(i)%elnum = elnum
      atoms(i)%coords = coords
      if (label_flag .and. label /= '') then
         read (label,*,iostat=stat) atoms(i)%typeid
         if (stat /= 0) then
            write (stderr, '(A,A)') 'Error: Invalid label', label
            stop 1
         end if
      end if
   end do
end subroutine

subroutine readfile_mol(unit, title, atoms, bonds)
   integer, intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, target, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   character(ll) :: buffer
   logical :: is_v3000
   integer :: stat

   ! Read header block (3 lines)
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Error: Invalid MOL format'
      stop 1
   end if
   title = trim(buffer)

   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Error: Invalid MOL format'
      stop 1
   end if

   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Error: Invalid MOL format'
      stop 1
   end if

   ! Read counts line and determine format version
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Error: Invalid MOL format'
      stop 1
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
      write (stderr, '(A)') 'Error: File contains no atoms'
      stop 1
   end if
end subroutine

subroutine readfile_sdf(unit, title, atoms, bonds)
   integer, intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   character(ll) :: buffer
   integer :: stat

   ! Read the molecule using MOL format reader
   call readfile_mol(unit, title, atoms, bonds)

   ! Continue reading until we reach $$$$ or end of file
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat < 0) then
         ! End of file reached, this is acceptable here
         exit
      else if (stat > 0) then
         write (stderr, '(A)') 'Error: Reading SDF property data'
         stop 1
      end if
      if (trim(buffer) == '$$$$') then
         ! End of molecule reached
         exit
      end if
   end do
end subroutine

subroutine read_v2000_format(unit, counts_line, atoms, bonds)
   integer, intent(in) :: unit
   character(*), intent(in) :: counts_line
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   real(rk) :: coords(3)
   character(3) :: elsym
   character(ll) :: buffer
   character(:), allocatable :: label
   integer :: elnum, atomidx1, atomidx2, typeid
   integer :: atoms_size, bonds_size, stat, i

   ! Parse counts from V2000 format
   read (counts_line(1:3), '(I3)', iostat=stat) atoms_size
   if (stat /= 0) then
      write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
      stop 1
   end if
   read (counts_line(4:6), '(I3)', iostat=stat) bonds_size
   if (stat /= 0) then
      write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
      stop 1
   end if

   allocate (atoms(atoms_size))
   ! Read atom block
   do i = 1, atoms_size
      read (unit,'(A)',iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if
      read (buffer(1:10), '(F10.4)', iostat=stat) coords(1)
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if
      read (buffer(11:20), '(F10.4)', iostat=stat) coords(2)
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if
      read (buffer(21:30), '(F10.4)', iostat=stat) coords(3)
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if

      elsym = adjustl(buffer(32:34))
      call split_symbol(elsym, elnum, label)
      atoms(i)%elnum = elnum
      atoms(i)%coords = coords
   end do

   ! Read bond block
   allocate (bonds(bonds_size))
   do i = 1, bonds_size
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if
      read (buffer(1:3), '(I3)', iostat=stat) atomidx1
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if
      read (buffer(4:6), '(I3)', iostat=stat) atomidx2
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if
      read (buffer(7:9), '(I3)', iostat=stat) typeid
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if

      bonds(i)%atomidx1 = atomidx1
      bonds(i)%atomidx2 = atomidx2
      bonds(i)%typeid = typeid
   end do
end subroutine

subroutine read_v3000_format(unit, counts_line, atoms, bonds)
   integer, intent(in) :: unit
   character(*), intent(in) :: counts_line
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   character(ll) :: buffer
   character(:), allocatable :: label
   real(rk) :: coords(3)
   integer :: atoms_size, bonds_size, stat, i, pos
   integer :: elnum, typeid, atomidx1, atomidx2, dummy_int
   character(3) :: elsym

   ! Find BEGIN CTAB
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if
      if (index(buffer, 'BEGIN CTAB') > 0) exit
   end do

   ! Find COUNTS line in V3000 format
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if
      if (index(buffer, 'COUNTS') > 0) then
         ! Parse: M  V30 COUNTS atoms_size bonds_size ...
         pos = index(buffer, 'COUNTS')
         read (buffer(pos+6:), *, iostat=stat) atoms_size, bonds_size
         if (stat /= 0) then
            write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
            stop 1
         end if
         exit
      end if
   end do

   allocate (atoms(atoms_size))

   ! Find BEGIN ATOM
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if
      if (index(buffer, 'BEGIN ATOM') > 0) exit
   end do

   ! Read atoms in V3000 format
   ! Format: M  V30 atom_id element_symbol x y z charge ...
   do i = 1, atoms_size
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if
      pos = index(buffer, 'V30')
      read (buffer(pos+3:), *, iostat=stat) dummy_int, elsym, coords
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if

      call split_symbol(elsym, elnum, label)
      atoms(i)%elnum = elnum
      atoms(i)%coords = coords
      if (label_flag .and. label /= '') then
         read (label,*,iostat=stat) atoms(i)%typeid
         if (stat /= 0) then
            write (stderr, '(A,A)') 'Error: Invalid label', label
            stop 1
         end if
      end if
   end do

   ! Find END ATOM
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if
      if (index(buffer, 'END ATOM') > 0) exit
   end do

   allocate (bonds(bonds_size))
   if (bonds_size > 0) then
      ! Find BEGIN BOND
      do
         read (unit, '(A)', iostat=stat) buffer
         if (stat /= 0) then
            write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
            stop 1
         end if
         if (index(buffer, 'BEGIN BOND') > 0) exit
      end do

      ! Read bonds in V3000 format
      ! Format: M  V30 bond_id type atom1 atom2 ...
      do i = 1, bonds_size
         read (unit, '(A)', iostat=stat) buffer
         if (stat /= 0) then
            write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
            stop 1
         end if
         pos = index(buffer, 'V30')
         read (buffer(pos+3:), *, iostat=stat) dummy_int, typeid, atomidx1, atomidx2
         if (stat /= 0) then
            write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
            stop 1
         end if

         bonds(i)%atomidx1 = atomidx1
         bonds(i)%atomidx2 = atomidx2
         bonds(i)%typeid = typeid
      end do

      ! Find END BOND
      do
         read (unit, '(A)', iostat=stat) buffer
         if (stat /= 0) then
            write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
            stop 1
         end if
         if (index(buffer, 'END BOND') > 0) exit
      end do
   end if

   ! Find END CTAB
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL/SDF format'
         stop 1
      end if
      if (index(buffer, 'END CTAB') > 0) exit
   end do
end subroutine

subroutine readfile_mol2(unit, title, atoms, bonds)
   integer, intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, target, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   real(rk) :: coords(3)
   character(wl) :: elsym, dummy, typestr
   character(ll) :: buffer
   character(:), allocatable :: label
   integer :: atoms_size, bonds_size, elnum, stat
   integer :: i, atomidx1, atomidx2

   ! Find @<TRIPOS>MOLECULE section
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL2 format'
         stop 1
      end if
      if (trim(buffer) == '@<TRIPOS>MOLECULE') exit
   end do

   ! Read molecule name
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Error: Invalid MOL2 format'
      stop 1
   end if
   title = trim(buffer)

   ! Read counts line
   read (unit, *, iostat=stat) atoms_size, bonds_size
   if (stat /= 0) then
      write (stderr, '(A)') 'Error: Invalid MOL2 format'
      stop 1
   end if

   ! Check for empty file
   if (atoms_size <= 0) then
      write (stderr, '(A)') 'Error: File contains no atoms'
      stop 1
   end if

   allocate (atoms(atoms_size))

   ! Find @<TRIPOS>ATOM section
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL2 format'
         stop 1
      end if
      if (trim(buffer) == '@<TRIPOS>ATOM') exit
   end do

   ! Read atom section
   ! Format: atom_id atom_name x y z typestr [subst_id subst_name charge]
   do i = 1, atoms_size
      read (unit,*,iostat=stat) dummy, elsym, coords, typestr
      if (stat /= 0) then
         write (stderr, '(A)') 'Error: Invalid MOL2 format'
         stop 1
      end if

      call split_symbol(elsym, elnum, label)
      atoms(i)%elnum = elnum
      atoms(i)%coords = coords
      if (label_flag .and. label /= '') then
         read (label,*,iostat=stat) atoms(i)%typeid
         if (stat /= 0) then
            write (stderr, '(A,A)') 'Error: Invalid label', label
            stop 1
         end if
      end if
   end do

   ! Read bonds
   allocate (bonds(bonds_size))
   if (bonds_size > 0) then
      ! Find @<TRIPOS>BOND section
      do
         read (unit, '(A)', iostat=stat) buffer
         if (stat /= 0) then
            write (stderr, '(A)') 'Error: Invalid MOL2 format'
            stop 1
         end if
         if (trim(buffer) == '@<TRIPOS>BOND') exit
      end do

      ! Read bond section
      ! Format: bond_id origin_atom_id target_atom_id typestr
      do i = 1, bonds_size
         read (unit, *, iostat=stat) dummy, atomidx1, atomidx2, typestr
         if (stat /= 0) then
            write (stderr, '(A)') 'Error: Invalid MOL2 format'
            stop 1
         end if

         bonds(i)%atomidx1 = atomidx1
         bonds(i)%atomidx2 = atomidx2
!         bonds(i)%typeid = type_index(typestr)
      end do
   end if
end subroutine

end module
