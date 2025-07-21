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

module file_read
use parameters
use options
use chemistry
use molecule
use permutation
implicit none

contains

subroutine open2read(filepath, unit)
   character(*), intent(in) :: filepath
   integer, intent(out) :: unit
   integer :: stat

   open(newunit=unit, file=filepath, action='read', status='old', iostat=stat)
   if (stat /= 0) then
      write (stderr, '(A,1X,A,1X,A)') 'Error opening', filepath, 'for reading'
      stop
   end if
end subroutine

subroutine readfile(unit, exten, title, atoms, bonds)
   integer, intent(in) :: unit
   character(*), intent(in) :: exten
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds

   select case (exten)
   case ('xyz')
      call readfile_xyz(unit, title, atoms, bonds)
   case ('mol')
      call readfile_mol(unit, title, atoms, bonds)
   case ('sdf')
      call readfile_sdf(unit, title, atoms, bonds)
   case ('mol2')
      call readfile_mol2(unit, title, atoms, bonds)
   case default
      write (stderr, '(A,1X,A)') 'Invalid file extension:', exten
      stop
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
   integer :: i, num_atoms, elnum, typeidx, stat
   real(rk) :: coords(3)

   ! Read number of atoms
   read (unit, *, iostat=stat) num_atoms
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid XYZ file'
      stop
   end if

   allocate (atoms(num_atoms))
   allocate (bonds(0))

   ! Read title line
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid XYZ file'
      stop
   end if
   title = trim(buffer)

   do i = 1, num_atoms
      read (unit, *, iostat=stat) elsym, coords
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid XYZ file'
         stop
      end if

      call split_symbol(elsym, elnum, typeidx)
      atoms(i)%elnum = elnum
      atoms(i)%typeidx = typeidx
      atoms(i)%coords = coords
   end do
end subroutine

subroutine readfile_mol(unit, title, atoms, bonds)
   integer, intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, target, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   real(rk) :: coords(3)
   character(3) :: elsym
   character(ll) :: buffer
   integer :: elnum, atomidx1, atomidx2, typeidx
   integer :: num_atoms, num_bonds, stat, i

   ! Read header block (3 lines)
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL file'
      stop
   end if
   title = trim(buffer)

   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL file'
      stop
   end if

   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL file'
      stop
   end if

   ! Read counts line
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL file'
      stop
   end if

   read (buffer(1:3), '(I3)', iostat=stat) num_atoms
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL file'
      stop
   end if
   read (buffer(4:6), '(I3)', iostat=stat) num_bonds
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL file'
      stop
   end if

   allocate (atoms(num_atoms))
   ! Read atom block
   do i = 1, num_atoms
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL file'
         stop
      end if
      read (buffer(1:10), '(F10.4)', iostat=stat) coords(1)
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL file'
         stop
      end if
      read (buffer(11:20), '(F10.4)', iostat=stat) coords(2)
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL file'
         stop
      end if
      read (buffer(21:30), '(F10.4)', iostat=stat) coords(3)
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL file'
         stop
      end if

      elsym = adjustl(buffer(32:34))
      call split_symbol(elsym, elnum, typeidx)
      atoms(i)%elnum = elnum
      atoms(i)%typeidx = typeidx
      atoms(i)%coords = coords
   end do

   ! Read bond block
   allocate (bonds(num_bonds))
   do i = 1, num_bonds
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL file'
         stop
      end if
      read (buffer(1:3), '(I3)', iostat=stat) atomidx1
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL file'
         stop
      end if
      read (buffer(4:6), '(I3)', iostat=stat) atomidx2
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL file'
         stop
      end if
      read (buffer(7:9), '(I3)', iostat=stat) typeidx
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL file'
         stop
      end if

      bonds(i)%atomidx1 = atomidx1
      bonds(i)%atomidx2 = atomidx2
      bonds(i)%typeidx = typeidx
   end do
end subroutine

subroutine readfile_sdf(unit, title, atoms, bonds)
   integer, intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   character(ll) :: buffer
   logical :: end_of_molecule
   integer :: stat

   ! Read the MOL block (SDF files contain MOL blocks followed by data)
   call readfile_mol(unit, title, atoms, bonds)

   ! Skip any property data until we reach $$ or end of file
   end_of_molecule = .false.
   do while (.not. end_of_molecule)
      read (unit, '(A)', iostat=stat) buffer
      if (stat < 0) then
         ! End of file reached - this is acceptable for SDF files
         exit
      else if (stat > 0) then
         write (stderr, '(A)') 'Error reading SDF property data!'
         stop
      end if

      if (trim(buffer) == '$$') then
         end_of_molecule = .true.
      end if
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
   integer :: num_atoms, num_bonds, elnum, typeidx, stat
   character(ll) :: buffer
   integer :: i, atomidx1, atomidx2

   ! Find @<TRIPOS>MOLECULE section
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL2 file'
         stop
      end if
      if (trim(buffer) == '@<TRIPOS>MOLECULE') exit
   end do

   ! Read molecule name
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL2 file'
      stop
   end if
   title = trim(buffer)

   ! Read counts line
   read (unit, *, iostat=stat) num_atoms, num_bonds
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL2 file'
      stop
   end if

   allocate (atoms(num_atoms))

   ! Find @<TRIPOS>ATOM section
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL2 file'
         stop
      end if
      if (trim(buffer) == '@<TRIPOS>ATOM') exit
   end do

   ! Read atom section
   ! Format: atom_id atom_name x y z typestr [subst_id subst_name charge]
   do i = 1, num_atoms
      read (unit, *, iostat=stat) dummy, elsym, coords, typestr
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL2 file'
         stop
      end if

      call split_symbol(elsym, elnum, typeidx)
      atoms(i)%elnum = elnum
      atoms(i)%typeidx = typeidx
      atoms(i)%coords = coords
      ! typestr is now available in local variable for future use
   end do

   ! Read bonds
   allocate (bonds(num_bonds))
   if (num_bonds > 0) then
      ! Find @<TRIPOS>BOND section
      do
         read (unit, '(A)', iostat=stat) buffer
         if (stat /= 0) then
            write (stderr, '(A)') 'Invalid MOL2 file'
            stop
         end if
         if (trim(buffer) == '@<TRIPOS>BOND') exit
      end do

      ! Read bond section
      ! Format: bond_id origin_atom_id target_atom_id typestr
      do i = 1, num_bonds
         read (unit, *, iostat=stat) dummy, atomidx1, atomidx2, typestr
         if (stat /= 0) then
            write (stderr, '(A)') 'Invalid MOL2 file'
            stop
         end if

         bonds(i)%atomidx1 = atomidx1
         bonds(i)%atomidx2 = atomidx2
!         bonds(i)%typeidx = type_index(typestr)
      end do
   end if
end subroutine

end module
