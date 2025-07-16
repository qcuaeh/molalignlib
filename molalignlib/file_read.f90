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

subroutine readmol(unit, format, title, atoms)
   integer, intent(in) :: unit
   character(*), intent(in) :: format
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms

   select case (format)
   case ('xyz')
      call readmol_xyz(unit, title, atoms)
   case ('mol')
      call readmol_mol(unit, title, atoms)
   case ('sdf')
      call readmol_sdf(unit, title, atoms)
   case ('mol2')
      call readmol_mol2(unit, title, atoms)
   case default
      write (stderr, '(A,1X,A)') 'Invalid format:', format
      stop
   end select
end subroutine

subroutine readmol_xyz(unit, title, atoms)
   integer, intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   ! Local variables
   character(ll) :: buffer
   character(wl) :: elsym
   integer :: i, num_atoms, elnum, label, stat
   real(rk) :: coords(3)

   ! Read number of atoms
   read (unit, *, iostat=stat) num_atoms
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid XYZ format'
      stop
   end if

   allocate (atoms(num_atoms))

   ! Read title line
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid XYZ format'
      stop
   end if
   title = trim(buffer)

   ! Read atom lines
   do i = 1, num_atoms
      read (unit, *, iostat=stat) elsym, coords
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid XYZ format'
         stop
      end if

      call split_symbol(elsym, elnum, label)
      atoms(i)%elnum = elnum
      atoms(i)%weight = atomic_weights(elnum)
      atoms(i)%label = label
      atoms(i)%coords = coords
   end do

   if (adjacency_flag) then
      if (bond_flag) then
         call bond_atoms(atoms)
      else
         write (stderr, '(A)') 'There are no bonds in this file, use the -bond option to use default bonds'
         stop
      end if
   end if
end subroutine

subroutine readmol_mol(unit, title, atoms)
   integer, intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, target, intent(out) :: atoms
   ! Local variables
   real(rk) :: coords(3)
   character(3) :: elsym
   character(ll) :: buffer
   integer :: num_atoms, num_bonds, elnum, label
   integer :: i, idx1, idx2, bond_type, stat

   ! Read header block (3 lines)
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL format'
      stop
   end if
   title = trim(buffer)

   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL format'
      stop
   end if

   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL format'
      stop
   end if

   ! Read counts line
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL format'
      stop
   end if

   read (buffer(1:3), '(I3)', iostat=stat) num_atoms
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL format'
      stop
   end if
   read (buffer(4:6), '(I3)', iostat=stat) num_bonds
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL format'
      stop
   end if

   allocate (atoms(num_atoms))

   ! Read atom block
   do i = 1, num_atoms
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL format'
         stop
      end if

      read (buffer(1:10), '(F10.4)', iostat=stat) coords(1)
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL format'
         stop
      end if
      read (buffer(11:20), '(F10.4)', iostat=stat) coords(2)
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL format'
         stop
      end if
      read (buffer(21:30), '(F10.4)', iostat=stat) coords(3)
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL format'
         stop
      end if

      elsym = adjustl(buffer(32:34))

      call split_symbol(elsym, elnum, label)
      atoms(i)%elnum = elnum
      atoms(i)%weight = atomic_weights(elnum)
      atoms(i)%label = label
      atoms(i)%coords = coords
   end do

   ! Read bond block if adjacency is needed
   if (adjacency_flag) then
      if (bond_flag) then
         call bond_atoms(atoms)
      else
         if (num_bonds > 0) then
         block
            integer, allocatable :: nadjs(:)
            allocate (nadjs(num_atoms))

            nadjs = 0
            do i = 1, num_bonds
               read (unit, '(A)', iostat=stat) buffer
               if (stat /= 0) then
                  write (stderr, '(A)') 'Invalid MOL format'
                  stop
               end if

               read (buffer(1:3), '(I3)', iostat=stat) idx1
               if (stat /= 0) then
                  write (stderr, '(A)') 'Invalid MOL format'
                  stop
               end if
               read (buffer(4:6), '(I3)', iostat=stat) idx2
               if (stat /= 0) then
                  write (stderr, '(A)') 'Invalid MOL format'
                  stop
               end if

               ! Read bond type (characters 7-9)
               if (len_trim(buffer) >= 9) then
                  read (buffer(7:9), '(I3)', iostat=stat) bond_type
                  if (stat /= 0) bond_type = 1
               else
                  bond_type = 1
               end if
               ! bond_type is now available in local variable for future use

               nadjs(idx1) = nadjs(idx1) + 1
               nadjs(idx2) = nadjs(idx2) + 1
               atoms(idx1)%adjlist_allocation(nadjs(idx1)) = idx2
               atoms(idx2)%adjlist_allocation(nadjs(idx2)) = idx1
            end do

            do i = 1, num_atoms
               atoms(i)%adjlist => atoms(i)%adjlist_allocation(1:nadjs(i))
            end do
         end block
         else if (num_bonds == 0) then
            write (stderr, '(A)') 'There are no bonds in this file, use the -bond option to use default bonds'
            stop
         else
            write (stderr, '(A)') 'Invalid MOL2 format'
            stop
         end if
      end if
   end if
end subroutine

subroutine readmol_sdf(unit, title, atoms)
   integer, intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   ! Local variables
   character(ll) :: buffer
   logical :: end_of_molecule
   integer :: stat

   ! Read the MOL block (SDF files contain MOL blocks followed by data)
   call readmol_mol(unit, title, atoms)

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

subroutine readmol_mol2(unit, title, atoms)
   integer, intent(in) :: unit
   character(:), allocatable, intent(out) :: title
   type(atom_t), dimension(:), allocatable, target, intent(out) :: atoms
   ! Local variables
   real(rk) :: coords(3)
   character(wl) :: elsym, dummy, atom_type, bond_type
   integer :: num_atoms, num_bonds, elnum, label, stat
   character(ll) :: buffer
   integer :: i, idx1, idx2

   ! Find @<TRIPOS>MOLECULE section
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL2 format'
         stop
      end if
      if (trim(buffer) == '@<TRIPOS>MOLECULE') exit
   end do

   ! Read molecule name
   read (unit, '(A)', iostat=stat) buffer
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL2 format'
      stop
   end if
   title = trim(buffer)

   ! Read counts line
   read (unit, *, iostat=stat) num_atoms, num_bonds
   if (stat /= 0) then
      write (stderr, '(A)') 'Invalid MOL2 format'
      stop
   end if

   allocate (atoms(num_atoms))

   ! Find @<TRIPOS>ATOM section
   do
      read (unit, '(A)', iostat=stat) buffer
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL2 format'
         stop
      end if
      if (trim(buffer) == '@<TRIPOS>ATOM') exit
   end do

   ! Read atom data
   ! Format: atom_id atom_name x y z atom_type [subst_id subst_name charge]
   do i = 1, num_atoms
      read (unit, *, iostat=stat) dummy, elsym, coords, atom_type
      if (stat /= 0) then
         write (stderr, '(A)') 'Invalid MOL2 format'
         stop
      end if

      call split_symbol(elsym, elnum, label)
      atoms(i)%elnum = elnum
      atoms(i)%weight = atomic_weights(elnum)
      atoms(i)%label = label
      atoms(i)%coords = coords
      ! atom_type is now available in local variable for future use
   end do

   ! Handle bonds if needed
   if (adjacency_flag) then
      if (bond_flag) then
         call bond_atoms(atoms)
      else
         if (num_bonds > 0) then
         block
            integer, allocatable :: nadjs(:)
            allocate (nadjs(num_atoms))

            ! Find @<TRIPOS>BOND section
            do
               read (unit, '(A)', iostat=stat) buffer
               if (stat /= 0) then
                  write (stderr, '(A)') 'Invalid MOL2 format'
                  stop
               end if
               if (trim(buffer) == '@<TRIPOS>BOND') exit
            end do

            ! Read bond data
            ! Format: bond_id origin_atom_id target_atom_id bond_type
            nadjs = 0
            do i = 1, num_bonds
               read (unit, *, iostat=stat) dummy, idx1, idx2, bond_type
               if (stat /= 0) then
                  write (stderr, '(A)') 'Invalid MOL2 format'
                  stop
               end if

               nadjs(idx1) = nadjs(idx1) + 1
               nadjs(idx2) = nadjs(idx2) + 1
               atoms(idx1)%adjlist_allocation(nadjs(idx1)) = idx2
               atoms(idx2)%adjlist_allocation(nadjs(idx2)) = idx1
               ! bond_type is now available in local variable for future use
            end do

            do i = 1, num_atoms
               atoms(i)%adjlist => atoms(i)%adjlist_allocation(1:nadjs(i))
            end do
         end block
         else if (num_bonds == 0) then
            write (stderr, '(A)') 'There are no bonds in this file, use the -bond option to use default bonds'
            stop
         else
            write (stderr, '(A)') 'Invalid MOL2 format'
            stop
         end if
      end if
   end if
end subroutine

end module
