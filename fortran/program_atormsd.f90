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

!> @defgroup atormsd AtoRMSD
!> @brief Program to calculate RMSDs between atom clusters
!> @{
program atormsd
use parameters
use str_utils
use chemdata
use molecule
use euclidean
use permutation
use assorting
use pruning_atoms
use recording
use alignment_atoms
use file_utils
use file_reading
use file_writing
use arg_parsing
use flags
implicit none

logical(lk) :: print_assignment
logical(lk) :: write_aligned
character(:), allocatable :: title1, title2
character(:), allocatable :: arg, aligned_path
character(:), allocatable :: in_format, out_format
type(strlist_type) :: posargs(2)
type(atom_t), dimension(:), allocatable :: atoms1, atoms2
type(bond_t), dimension(:), allocatable :: bonds1, bonds2
type(bool_matrix), dimension(:), allocatable :: prunes
type(partition_t) :: atomtypes
type(registry_t) :: registry
real(rk) :: rmsd
real(rk) :: center1(3), center2(3), rotquat(4)
real(rk), dimension(:), allocatable :: weights1, weights2
real(rk), dimension(:,:), allocatable :: coords1, coords2, coords1w, coords2w, coords2r
integer(ik), dimension(:), allocatable :: atomset1, atomset2
integer(ik), dimension(:), allocatable :: atomperm1
integer(ik) :: num_records, max_trials, conv_freq
integer(ik) :: in_unit, aligned_unit
integer(ik) :: error_code
integer(ik) :: i

! Set default options

heavy_flag = .FALSE.
mirror_flag = .FALSE.
align_flag = .FALSE.
remap_flag = .FALSE.
mass_flag = .FALSE.
label_flag = .FALSE.
random_flag = .FALSE.
print_stats = .FALSE.
print_assignment = .FALSE.
write_aligned = .FALSE.

num_records = 1
conv_freq = 10
max_trials = MAX_TRIALS_DEFAULT
prune_procedure => prune_none

! Get user options

call init_args()

do while (get_arg(arg))
   select case (lowercase(arg))
   case ('-align')
      align_flag = .TRUE.
   case ('-remap')
      remap_flag = .TRUE.
   case ('-prune')
      prune_procedure => prune_rd
      call read_optarg(arg, prune_tol)
   case ('-label')
      label_flag = .TRUE.
   case ('-heavy')
      heavy_flag = .TRUE.
   case ('-mass')
      mass_flag = .TRUE.
   case ('-mirror')
      mirror_flag = .TRUE.
   case ('-freq')
      call read_optarg(arg, conv_freq)
   case ('-trials')
      call read_optarg( arg, max_trials)
   case ('-records')
      call read_optarg( arg, num_records)
   case ('-aligned')
      write_aligned = .TRUE.
      call read_optarg( arg, aligned_path)
   case ('-assignment')
      print_assignment = .TRUE.
   case ('-stats')
      print_stats = .TRUE.
   case ('-random')
      random_flag = .TRUE.
   case default
      call read_posarg( arg, posargs)
   end select
end do

select case (ipos)
case (0)
   stop 'File paths are missing'
case (1)
   stop 'Too few file paths'
case (2)
   call open2read( posargs(1)%arg, in_format, in_unit)
   call read_file( in_unit, in_format, title1, atoms1, bonds1)
   close (in_unit)
   call open2read( posargs(2)%arg, in_format, in_unit)
   call read_file( in_unit, in_format, title2, atoms2, bonds2)
   close (in_unit)
case default
   stop 'Too many file paths'
end select

if (heavy_flag) then
   ! Include only heavy atoms
   call include_heavy_atoms( atoms1, atomset1)
   call include_heavy_atoms( atoms2, atomset2)
else
   ! Include all atoms
   call include_all_atoms( atoms1, atomset1)
   call include_all_atoms( atoms2, atomset2)
end if

! Collect atom types in a partition
call collect_atomtypes( atomset1, atomset2, atoms1, atoms2, atomtypes)

! Abort if molecules are not isomers
if (any(atomtypes%parts%num_items1 /= atomtypes%parts%num_items2)) then
   stop 'These molecules are not isomers'
end if

if (mass_flag) then
   weights1 = atomic_masses(atoms1%elnum)
   weights2 = atomic_masses(atoms2%elnum)
else
   allocate (weights1(size(atoms1)), source=1.0_rk)
   allocate (weights2(size(atoms2)), source=1.0_rk)
end if

! Get mol1 coordinates
coords1 = get_coords( atoms1)

! Get mol2 coordinates
if (mirror_flag) then
   coords2 = get_mirrored_coords( atoms2)
else
   coords2 = get_coords( atoms2)
end if

if (align_flag) then

   if (write_aligned) then
      call open2write( aligned_path, out_format, aligned_unit)
   end if

   center1 = get_centroid( atomset1, atoms1, weights1)
   center2 = get_centroid( atomset2, atoms2, weights2)
   call translate_coords( coords2, center1 - center2)

   coords1w = get_weighted_coords( atoms1, weights1, center1)
   coords2w = get_weighted_coords( atoms2, weights2, center2)

else

   coords1w = get_weighted_coords( atoms1, weights1)
   coords2w = get_weighted_coords( atoms2, weights2)

end if

if (remap_flag) then

   call prune_procedure( atomtypes, coords1, coords2, prunes)

   if (align_flag) then

      call allocate_registry( registry, num_records)
      call optimize_atomperm_atoms( atomset1, atomset2, atomtypes, prunes, &
            coords1w, coords2w, conv_freq, max_trials, registry, error_code)
      if (error_code /= 0) stop 'Error: Assignment failed'

      if (print_stats) call print_records( registry)

      do i = 1, registry%occ_records
         atomperm1 = registry%records(i)%atomperm1
         rotquat = least_rotquat( atomset1, atomperm1, coords1w, coords2w)
         coords2r = rotated_coords( coords2, rotquat, center1)
         rmsd = sqrt( sqdistmean( atomset1, atomperm1, weights1, coords1, coords2r))

         if (write_aligned) then
            title2 = 'rmsd=' // str( rmsd)
            call set_coords( atoms2, coords2r)
            call write_file( aligned_unit, out_format, title2, atoms2, bonds2, atomperm1)
         else
            write (stdout,'(A)',advance='no') str( rmsd)
            if (print_assignment) then
               write (stdout,'(1X)',advance='no')
               call print_permutation(atomperm1)
            end if
            write (stdout,*)
         end if
      end do

   else

      call assign_atoms_pruned( atomtypes, coords1w, coords2w, prunes, atomperm1, error_code)
      if (error_code /= 0) stop 'Error: Assignment failed'

      coords2r = coords2
      rmsd = sqrt( sqdistmean( atomset1, atomperm1, weights1, coords1, coords2r))
      write (stdout,'(A)',advance='no') str( rmsd)
      if (print_assignment) then
         write (stdout,'(1X)',advance='no')
         call print_permutation(atomperm1)
      end if
      write (stdout,*)

   end if

else

   ! Maintain input atom order
   call init_array( atomperm1, size(atoms1), identity)

   if (any(atomtypes%itemdir1 /= atomtypes%itemdir2)) then
      stop 'Atoms do not match'
   end if

   if (align_flag) then
      rotquat = least_rotquat( atomset1, atomperm1, coords1w, coords2w)
      coords2r = rotated_coords( coords2, rotquat, center1)
   else
      coords2r = coords2
   end if

   rmsd = sqrt( sqdistmean( atomset1, atomperm1, weights1, coords1, coords2r))

   if (align_flag .and. write_aligned) then
      title2 = 'rmsd=' // str( rmsd)
      call set_coords( atoms2, coords2r)
      call write_file( aligned_unit, out_format, title2, atoms2, bonds2, atomperm1)
   else
      write (stdout,'(A)') str( rmsd)
   end if

end if

end program
!> @}
