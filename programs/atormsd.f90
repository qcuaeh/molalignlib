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
use molecule
use euclidean
use utils
use chemistry
use adjacency
use permutation
use file_reading
use file_writing
use argparse
use assorting
use biasing_isomer
use pruning_atoms
use recording
use alignment_atoms
use options
implicit none

character(:), allocatable :: title1, title2
character(:), allocatable :: arg, coords_path
character(:), allocatable :: typein, typeout
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
integer, dimension(:), allocatable :: atomset1, atomset2
integer, dimension(:), allocatable :: atomperm1
integer :: unitin, unitout
integer :: i

! Set default options

stats_flag = .FALSE.
mirror_flag = .FALSE.
align_flag = .FALSE.
remap_flag = .FALSE.
test_flag = .FALSE.
aligned_flag = .FALSE.
write_aligned = .FALSE.
mass_flag = .FALSE.
label_flag = .FALSE.
random_flag = .FALSE.
atomorder_flag = .FALSE.

num_records = 1
ato_thres = 10
max_trials = MAX_TRIALS_DEFAULT
unitout = stdout
prune_procedure => prune_none

! Get user options

call init_args()

do while (get_arg(arg))
   select case (lowercase(arg))
   case ('-align')
      align_flag = .TRUE.
   case ('-remap')
      remap_flag = .TRUE.
   case ('-atomorder')
      atomorder_flag = .TRUE.
   case ('-near')
      prune_procedure => prune_none
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
   case ('-thres')
      call read_optarg(arg, ato_thres)
   case ('-trials')
      call read_optarg( arg, max_trials)
   case ('-records')
      call read_optarg( arg, num_records)
   case ('-aligned')
      aligned_flag = .TRUE.
      call read_optarg( arg, coords_path)
   case ('-stats')
      stats_flag = .TRUE.
   case ('-random')
      random_flag = .TRUE.
   case ('-test')
      test_flag = .TRUE.
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
   call open2read( posargs(1)%arg, typein, unitin)
   call read_file( unitin, typein, title1, atoms1, bonds1)
   close (unitin)
   call open2read( posargs(2)%arg, typein, unitin)
   call read_file( unitin, typein, title2, atoms2, bonds2)
   close (unitin)
case default
   stop 'Too many file paths'
end select

if (aligned_flag) then
   write_aligned = .TRUE.
   call parse_path( coords_path, typeout)
   call open2write( coords_path, unitout)
end if

if (test_flag) then
   write_aligned = .TRUE.
   typeout = 'xyz'
   unitout = stdout
end if

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
   weights1 = uniform_weights( 1._rk, size(atoms1))
   weights2 = uniform_weights( 1._rk, size(atoms2))
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

   center1 = get_centroid( atomset1, atoms1, weights1)
   center2 = get_centroid( atomset2, atoms2, weights2)
   call translate_coords( coords2, center1 - center2)

   ! Get weighted-centered coordinates
   coords1w = get_weighted_coords( atoms1, weights1, center1)
   coords2w = get_weighted_coords( atoms2, weights2, center2)

   if (remap_flag) then

      ! Remap atoms to minimize the MSD
      call prune_procedure( atomtypes, coords1, coords2, prunes)
      call allocate_registry( registry, num_records)
      call optimize_atomperm_atoms( atomset1, atomset2, atomtypes, prunes, coords1w, coords2w, registry)

      ! Print optimization stats
      if (stats_flag) then
         call print_records( registry)
      end if

      do i = 1, registry%occ_records
         atomperm1 = registry%records(i)%atomperm1
!         rotquat = registry%records(i)%rotquat
         rotquat = least_rotquat( atomset1, atomperm1, coords1w, coords2w)
         coords2r = rotated_coords( coords2, rotquat, center1)
         rmsd = sqrt( sqdistmean( atomset1, atomperm1, weights1, coords1, coords2r))

         if (write_aligned) then
            title2 = 'RMSD=' // str( rmsd)
            coords2r = rotated_coords( coords2, rotquat, center1)
            call set_coords( atoms2, coords2r)
            call writefile( unitout, typeout, title2, atoms2, bonds2, atomperm1)
         else
            write (stdout,'(A)',advance='no') str( rmsd)
            if (atomorder_flag) then
               write (stdout,'(1X)',advance='no')
               call print_permutation(atomperm1)
            end if
            write (stdout, *)
         end if
      end do

   else

      ! Abort if atom types do not match
      if (any(atomtypes%itemdir1 /= atomtypes%itemdir2)) then
         stop 'Atom types do not match'
      end if

      allocate (atomperm1(size( coords1w, 2)))
      call init_identity_permutation( atomperm1)
      rotquat = least_rotquat( coords1w, coords2w)
      coords2r = rotated_coords( coords2, rotquat, center1)
      rmsd = sqrt( sqdistmean( atomset1, atomperm1, weights1, coords1, coords2r))

      if (write_aligned) then
         title2 = 'RMSD=' // str( rmsd)
         coords2r = rotated_coords( coords2, rotquat, center1)
         call set_coords( atoms2, coords2r)
         call writefile( unitout, typeout, title2, atoms2, bonds2, atomperm1)
      else
         write (unitout,'(A)') str( rmsd)
      end if

   end if

else

   ! Get weighted coordinates
   coords1w = get_weighted_coords( atoms1, weights1)
   coords2w = get_weighted_coords( atoms2, weights2)

   if (remap_flag) then
      call prune_procedure( atomtypes, coords1, coords2, prunes)
      call assign_atoms_pruned( atomtypes, coords1w, coords2w, prunes, atomperm1)
      rmsd = sqrt( sqdistmean( atomset1, atomperm1, weights1, coords1, coords2))
   else
      allocate (atomperm1(size( coords1w, 2)))
      call init_identity_permutation( atomperm1)
      rmsd = sqrt( sqdistmean( atomset1, atomperm1, weights1, coords1, coords2))
   end if

   if (write_aligned) then
      title2 = 'RMSD=' // str( rmsd)
      coords2r = rotated_coords( coords2, rotquat, center1)
      call set_coords( atoms2, coords2r)
      call writefile( unitout, typeout, title2, atoms2, bonds2, atomperm1)
   else
      write (stdout,'(A)',advance='no') str( rmsd)
      if (remap_flag .and. atomorder_flag) then
         write (stdout,'(1X)',advance='no')
         call print_permutation(atomperm1)
      end if
      write (stdout, *)
   end if

end if

end program
!> @}
