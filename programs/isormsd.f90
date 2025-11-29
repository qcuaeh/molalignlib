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

!> @defgroup isormsd IsoRMSD
!> @brief Program to calculate RMSDs between isomers
!> @{
program isormsd
use parameters
use molecule
use euclidean
use utils
use chemistry
use permutation
use adjacency
use file_reading
use file_writing
use argparse
use assorting
use biasing_isomer
use recording
use alignment_isomer
use assignment_conformer
use alignment_conformer
use pruning_atoms
use options
implicit none

character(:), allocatable :: title1, title2
character(:), allocatable :: arg, coords_path
character(:), allocatable :: typein, typeout
type(strlist_type) :: posargs(2)
type(atom_t), dimension(:), allocatable :: atoms1, atoms2
type(bond_t), dimension(:), allocatable :: bonds1, bonds2
type(adjc_t), dimension(:), allocatable :: adjcs1, adjcs2
type(adjc_t), dimension(:), allocatable :: adjcs1_mod, adjcs2_mod
type(partition_t) :: atomtypes
type(registry_t) :: temp_registry
type(registry_t), target :: iso_registry, confo_registry
type(registry_t), pointer :: output_registry
real(rk) :: center1(3), center2(3), rotquat(4)
real(rk), dimension(:), allocatable :: weights1, weights2
real(rk), dimension(:,:), allocatable :: coords1, coords2, coords1w, coords2w, coords2r
integer, dimension(:), allocatable :: atomset1, atomset2
integer, dimension(:), allocatable :: atomperm1
integer :: unitin, unitout
integer :: adjd
real(rk) :: rmsd
integer :: i

procedure(bond_modifier_interface), pointer :: match_diff_bonds => null()

! Set default options

stats_flag = .FALSE.
heavy_flag = .FALSE.
mirror_flag = .FALSE.
align_flag = .FALSE.
remap_flag = .FALSE.
coords_flag = .FALSE.
mass_flag = .FALSE.
stoch_flag = .TRUE.
adaptive_flag = .TRUE.
label_flag = .FALSE.
random_flag = .FALSE.
rebond_flag = .FALSE.
atomorder_flag = .FALSE.

num_records = 1
iso_thres = 100
confo_thres = 100
max_trials = MAX_TRIALS_DEFAULT
unitout = stdout

! Read command options

call init_args()

do while (get_arg(arg))
   select case (lowercase(arg))
   case ('-align')
      align_flag = .TRUE.
   case ('-remap')
      remap_flag = .TRUE.
   case ('-atomorder')
      atomorder_flag = .TRUE.
   case ('-label')
      label_flag = .TRUE.
   case ('-heavy')
      heavy_flag = .TRUE.
   case ('-mass')
      mass_flag = .TRUE.
   case ('-mirror')
      mirror_flag = .TRUE.
   case ('-thres')
      call read_optarg( arg, iso_thres)
   case ('-trials')
      call read_optarg( arg, max_trials)
   case ('-records')
      call read_optarg( arg, num_records)
   case ('-coords')
      coords_flag = .TRUE.
      call read_optarg( arg, coords_path)
   case ('-stats')
      stats_flag = .TRUE.
   case ('-random')
      random_flag = .TRUE.
   case ('-rebond')
      rebond_flag = .TRUE.
   case ('-confo1')
      match_diff_bonds => toggle_bonds1
   case ('-confo2')
      match_diff_bonds => toggle_bonds2
   case ('-confoadd')
      match_diff_bonds => add_missing_bonds
   case ('-confodel')
      match_diff_bonds => delete_extra_bonds
   case default
      call read_posarg( arg, posargs)
   end select
end do

select case (ipos)
case (0)
   write (stderr, '(A)') 'Error: Missing file paths'
   stop 1
case (1)
   write (stderr, '(A)') 'Error: Too few file paths'
   stop 1
case (2)
   call open2read( posargs(1)%arg, typein, unitin)
   call read_file( unitin, typein, title1, atoms1, bonds1)
   close (unitin)
   call open2read( posargs(2)%arg, typein, unitin)
   call read_file( unitin, typein, title2, atoms2, bonds2)
   close (unitin)
case default
   write (stderr, '(A)') 'Error: Too many file paths'
   stop 1
end select

if (coords_flag) then
   call parse_path( coords_path, typeout)
   call open2write( coords_path, unitout)
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

! Abort if atom types do not match
if (any(atomtypes%parts%num_items1 /= atomtypes%parts%num_items2)) then
   write (stderr, '(A)') 'Error: These molecules are not isomers'
   stop 1
end if

! Set adjacency lists
if (rebond_flag) then
   call adjacency_from_atoms( atomset1, atoms1, adjcs1)
   call adjacency_from_atoms( atomset2, atoms2, adjcs2)
else
   if (size(bonds1) < 1 .or. size(bonds2) < 1) then
      if (size(bonds1) < 1 .and. size(bonds2) < 1) then
         write (stdout,'(A)') 'Error: Molecules have no bonds!'
         stop 1
      else if (size(bonds1) < 1) then
         write (stdout,'(A)') 'Error: First molecule has no bonds!'
         stop 1
      else if (size(bonds2) < 1) then
         write (stdout,'(A)') 'Error: Second molecule has no bonds!'
         stop 1
      end if
   end if
   call adjacency_from_bonds( atomset1, bonds1, size(atoms1), adjcs1)
   call adjacency_from_bonds( atomset2, bonds2, size(atoms2), adjcs2)
end if

! Get user defined atom weights
if (mass_flag) then
   weights1 = atomic_masses(atoms1%elnum)
   weights2 = atomic_masses(atoms2%elnum)
else
   weights1 = uniform_weights( 1._rk, size(atoms1))
   weights2 = uniform_weights( 1._rk, size(atoms2))
end if

! Get mol1 coordinates
coords1 = get_coords(atoms1)

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

      ! Allocate registries
      call allocate_registry( iso_registry, num_records)

      ! Remap atoms to minimize adjacency difference and MSD
      call optimize_atomperm_isomer( atomset1, atomset2, atomtypes, adjcs1, adjcs2, &
            coords1w, coords2w, iso_registry)

      ! Print isomer optimization stats
      if (stats_flag) then
         call print_records( iso_registry)
      end if

      ! Check if conformer optimization is needed
      if (associated(match_diff_bonds)) then

         call allocate_registry( confo_registry, num_records)
         call allocate_registry( temp_registry, 1)

         do i = 1, iso_registry%occ_records

            atomperm1 = iso_registry%records(i)%atomperm1

            ! Apply bond modification strategy
            call match_diff_bonds( adjcs1, adjcs2, atomperm1, &
                  iso_registry%records(i)%moldiffs, adjcs1_mod, adjcs2_mod)
            call optimize_atomperm_conformer( atomset1, atomset2, adjcs1_mod, &
                  adjcs2_mod, atomtypes, coords1w, coords2w, temp_registry)

            call insert_record_atomperm( &
               confo_registry, &
               temp_registry%records(1)%atomperm1, &
               temp_registry%records(1)%steps, &
               temp_registry%records(1)%rotation, &
               adjacencydiff(atomset1, temp_registry%records(1)%atomperm1, adjcs1, adjcs2), &
               temp_registry%records(1)%permdist &
            )

         end do

         ! Print conformer optimization stats
         if (stats_flag) then
            call print_records( confo_registry)
         end if

         ! Point to conformer registry for output
         output_registry => confo_registry

      else

         ! Point to isomer registry for output
         output_registry => iso_registry

      end if

      ! Unified printing loop using stored rotations
      do i = 1, output_registry%occ_records
         atomperm1 = output_registry%records(i)%atomperm1
         rotquat = output_registry%records(i)%rotation
         coords2r = rotated_coords( coords2, rotquat, center1)
         rmsd = sqrt( sqdistmean( atomset1, atomperm1, weights1, coords1, coords2r))
         adjd = adjacencydiff( atomset1, atomperm1, adjcs1, adjcs2)

         if (coords_flag) then
            title2 = 'RMSD=' // str(rmsd) // ' Δadj=' // str(adjd)
            call set_coords( atoms2, coords2r)
            call writefile( unitout, typeout, title2, atoms2, bonds2, atomperm1)
         else
            write (stdout,'(A,"(",I0,")")',advance='no') str(rmsd), adjd
            if (atomorder_flag) then
               write (stdout,'(1X)',advance='no')
               call print_permutation( atomperm1)
            end if
            write (stdout, *)
         end if
      end do

   else

      allocate (atomperm1(size( coords1w, 2)))
      call init_identity_permutation( atomperm1)
      rotquat = least_rotquat( atomset1, atomperm1, coords1w, coords2w)
      coords2r = rotated_coords( coords2, rotquat, center1)
      rmsd = sqrt(sqdistmean( atomset1, atomperm1, weights1, coords1, coords2r))
      adjd = adjacencydiff( atomset1, atomperm1, adjcs1, adjcs2)

      if (coords_flag) then
         title2 = 'RMSD=' // str(rmsd) // ' Δadj=' // str(adjd)
         call set_coords( atoms2, coords2r)
         call writefile( unitout, typeout, title2, atoms2, bonds2, atomperm1)
      else
         write (stdout,'(A,"(",I0,")")') str(rmsd), adjd
      end if

   end if

else

   ! Get weighted coordinates
   coords1w = get_weighted_coords( atoms1, weights1)
   coords2w = get_weighted_coords( atoms2, weights2)

   if (remap_flag) then
      write (stdout, '(A)') 'Error: Remapping without alignment is not implemented'
      stop 1
   else
      allocate (atomperm1(size(coords1w, 2)))
      call init_identity_permutation( atomperm1)
      rmsd = sqrt(sqdistmean( atomset1, atomperm1, weights1, coords1, coords2))
      adjd = adjacencydiff( atomset1, atomperm1, adjcs1, adjcs2)
   end if

   if (coords_flag) then
      title2 = 'RMSD=' // str(rmsd) // ' Δadj=' // str(adjd)
      call set_coords( atoms2, coords2)
      call writefile( unitout, typeout, title2, atoms2, bonds2, atomperm1)
   else
      write (stdout,'(A,A,I0,A)') str(rmsd), '(', adjd, ')'
   end if

end if

end program
!> @}
