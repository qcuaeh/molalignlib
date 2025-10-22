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

!> @defgroup conformsd ConfoRMSD
!> @brief Program to calculate RMSDs between conformers
!> @{
program conformsd
use parameters
use molecule
use euclidean
use utils
use chemistry
use permutation
use file_path
use file_read
use file_write
use argparse
use biasing_isomer
use pruning_atoms
use recording
use assignment_conformer
use alignment_conformer
use options
implicit none

character(:), allocatable :: title1, title2
character(:), allocatable :: arg, fileout_path, dummy
character(:), allocatable :: extin1, extin2, extout
type(strlist_type) :: posargs(2)
type(atom_t), dimension(:), allocatable :: atoms1, atoms2
type(bond_t), dimension(:), allocatable :: bonds1, bonds2
type(adjc_t), dimension(:), allocatable :: adjcs1, adjcs2
type(partition_t) :: atomtypes
type(registry_t) :: registry
real(rk) :: rmsd
real(rk) :: center1(3), center2(3), rotquat(4)
real(rk), dimension(:), allocatable :: weights1, weights2
real(rk), dimension(:,:), allocatable :: coords1, coords2, coords1w, coords2w, coords2r
integer, dimension(:), pointer :: atomset1, atomset2
integer, dimension(:), allocatable :: atomset1_alloc, atomset2_alloc
integer, dimension(:), allocatable :: atomperm1
integer :: unitin1, unitin2, unitout
integer :: i

! Set default options

stats_flag = .false.
heavy_flag = .false.
mirror_flag = .false.
align_flag = .false.
remap_flag = .false.
coords_flag = .false.
tree_flag = .false.
mass_flag = .false.
stoch_flag = .true.
adaptive_flag = .true.
bond_flag = .false.
label_flag = .false.
random_flag = .false.
full_flag = .false.
permutation_flag = .false.

num_records = 1
count_thres = 100
unitout = stdout
max_trials = huge( ik)

! Read command options

call init_args()

do while (get_arg(arg))
   select case (lowercase(arg))
   case ('-align')
      align_flag = .true.
   case ('-remap')
      remap_flag = .true.
   case ('-permutation')
      permutation_flag = .true.
   case ('-full')
      full_flag = .true.
   case ('-exhaustive')
      stoch_flag = .false.
   case ('-stochastic')
      stoch_flag = .true.
      adaptive_flag = .false.
   case ('-label')
      label_flag = .true.
   case ('-heavy')
      heavy_flag = .true.
   case ('-mass')
      mass_flag = .true.
   case ('-mirror')
      mirror_flag = .true.
   case ('-count')
      call read_optarg( arg, count_thres)
   case ('-trials')
      call read_optarg( arg, max_trials)
   case ('-records')
      call read_optarg( arg, num_records)
   case ('-coords')
      coords_flag = .true.
      call read_optarg( arg, fileout_path)
   case ('-tree')
      tree_flag = .true.
   case ('-stats')
      stats_flag = .true.
   case ('-random')
      random_flag = .true.
   case ('-bond')
      bond_flag = .true.
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
   call split_path( posargs(1)%arg, dummy, dummy, extin1)
   call split_path( posargs(2)%arg, dummy, dummy, extin2)
   call open2read( posargs(1)%arg, unitin1)
   call open2read( posargs(2)%arg, unitin2)
case default
   write (stderr, '(A)') 'Error: Too many file paths'
   stop 1
end select

if (coords_flag) then
   call split_path( fileout_path, dummy, dummy, extout)
   call open2write( fileout_path, unitout)
end if

! Read coordinates
call readfile( unitin1, extin1, title1, atoms1, bonds1)
call readfile( unitin2, extin2, title2, atoms2, bonds2)

if (heavy_flag) then
   ! Include only heavy atoms
   call include_heavy_atoms( atoms1, atomset1, atomset1_alloc)
   call include_heavy_atoms( atoms2, atomset2, atomset2_alloc)
else
   ! Include all atoms
   call include_all_atoms( atoms1, atomset1, atomset1_alloc)
   call include_all_atoms( atoms2, atomset2, atomset2_alloc)
end if

! Collect atom types in a partition
call collect_atomtypes( atomset1, atomset2, atoms1, atoms2, atomtypes)

! Abort if atom types do not match
if (any(atomtypes%parts%num_items1 /= atomtypes%parts%num_items2)) then
   write (stderr, '(A)') 'Error: These molecules are not isomers'
   stop 1
end if

if (bond_flag) then
   call adjacency_from_distance( atomset1, atoms1, adjcs1)
   call adjacency_from_distance( atomset2, atoms2, adjcs2)
else
   if (size(bonds1) < 1 .or. size(bonds2) < 1) then
      if (size(bonds1) < 1 .and. size(bonds2) < 1) then
         write (stdout,'(A)') 'Error: Molecules have no bonds!'
         stop
      else if (size(bonds1) < 1) then
         write (stdout,'(A)') 'Error: First molecule has no bonds!'
         stop
      else if (size(bonds2) < 1) then
         write (stdout,'(A)') 'Error: Second molecule has no bonds!'
         stop
      end if
   end if
   call adjacency_from_bonds( atomset1, atoms1, bonds1, adjcs1)
   call adjacency_from_bonds( atomset2, atoms2, bonds2, adjcs2)
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
      call allocate_registry( registry, num_records)
      call optimize_atomperm_conformer( atomset1, atomset2, adjcs1, adjcs2, atomtypes, &
            coords1w, coords2w, registry)

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

         if (coords_flag) then
            title2 = 'RMSD=' // str( rmsd)
            coords2r = rotated_coords( coords2, rotquat, center1)
            call set_coords( atoms2, coords2r)
            call writefile( unitout, extout, title2, atoms2, bonds2, atomperm1)
         else
            write (stdout,'(A)') str( rmsd)
            if (permutation_flag) then
               call print_permutation(atomperm1)
            end if
         end if

      end do

   else

      allocate (atomperm1(size( coords1w, 2)))
      call init_identity_permutation( atomperm1)
      rotquat = least_rotquat( atomset1, atomperm1, coords1w, coords2w)
      coords2r = rotated_coords( coords2, rotquat, center1)
      rmsd = sqrt( sqdistmean( atomset1, atomperm1, weights1, coords1, coords2r))

      if (coords_flag) then
         title2 = 'RMSD=' // str( rmsd)
         coords2r = rotated_coords( coords2, rotquat, center1)
         call set_coords( atoms2, coords2r)
         call writefile( unitout, extout, title2, atoms2, bonds2, atomperm1)
      else
         write (stdout,'(A)') str( rmsd)
      end if

   end if

else

   ! Get weighted coordinates
   coords1w = get_weighted_coords( atoms1, weights1)
   coords2w = get_weighted_coords( atoms2, weights2)

   if (remap_flag) then
      call assign_atomperm_conformer( adjcs1, adjcs2, atomtypes, coords1w, coords2w, atomperm1)
      rmsd = sqrt( sqdistmean( atomset1, atomperm1, weights1, coords1, coords2))
   else
      allocate (atomperm1(size( coords1w, 2)))
      call init_identity_permutation( atomperm1)
      rmsd = sqrt( sqdistmean( atomset1, atomperm1, weights1, coords1, coords2))
   end if

   if (coords_flag) then
      title2 = 'RMSD=' // str( rmsd)
      coords2r = rotated_coords( coords2, rotquat, center1)
      call set_coords( atoms2, coords2r)
      call writefile( unitout, extout, title2, atoms2, bonds2, atomperm1)
   else
      write (stdout,'(A)') str( rmsd)
      if (remap_flag .and. permutation_flag) then
         call print_permutation(atomperm1)
      end if
   end if

end if

end program
!> @}
