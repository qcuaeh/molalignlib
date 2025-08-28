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

!> @defgroup clusrmsd ClusRMSD
!> @brief Program to calculate RMSDs between atom clusters
!> @{
program clusrmsd
use parameters
use options
use molecule
use spatial_transforms
use utils
use chemistry
use adjacency
use permutation
use file_path
use file_read
use file_write
use argparse
use biasing
use pruning
use registration
use assignment_cluster
use assignment_default

implicit none

character(:), allocatable :: title1, title2
character(:), allocatable :: arg, pathout, dummy
character(:), allocatable :: extin1, extin2, extout, extpipe
logical :: heavy_flag, mass_flag, align_flag, remap_flag, write_flag, pipe_flag, stats_flag
type(strlist_type) :: posargs(2)
type(atom_t), dimension(:), allocatable :: atoms1, atoms2
type(bond_t), dimension(:), allocatable :: bonds1, bonds2
type(partition_t) :: atomtypes
type(registry_t) :: registry
type(subperm_t) :: atomperm
real(rk) :: rmsd
real(rk) :: center1(3), center2(3), rotquat(4)
real(rk), dimension(:), allocatable :: weights1, weights2
real(rk), dimension(:,:), allocatable :: coords1, coords2, coords1w, coords2w, coords2r
type(bool_matrix), dimension(:), allocatable :: prunes
integer :: unitin1, unitin2, unitout
integer :: i

! Set default options

test_flag = .false.
stats_flag = .false.
mirror_flag = .false.
align_flag = .false.
remap_flag = .false.
write_flag = .false.
pipe_flag = .false.
mass_flag = .false.

max_records = 1
max_count = 10
max_trials = huge( max_trials)

prune_tol = 0.5
prune_procedure => prune_none

! Get user options

call init_args()

do while (get_arg(arg))
   select case (arg)
   case ('-align')
      align_flag = .true.
   case ('-remap')
      remap_flag = .true.
   case ('-near')
      prune_procedure => prune_none
   case ('-prune')
      prune_procedure => prune_rd
   case ('-tol')
      call read_optarg(arg, prune_tol)
   case ('-heavy')
      heavy_flag = .true.
   case ('-mass')
      mass_flag = .true.
   case ('-mirror')
      mirror_flag = .true.
   case ('-count')
      call read_optarg( arg, max_count)
   case ('-trials')
      call read_optarg( arg, max_trials)
   case ('-N')
      call read_optarg( arg, max_records)
   case ('-out')
      write_flag = .true.
      call read_optarg( arg, pathout)
   case ('-pipe')
      pipe_flag = .true.
      write_flag = .true.
      call read_optarg( arg, extpipe)
   case ('-stats')
      stats_flag = .true.
   case ('-test')
      test_flag = .true.
   case default
      call read_posarg( arg, posargs)
   end select
end do

if (pipe_flag) then
   extin1 = extpipe
   extin2 = extpipe
   extout = extpipe
   unitin1 = stdin
   unitin2 = stdin
   unitout = stdout
else
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
   if (write_flag) then
      call split_path( pathout, dummy, dummy, extout)
      call open2write( pathout, unitout)
   end if
end if

! Read coordinates
call readfile( unitin1, extin1, title1, atoms1, bonds1)
call readfile( unitin2, extin2, title2, atoms2, bonds2)

if (heavy_flag) then
   ! Include heavy atoms only
   call include_heavy_atoms( atoms1)
   call include_heavy_atoms( atoms2)
else
   ! Include all atoms
   atoms1%mask = .true.
   atoms2%mask = .true.
end if

! Collect atom types in a partition
call collect_atomtypes( atoms1, atoms2, atomtypes)

! Abort if there are conflicting atomic types
if (any(atomtypes%parts%num_items1 /= atomtypes%parts%num_items2)) then
   write (stderr, '(A)') 'Error: These molecules are not isomers'
   stop 1
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

   center1 = get_centroid( atoms1, weights1)
   center2 = get_centroid( atoms2, weights2)
   call translate_coords( coords2, center1 - center2)

   ! Get weighted-centered coordinates
   coords1w = get_weighted_coords( atoms1, weights1, center1)
   coords2w = get_weighted_coords( atoms2, weights2, center2)

   if (remap_flag) then

      ! Remap atoms to minimize the MSD
      call prune_procedure( atomtypes, atoms1, atoms2, prunes)
      call optimize_atomperm_cluster( coords1w, coords2w, atomtypes, prunes, registry)

      ! Print optimization stats
      if (stats_flag) then
         call print_records( registry)
      end if

      do i = 1, registry%num_records
         atomperm = registry%records(i)%atomperm
!         rotquat = registry%records(i)%rotquat
         rotquat = least_rotquat( atomperm, coords1w, coords2w)
         coords2r = rotated_coords( coords2w, rotquat)
         rmsd = sqrt( mean_sqdist( atomperm, weights1, coords1w, coords2r))

         write (stdout,'(A)') str( rmsd)

         if (write_flag) then
            title2 = 'RMSD=' // str( rmsd)
            coords2r = rotated_coords( coords2, rotquat, center1)
            call set_coords( atoms2, coords2r)
            call writefile( unitout, extout, title2, atoms2, bonds2, atomperm)
         end if
      end do

   else

      atomperm = default_atomperm( atoms1, atoms2)
      rotquat = least_rotquat( coords1w, coords2w)
      coords2r = rotated_coords( coords2w, rotquat)
      rmsd = sqrt( mean_sqdist( atomperm, weights1, coords1w, coords2r))

      write (stdout,'(A)') str( rmsd)

      if (write_flag) then
         title2 = 'RMSD=' // str( rmsd)
         coords2r = rotated_coords( coords2, rotquat, center1)
         call set_coords( atoms2, coords2r)
         call writefile( unitout, extout, title2, atoms2, bonds2, atomperm)
      end if

   end if

else

   ! Get weighted coordinates
   coords1w = get_weighted_coords( atoms1, weights1)
   coords2w = get_weighted_coords( atoms2, weights2)

   if (remap_flag) then
      call prune_procedure( atomtypes, atoms1, atoms2, prunes)
      call assign_atoms_pruned( atomtypes, coords1w, coords2w, prunes, atomperm)
      rmsd = sqrt( mean_sqdist( atomperm, weights1, coords1w, coords2w))
   else
      atomperm = default_atomperm( atoms1, atoms2)
      rmsd = sqrt( mean_sqdist( atomperm, weights1, coords1w, coords2w))
   end if

   write (stdout,'(A)') str( rmsd)

   if (write_flag) then
      title2 = 'RMSD=' // str( rmsd)
      coords2r = rotated_coords( coords2, rotquat, center1)
      call set_coords( atoms2, coords2r)
      call writefile( unitout, extout, title2, atoms2, bonds2, atomperm)
   end if

end if

end program
!> @}
