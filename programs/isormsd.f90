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
use file_path
use file_read
use file_write
use argparse
use biasing_isomer
use recording
use alignment_isomer
use assignment_conformer
use alignment_conformer
use pruning_atoms
use options
implicit none

character(:), allocatable :: title1, title2
character(:), allocatable :: arg, fileout_path, dummy
character(:), allocatable :: extin1, extin2, extout, extin
type(strlist_type) :: posargs(2)
type(atom_t), dimension(:), allocatable :: atoms1, atoms2
type(bond_t), dimension(:), allocatable :: bonds1, bonds2
type(adjc_t), dimension(:), allocatable :: adjcs1, adjcs2
type(adjc_t), dimension(:), allocatable :: adjcs1_mod, adjcs2_mod
logical, dimension(:,:), allocatable :: adjmat1, adjmat2
type(partition_t) :: atomtypes
type(registry_t) :: iso_registry, confo_registry, temp_registry
real(rk) :: center1(3), center2(3), rotquat(4)
real(rk), dimension(:), allocatable :: weights1, weights2
real(rk), dimension(:,:), allocatable :: coords1, coords2, coords1w, coords2w, coords2r
integer, dimension(:), pointer :: atomset1, atomset2
integer, dimension(:), allocatable :: atomset1_alloc, atomset2_alloc
integer, dimension(:), allocatable :: atomperm1
integer :: unitin1, unitin2, unitout
integer :: adjd
real(rk) :: rmsd
integer :: i
logical :: union_flag
logical :: intersection_flag

! Set default options

stats_flag = .false.
heavy_flag = .false.
mirror_flag = .false.
align_flag = .false.
remap_flag = .false.
coords_flag = .false.
stdin_flag = .false.
mass_flag = .false.
stoch_flag = .true.
adaptive_flag = .true.
mapping_flag = .false.
label_flag = .false.
random_flag = .false.
iterate_flag = .true.
bond_flag = .false.
union_flag = .false.
intersection_flag = .false.

extin = 'xyz'
extout = 'xyz'
num_records = 10
count_thres = 100
unitout = stdout
max_trials = huge(ik)

! Read command options

call init_args()

do while (get_arg(arg))
   select case (lowercase(arg))
   case ('-align')
      align_flag = .true.
   case ('-remap')
      remap_flag = .true.
   case ('-mapping')
      mapping_flag = .true.
   case ('-label')
      label_flag = .true.
   case ('-heavy')
      heavy_flag = .true.
   case ('-mass')
      mass_flag = .true.
   case ('-mirror')
      mirror_flag = .true.
   case ('-count')
      call read_optarg(arg, count_thres)
   case ('-trials')
      call read_optarg(arg, max_trials)
   case ('-records')
      call read_optarg(arg, num_records)
   case ('-coords')
      coords_flag = .true.
   case ('-out')
      fileout_flag = .true.
      call read_optarg(arg, fileout_path)
   case ('-stdin')
      stdin_flag = .true.
   case ('-extin')
      call read_optarg(arg, extin)
   case ('-stats')
      stats_flag = .true.
   case ('-random')
      random_flag = .true.
   case ('-bond')
      bond_flag = .true.
   case ('-union')
      union_flag = .true.
   case ('-intersection')
      intersection_flag = .true.
   case default
      call read_posarg(arg, posargs)
   end select
end do

if (stdin_flag) then
   extin1 = extin
   extin2 = extin
   unitin1 = stdin
   unitin2 = stdin
else
   select case (ipos)
   case (0)
      write (stderr, '(A)') 'Error: Missing file paths'
      stop 1
   case (1)
      write (stderr, '(A)') 'Error: Too few file paths'
      stop 1
   case (2)
      call split_path(posargs(1)%arg, dummy, dummy, extin1)
      call split_path(posargs(2)%arg, dummy, dummy, extin2)
      call open2read(posargs(1)%arg, unitin1)
      call open2read(posargs(2)%arg, unitin2)
   case default
      write (stderr, '(A)') 'Error: Too many file paths'
      stop 1
   end select
end if

if (fileout_flag) then
   call split_path(fileout_path, dummy, dummy, extout)
   call open2write(fileout_path, unitout)
end if

! Read coordinates
call readfile(unitin1, extin1, title1, atoms1, bonds1)
call readfile(unitin2, extin2, title2, atoms2, bonds2)

if (heavy_flag) then
   ! Include only heavy atoms
   call include_heavy_atoms(atoms1, atomset1, atomset1_alloc)
   call include_heavy_atoms(atoms2, atomset2, atomset2_alloc)
else
   ! Include all atoms
   call include_all_atoms(atoms1, atomset1, atomset1_alloc)
   call include_all_atoms(atoms2, atomset2, atomset2_alloc)
end if

! Collect atom types in a partition
call collect_atomtypes(atomset1, atomset2, atoms1, atoms2, atomtypes)

! Abort if atom types do not match
if (any(atomtypes%parts%num_items1 /= atomtypes%parts%num_items2)) then
   write (stderr, '(A)') 'Error: These molecules are not isomers'
   stop 1
end if

! Get adjacency information
if (bond_flag) then
   call adjacency_from_distance(atomset1, atoms1, adjcs1)
   call adjacency_from_distance(atomset2, atoms2, adjcs2)
else
   if (size(bonds1) < 1 .or. size(bonds2) < 1) then
      if (size(bonds1) < 1 .and. size(bonds2) < 1) then
         write (stdout,'(A)') 'ERROR: Molecules have no bonds!'
         stop
      else if (size(bonds1) < 1) then
         write (stdout,'(A)') 'ERROR: First molecule has no bonds!'
         stop
      else if (size(bonds2) < 1) then
         write (stdout,'(A)') 'ERROR: Second molecule has no bonds!'
         stop
      end if
   end if
   call adjacency_from_bonds(atomset1, atoms1, bonds1, adjcs1)
   call adjacency_from_bonds(atomset2, atoms2, bonds2, adjcs2)
end if

! Get user defined atom weights
if (mass_flag) then
   weights1 = atomic_masses(atoms1%elnum)
   weights2 = atomic_masses(atoms2%elnum)
else
   weights1 = uniform_weights(1._rk, size(atoms1))
   weights2 = uniform_weights(1._rk, size(atoms2))
end if

! Get mol1 coordinates
coords1 = get_coords(atoms1)

! Get mol2 coordinates
if (mirror_flag) then
   coords2 = get_mirrored_coords(atoms2)
else
   coords2 = get_coords(atoms2)
end if

if (align_flag) then

   center1 = get_centroid(atomset1, atoms1, weights1)
   center2 = get_centroid(atomset2, atoms2, weights2)
   call translate_coords(coords2, center1 - center2)

   ! Get weighted-centered coordinates
   coords1w = get_weighted_coords(atoms1, weights1, center1)
   coords2w = get_weighted_coords(atoms2, weights2, center2)

   if (remap_flag) then

      ! Allocate registries
      call allocate_registry(iso_registry, num_records)
      call allocate_registry(confo_registry, num_records)
      call allocate_registry(temp_registry, 1)

      ! Remap atoms to minimize adjacency difference and MSD
      call optimize_atomperm_isomer(atomset1, atomset2, atomtypes, adjcs1, adjcs2, &
                                     coords1w, coords2w, iso_registry)

      ! Print isomer optimization stats
      if (stats_flag) then
         call print_records(iso_registry)
      end if

      do i = 1, iso_registry%occ_records

         atomperm1 = iso_registry%records(i)%atomperm1

         ! Modify bonds according to selected strategy
         if (iso_registry%records(i)%permdiff == 0) then
            ! Already conformers, use original adjacencies
            call optimize_atomperm_conformer(atomset1, atomset2, adjcs1, adjcs2, atomtypes, &
                                          coords1w, coords2w, temp_registry)
         else
            if (union_flag) then
               ! Add all differing bonds to both molecules
               adjmat1 = adjcs_to_adjmat(adjcs1)
               adjmat2 = adjcs_to_adjmat(adjcs2)
               call bonds_union(adjmat1, adjmat2, atomperm1, iso_registry%records(i)%moldiffs)
               call adjmat_to_adjcs(adjmat1, adjcs1_mod)
               call adjmat_to_adjcs(adjmat2, adjcs2_mod)
               call optimize_atomperm_conformer(atomset1, atomset2, adjcs1_mod, adjcs2_mod, atomtypes, &
                                             coords1w, coords2w, temp_registry)
            else if (intersection_flag) then
               ! Remove all differing bonds from both molecules
               adjmat1 = adjcs_to_adjmat(adjcs1)
               adjmat2 = adjcs_to_adjmat(adjcs2)
               call bonds_intersection(adjmat1, adjmat2, atomperm1, iso_registry%records(i)%moldiffs)
               call adjmat_to_adjcs(adjmat1, adjcs1_mod)
               call adjmat_to_adjcs(adjmat2, adjcs2_mod)
               call optimize_atomperm_conformer(atomset1, atomset2, adjcs1_mod, adjcs2_mod, atomtypes, &
                                             coords1w, coords2w, temp_registry)
            else
               ! Default: match mol2 to mol1
               adjmat2 = adjcs_to_adjmat(adjcs2)
               call match_bonds_to_mol1(adjmat2, iso_registry%records(i)%moldiffs)
               call adjmat_to_adjcs(adjmat2, adjcs2_mod)
               call optimize_atomperm_conformer(atomset1, atomset2, adjcs1, adjcs2_mod, atomtypes, &
                                             coords1w, coords2w, temp_registry)
            end if
         end if

         call insert_record_atomperm(confo_registry, &
               temp_registry%records(1)%atomperm1, &
               temp_registry%records(1)%steps, &
               temp_registry%records(1)%rotation, &
               adjacencydiff(atomset1, temp_registry%records(1)%atomperm1, adjcs1, adjcs2), &
               temp_registry%records(1)%permdist)

      end do

      ! Print conformer optimization stats
      if (stats_flag) then
         call print_records(confo_registry)
      end if

      ! Process conformer optimization results
      do i = 1, confo_registry%occ_records
         atomperm1 = confo_registry%records(i)%atomperm1
         rotquat = least_rotquat(atomset1, atomperm1, coords1w, coords2w)
         coords2r = rotated_coords(coords2, rotquat, center1)
         rmsd = sqrt(sqdistmean(atomset1, atomperm1, weights1, coords1, coords2r))
         adjd = adjacencydiff(atomset1, atomperm1, adjcs1, adjcs2)

         if (coords_flag) then
            title2 = 'RMSD=' // str(rmsd) // ' Δadj=' // str(adjd)
            call set_coords(atoms2, coords2r)
            call writefile(unitout, extout, title2, atoms2, bonds2, atomperm1)
         else
            write (stdout,'(A,A,I0,A)') str(rmsd), '(', adjd, ')'
         end if
      end do

   else

      allocate (atomperm1(size( coords1w, 2)))
      call init_identity_permutation( atomperm1)
      rotquat = least_rotquat(atomset1, atomperm1, coords1w, coords2w)
      coords2r = rotated_coords(coords2, rotquat, center1)
      rmsd = sqrt(sqdistmean(atomset1, atomperm1, weights1, coords1, coords2r))
      adjd = adjacencydiff(atomset1, atomperm1, adjcs1, adjcs2)

      if (coords_flag) then
         title2 = 'RMSD=' // str(rmsd) // ' Δadj=' // str(adjd)
         call set_coords(atoms2, coords2r)
         call writefile(unitout, extout, title2, atoms2, bonds2, atomperm1)
      else
         write (stdout,'(A,A,I0,A)') str(rmsd), '(', adjd, ')'
      end if

   end if

else

   ! Get weighted coordinates
   coords1w = get_weighted_coords(atoms1, weights1)
   coords2w = get_weighted_coords(atoms2, weights2)

   if (remap_flag) then
      write (stdout, '(A)') 'ERROR: Remapping without alignment is not implemented'
      stop
   else
      allocate (atomperm1(size( coords1w, 2)))
      call init_identity_permutation( atomperm1)
      rmsd = sqrt(sqdistmean(atomset1, atomperm1, weights1, coords1, coords2))
      adjd = adjacencydiff(atomset1, atomperm1, adjcs1, adjcs2)
   end if

   if (coords_flag) then
      title2 = 'RMSD=' // str(rmsd) // ' Δadj=' // str(adjd)
      call set_coords(atoms2, coords2)
      call writefile(unitout, extout, title2, atoms2, bonds2, atomperm1)
   else
      write (stdout,'(A,A,I0,A)') str(rmsd), '(', adjd, ')'
   end if

end if

end program
!> @}
