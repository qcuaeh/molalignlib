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

!> @defgroup conformsd ConfoRMSD
!> @brief Program to calculate RMSDs between conformers
!> @{
program conformsd
use parameters
use str_utils
use chemdata
use molecule
use adjacency
use euclidean
use permutation
use assorting
use recording
use assignment_conformer
use alignment_conformer
use file_utils
use file_reading
use file_writing
use arg_parsing
use flags
implicit none

logical(lk) :: printmapping_flag
logical(lk) :: write_aligned
character(:), allocatable :: title1, title2
character(:), allocatable :: arg, aligned_path
character(:), allocatable :: in_format1, in_format2, out_format
type(strlist_type) :: posargs(2)
type(atom_t), dimension(:), allocatable :: atoms1, atoms2
type(bond_t), dimension(:), allocatable :: bonds1, bonds2
type(adjc_t), dimension(:), allocatable :: adjcs1, adjcs2
type(partition_t) :: atomtypes
type(registry_t) :: registry
real(rk) :: rmsd
real(rk) :: center1(3), center2(3), rotquat(4)
real(rk), dimension(:), allocatable :: weights1, weights2, unit_weights
real(rk) :: transmat2(3,3)
real(rk), dimension(:,:), allocatable :: coords1, coords2, coords1w, coords2w, coords2r
real(rk), dimension(:,:), allocatable :: full_coords1, full_coords2, full_coords2r
integer(ik), dimension(:), allocatable :: atomset1, atomset2
integer(ik), dimension(:), allocatable :: bondtypes
integer(ik), dimension(:), allocatable :: mapping1, full_atomperm1
integer(ik) :: max_records, max_trials, confo_freq, max_fragments
integer(ik) :: n_frags1, n_frags2
integer(ik) :: in_unit, aligned_unit
integer(ik) :: error_code
integer(ik) :: n_atoms1, n_atoms2, n_padding
integer(ik) :: i

! Set default options

heavy_flag = .FALSE.
mirror_flag = .FALSE.
align_flag = .FALSE.
remap_flag = .FALSE.
massweight_flag = .FALSE.
bonding_flag = .FALSE.
usebondtype_flag = .FALSE.
useatomtype_flag = .FALSE.
random_flag = .FALSE.
printstats_flag = .FALSE.
printassigntree_flag = .FALSE.
printmapping_flag = .FALSE.
write_aligned = .FALSE.

max_records = 1
confo_freq = CONFO_FREQ_DEFAULT
max_trials = MAX_TRIALS_DEFAULT
max_fragments = MAX_FRAGS_DEFAULT

! Read command line options

call init_args()

do while (get_arg(arg))
   select case (lowercase(arg))
   case ('-align')
      align_flag = .TRUE.
   case ('-remap')
      remap_flag = .TRUE.
   case ('-bondtol')
      bonding_flag = .TRUE.
      call read_optarg(arg, bond_tol)
   case ('-bondtype')
      usebondtype_flag = .TRUE.
   case ('-atomtype')
      useatomtype_flag = .TRUE.
   case ('-heavy')
      heavy_flag = .TRUE.
   case ('-massweight')
      massweight_flag = .TRUE.
   case ('-mirror')
      mirror_flag = .TRUE.
   case ('-confofreq')
      call read_optarg( arg, confo_freq, 1_ik)
   case ('-maxtrials')
      call read_optarg( arg, max_trials, 1_ik)
   case ('-maxfrags')
      call read_optarg( arg, max_fragments, 1_ik)
   case ('-maxrecs')
      call read_optarg( arg, max_records, 1_ik)
   case ('-aligned')
      write_aligned = .TRUE.
      call read_optarg( arg, aligned_path)
   case ('-mapping')
      printmapping_flag = .TRUE.
   case ('-assigntree')
      printassigntree_flag = .TRUE.
   case ('-stats')
      printstats_flag = .TRUE.
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
   call open2read( posargs(1)%arg, in_format1, in_unit)
   call read_file( in_unit, in_format1, title1, atoms1, bonds1)
   close (in_unit)
   call open2read( posargs(2)%arg, in_format2, in_unit)
   call read_file( in_unit, in_format2, title2, atoms2, bonds2)
   close (in_unit)
case default
   stop 'Too many file paths'
end select

! Bond types are compared, not interpreted, so they must come from the
! same parser: both files must have the same format
if (usebondtype_flag .and. .not. bonding_flag) then
   if (in_format1 /= in_format2) then
      stop 'Bond types can only be compared between files of the same format'
   end if
end if

! Pad the smaller molecule with dummy atoms appended at the end. Real atoms
! keep their indices; padding atoms are recognised by index from here on.
n_atoms1 = size(atoms1)
n_atoms2 = size(atoms2)
n_padding = max(n_atoms1, n_atoms2)
call pad_atoms( atoms1, n_padding)
call pad_atoms( atoms2, n_padding)

! Atom sets only ever contain real atoms
if (heavy_flag) then
   call include_heavy_atoms( atoms1(1:n_atoms1), atomset1)
   call include_heavy_atoms( atoms2(1:n_atoms2), atomset2)
else
   call include_all_atoms( atoms1(1:n_atoms1), atomset1)
   call include_all_atoms( atoms2(1:n_atoms2), atomset2)
end if

! From here on the comparison works on the included atoms only, numbered
! 1..size(atomset1) in molecule 1 and 1..size(atomset2) in molecule 2.
! mapping1 maps included atoms to included atoms in that compact numbering;
! complete_mapping turns it into full_atomperm1 over all (padded) atoms.

! Collect atom types in a partition
call collect_atomtypes( atoms1(atomset1), atoms2(atomset2), atomtypes)

! Abort if molecules are not isomers
if (any(atomtypes%parts%n_items1 /= atomtypes%parts%n_items2)) then
   stop 'These molecules are not isomers'
end if

! With -bondtol, perceive bonds from geometry instead of the files
if (bonding_flag) then
   call bonds_from_atoms( atoms1(1:n_atoms1), bonds1)
   call bonds_from_atoms( atoms2(1:n_atoms2), bonds2)
end if

! Set adjacency lists of the included atoms (bonds to excluded atoms are
! dropped). With -bondtype, the edges used by the refinement carry the bond
! types, compacted jointly for both molecules.
if (usebondtype_flag) then
   bondtypes = distinct_bondtypes( bonds1, bonds2)
else
   ! No bond types: every bond is GENERIC_BOND
   allocate (bondtypes(0))
end if
call adjacency_from_bonds( atoms1(atomset1), extract_bonds( atomset1, n_padding, bonds1), adjcs1, bondtypes)
call adjacency_from_bonds( atoms2(atomset2), extract_bonds( atomset2, n_padding, bonds2), adjcs2, bondtypes)

! Abort if either molecule has more fragments than allowed (counted on the
! bond graph of the included atoms, which is what the search uses)
n_frags1 = count_fragments( adjcs1)
n_frags2 = count_fragments( adjcs2)
if (n_frags1 > max_fragments) then
   write (stderr, '(A,1X,I0,1X,A,1X,I0,A)') 'First molecule has', n_frags1, &
         'fragments, more than the maximum of', max_fragments, ' (see -maxfrags)'
   stop
end if
if (n_frags2 > max_fragments) then
   write (stderr, '(A,1X,I0,1X,A,1X,I0,A)') 'Second molecule has', n_frags2, &
         'fragments, more than the maximum of', max_fragments, ' (see -maxfrags)'
   stop
end if

! Weights of the included atoms, normalised to sum 1 over them so both
! molecules are scaled by the same factor
if (massweight_flag) then
   weights1 = atomic_masses(atoms1(atomset1)%elnum)
   weights2 = atomic_masses(atoms2(atomset2)%elnum)
else
   allocate (weights1(size(atomset1)), source=1.0_rk)
   allocate (weights2(size(atomset2)), source=1.0_rk)
end if
weights1 = weights1/sum(weights1)
weights2 = weights2/sum(weights2)

! Linear transformation applied to molecule 2
if (mirror_flag) then
   transmat2 = MIRROR_MATRIX
else
   transmat2 = IDENTITY_MATRIX
end if

! Coordinates of all (padded) atoms, molecule 2 transformed. They are only
! used to align the whole molecule 2 and to pair the excluded atoms.
allocate (unit_weights(n_padding), source=1.0_rk)
full_coords1 = get_coords( atoms1, unit_weights, ORIGIN, IDENTITY_MATRIX)
full_coords2 = get_coords( atoms2, unit_weights, ORIGIN, transmat2)

! Coordinates of the included atoms, in the compact numbering of mapping1
coords1 = full_coords1(:, atomset1)
coords2 = full_coords2(:, atomset2)

! Centers of the included atoms. They are taken from the coordinates, so
! center2 is in the transformed (e.g. mirrored) frame of molecule 2.
if (align_flag) then

   if (write_aligned) then
      call open2write( aligned_path, out_format, aligned_unit)
   end if

   center1 = get_centroid( coords1, weights1)
   center2 = get_centroid( coords2, weights2)
else
   center1 = ORIGIN
   center2 = ORIGIN
end if

! Weighted coordinates of the included atoms.
! coords2 is already transformed, hence IDENTITY_MATRIX for both.
coords1w = get_coords( coords1, weights1, center1, IDENTITY_MATRIX)
coords2w = get_coords( coords2, weights2, center2, IDENTITY_MATRIX)

! Move molecule 2 onto the center of molecule 1
if (align_flag) then
   call translate_coords( full_coords2, center1 - center2)
   call translate_coords( coords2, center1 - center2)
end if

if (remap_flag) then

   if (align_flag) then

      call allocate_registry( registry, max_records)
      call optimize_mapping_conformer( adjcs1, adjcs2, atomtypes, &
            coords1w, coords2w, confo_freq, max_trials, registry, error_code)
      if (error_code /= 0) stop 'These molecules are not conformers'

      if (printstats_flag) call print_records( registry)

      do i = 1, registry%n_records
         mapping1 = registry%records(i)%mapping1
         rotquat = least_rotquat( mapping1, coords1w, coords2w)
         full_coords2r = rotated_coords( full_coords2, rotquat, center1)
         coords2r = rotated_coords( coords2, rotquat, center1)
         rmsd = sqrt( sqdistmean( mapping1, weights1, coords1, coords2r))

         ! Pair the excluded atoms in the aligned frame
         call complete_mapping( atomset1, atomset2, mapping1, atoms1, atoms2, &
               n_atoms1, n_atoms2, bonds1, bonds2, full_coords1, full_coords2r, full_atomperm1)

         if (write_aligned) then
            title2 = 'rmsd=' // str( rmsd)
            call set_coords( atoms2, full_coords2r)
            call write_file( aligned_unit, out_format, title2, atoms2, bonds2, full_atomperm1)
         else
            write (stdout,'(A)',advance='no') str( rmsd)
            if (printmapping_flag) then
               write (stdout,'(1X)',advance='no')
               call print_permutation(full_atomperm1)
            end if
            write (stdout,*)
         end if
      end do

   else

      call assign_mapping_conformer( adjcs1, adjcs2, atomtypes, coords1w, coords2w, mapping1, error_code)
      if (error_code /= 0) stop 'These molecules are not conformers'

      full_coords2r = full_coords2
      coords2r = coords2
      rmsd = sqrt( sqdistmean( mapping1, weights1, coords1, coords2r))

      call complete_mapping( atomset1, atomset2, mapping1, atoms1, atoms2, &
            n_atoms1, n_atoms2, bonds1, bonds2, full_coords1, full_coords2r, full_atomperm1)

      write (stdout,'(A)',advance='no') str( rmsd)
      if (printmapping_flag) then
         write (stdout,'(1X)',advance='no')
         call print_permutation(full_atomperm1)
      end if
      write (stdout,*)

   end if

else

   ! Abort if atoms do not match: the same atoms must be included in both
   ! molecules, with the same types (the isomer check above guarantees
   ! that both atom sets have the same size)
   if (any(atomset1 /= atomset2)) then
      stop 'Atoms do not match'
   end if
   if (any(atomtypes%itemdir1 /= atomtypes%itemdir2)) then
      stop 'Atoms do not match'
   end if

   ! Maintain input atom order, both among the included atoms and in the
   ! full (padded) permutation
   call init_array( mapping1, size(atomset1), identity)
   call init_array( full_atomperm1, size(atoms1), identity)

   ! Abort if bonds do not match
   if (adjacencydiff( mapping1, adjcs1, adjcs2) > 0) then
      stop 'Bonds do not match'
   end if

   if (align_flag) then
      rotquat = least_rotquat( mapping1, coords1w, coords2w)
      full_coords2r = rotated_coords( full_coords2, rotquat, center1)
      coords2r = rotated_coords( coords2, rotquat, center1)
   else
      full_coords2r = full_coords2
      coords2r = coords2
   end if

   rmsd = sqrt( sqdistmean( mapping1, weights1, coords1, coords2r))

   if (align_flag .and. write_aligned) then
      title2 = 'rmsd=' // str( rmsd)
      call set_coords( atoms2, full_coords2r)
      call write_file( aligned_unit, out_format, title2, atoms2, bonds2, full_atomperm1)
   else
      write (stdout,'(A)') str( rmsd)
   end if

end if

end program
!> @}
