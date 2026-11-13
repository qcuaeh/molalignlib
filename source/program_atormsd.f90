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
use registry
use alignment_atoms
use file_utils
use file_reading
use file_writing
use arg_parsing
use flags
implicit none

logical(LK) :: printmapping_flag
logical(LK) :: coords_flag, savecoords_flag
character(:), allocatable :: title1, title2
character(:), allocatable :: arg, pairs_prefix
character(:), allocatable :: in_format, out_format
character(LINE_LEN) :: posargs(2)
type(atom_t), dimension(:), allocatable :: atoms1, atoms2
type(bond_t), dimension(:), allocatable :: bonds1, bonds2
type(bool_matrix), dimension(:), allocatable :: prunes
type(partition_t) :: atomtypes
type(register_t) :: register
real(RK) :: rmsd
real(RK) :: center1(3), center2(3), rotquat(4)
real(RK), dimension(:), allocatable :: weights1, weights2, unit_weights
real(RK) :: transmat2(3,3)
real(RK), dimension(:,:), allocatable :: coords1, coords2, coords1w, coords2w, coords2r
real(RK), dimension(:,:), allocatable :: full_coords1, full_coords2, full_coords2r
integer(IK), dimension(:), allocatable :: atomset1, atomset2
integer(IK), dimension(:), allocatable :: mapping1, full_mapping1, ref_mapping
integer(IK) :: max_recs, max_trials, min_freq
integer(IK) :: in_unit
integer(IK) :: error_code
integer(IK) :: n_atoms1, n_atoms2, n_padding
integer(IK) :: i

! Set default options

heavy_flag = .FALSE.
mirror_flag = .FALSE.
align_flag = .FALSE.
remap_flag = .FALSE.
massweight_flag = .FALSE.
useatomtype_flag = .FALSE.
random_flag = .FALSE.
printstats_flag = .FALSE.
printmapping_flag = .FALSE.
coords_flag = .FALSE.
savecoords_flag = .FALSE.

max_recs = 1
min_freq = ATO_FREQ_DEFAULT
max_trials = MAX_TRIALS_DEFAULT
prune_procedure => prune_none

! Print the usage when no arguments are given

if (command_argument_count() == 0) then
   call print_usage()
   stop
end if

! Read command line options

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
   case ('-useatomtype')
      useatomtype_flag = .TRUE.
   case ('-heavy')
      heavy_flag = .TRUE.
   case ('-massweight')
      massweight_flag = .TRUE.
   case ('-mirror')
      mirror_flag = .TRUE.
   case ('-minfreq')
      call read_optarg(arg, min_freq, 1)
   case ('-maxtrials')
      call read_optarg( arg, max_trials, 1)
   case ('-maxrecs')
      call read_optarg( arg, max_recs, 1)
   case ('-coords')
      coords_flag = .TRUE.
      call read_optarg( arg, out_format)
      call check_out_format( arg, out_format)
   case ('-savecoords')
      savecoords_flag = .TRUE.
      call read_optarg( arg, pairs_prefix)
   case ('-printmapping')
      printmapping_flag = .TRUE.
   case ('-printstats')
      printstats_flag = .TRUE.
   case ('-random')
      random_flag = .TRUE.
   case default
      call read_posarg( arg, posargs)
   end select
end do

! -savecoords takes the output format from -coords
if (savecoords_flag .and. .not. coords_flag) then
   write (stderr, '(A)') 'Option -savecoords requires -coords'
   stop
end if

select case (ipos)
case (0)
   write (stderr, '(A)') 'File paths are missing'
   stop
case (1)
   write (stderr, '(A)') 'Too few file paths'
   stop
case (2)
   call open2read( trim(posargs(1)), in_format, in_unit)
   call read_file( in_unit, in_format, title1, atoms1, bonds1)
   close (in_unit)
   call open2read( trim(posargs(2)), in_format, in_unit)
   call read_file( in_unit, in_format, title2, atoms2, bonds2)
   close (in_unit)
case default
   write (stderr, '(A)') 'Too many file paths'
   stop
end select

! Pad the smaller cluster with dummy atoms appended at the end. Real atoms
! keep their indices; padding atoms are recognised by index from here on.
n_atoms1 = size(atoms1)
n_atoms2 = size(atoms2)
n_padding = max(n_atoms1, n_atoms2)
call pad_atoms( atoms1, n_padding)
call pad_atoms( atoms2, n_padding)

! Molecule 1 is written in input order (with -savecoords)
call init_array( ref_mapping, n_padding, identity)

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
! complete_mapping turns it into full_mapping1 over all (padded) atoms.

! Collect atom types in a partition
call collect_atomtypes( atoms1(atomset1), atoms2(atomset2), atomtypes)

! Abort if molecules are not isomers
if (any(atomtypes%parts%n_items1 /= atomtypes%parts%n_items2)) then
   write (stderr, '(A)') 'These molecules are not isomers'
   stop
end if

! Weights of the included atoms, normalised to sum 1 over them so both
! molecules are scaled by the same factor
if (massweight_flag) then
   weights1 = atomic_masses(atoms1(atomset1)%elnum)
   weights2 = atomic_masses(atoms2(atomset2)%elnum)
else
   allocate (weights1(size(atomset1)), weights2(size(atomset2)))
   weights1 = 1.0_RK
   weights2 = 1.0_RK
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
allocate (unit_weights(n_padding))
unit_weights = 1.0_RK
full_coords1 = get_coords( atoms1, unit_weights, ORIGIN, IDENTITY_MATRIX)
full_coords2 = get_coords( atoms2, unit_weights, ORIGIN, transmat2)

! Coordinates of the included atoms, in the compact numbering of mapping1
coords1 = full_coords1(:, atomset1)
coords2 = full_coords2(:, atomset2)

! Centers of the included atoms. They are taken from the coordinates, so
! center2 is in the transformed (e.g. mirrored) frame of molecule 2.
if (align_flag) then
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

   call prune_procedure( atomtypes, coords1, coords2, prunes)

   if (align_flag) then

      call allocate_register( register, max_recs)
      call optimize_mapping_atoms( atomtypes, prunes, &
            coords1w, coords2w, min_freq, max_trials, register, error_code)
      if (error_code /= 0) then
         write (stderr, '(A)') 'Error: Assignment failed'
         stop
      end if

      if (printstats_flag) call print_records( register)

      do i = 1, register%n_records
         mapping1 = register%records(i)%mapping1
         rotquat = least_rotquat( mapping1, coords1w, coords2w)
         full_coords2r = rotated_coords( full_coords2, rotquat, center1)
         coords2r = rotated_coords( coords2, rotquat, center1)
         rmsd = sqrt( sqdistmean( mapping1, weights1, coords1, coords2r))

         ! Pair the excluded atoms in the aligned frame (bonds, if the input
         ! format had any, are used as in conformsd)
         call complete_mapping( atomset1, atomset2, mapping1, atoms1, atoms2, &
               n_atoms1, n_atoms2, bonds1, bonds2, full_coords1, full_coords2r, full_mapping1)

         call report_result( i)
      end do

   else

      allocate (mapping1(size(atomset1)))
      call assign_atoms_pruned( atomtypes, coords1w, coords2w, prunes, mapping1, error_code)
      if (error_code /= 0) then
         write (stderr, '(A)') 'Error: Assignment failed'
         stop
      end if

      full_coords2r = full_coords2
      coords2r = coords2
      rmsd = sqrt( sqdistmean( mapping1, weights1, coords1, coords2r))

      call complete_mapping( atomset1, atomset2, mapping1, atoms1, atoms2, &
            n_atoms1, n_atoms2, bonds1, bonds2, full_coords1, full_coords2r, full_mapping1)

      call report_result( 1)

   end if

else

   ! Abort if atoms do not match: the same atoms must be included in both
   ! clusters, with the same types (the isomer check above guarantees that
   ! both atom sets have the same size)
   if (any(atomset1 /= atomset2)) then
      write (stderr, '(A)') 'Atoms do not match'
      stop
   end if
   if (any(atomtypes%itemdir1 /= atomtypes%itemdir2)) then
      write (stderr, '(A)') 'Atoms do not match'
      stop
   end if

   ! Maintain input atom order, both among the included atoms and in the
   ! full (padded) permutation
   call init_array( mapping1, size(atomset1), identity)
   call init_array( full_mapping1, size(atoms1), identity)

   if (align_flag) then
      rotquat = least_rotquat( mapping1, coords1w, coords2w)
      full_coords2r = rotated_coords( full_coords2, rotquat, center1)
      coords2r = rotated_coords( coords2, rotquat, center1)
   else
      full_coords2r = full_coords2
      coords2r = coords2
   end if

   rmsd = sqrt( sqdistmean( mapping1, weights1, coords1, coords2r))

   call report_result( 1)

end if

contains

subroutine print_usage()
   write (stderr, '(A)') 'Usage: atormsd file1 file2 [options]'
   write (stderr, '(A)') ''
   write (stderr, '(A)') 'Calculate the RMSD between two atom clusters.'
   write (stderr, '(A)') ''
   write (stderr, '(A)') 'Files are read by extension (.xyz, .mol, .sdf or .mol2); a path that is'
   write (stderr, '(A)') 'just a format name (e.g. xyz) reads standard input.'
   write (stderr, '(A)') ''
   write (stderr, '(A)') 'Options:'
   call print_option('-align', 'Align atoms to minimise the RMSD')
   call print_option('-remap', 'Remap atoms to minimise the RMSD')
   call print_option('-prune TOL', 'Prune atom pairs with tolerance TOL (Angstrom)')
   call print_option('-useatomtype', 'Only match atoms with the same label (e.g. C1, C2)')
   call print_option('-heavy', 'Ignore hydrogen atoms')
   call print_option('-massweight', 'Use mass-weighted coordinates')
   call print_option('-mirror', 'Reflect cluster 2 (x -> -x) before comparison')
   call print_option('-minfreq N', 'Stop once the best solution is found N times (default: ' &
         // str(ATO_FREQ_DEFAULT) // ')')
   call print_option('-maxtrials N', 'Stop after at most N random orientations (default: ' &
         // str(MAX_TRIALS_DEFAULT) // ')')
   call print_option('-maxrecs N', 'Record the N lowest RMSDs found (default: 1)')
   call print_option('-coords FORMAT', 'Write the coordinates of cluster 2 to stdout in FORMAT (xyz, sdf or mol2)')
   call print_option('-savecoords PREFIX', 'Write each cluster 1/cluster 2 pair to PREFIX_<n>.FORMAT (needs -coords)')
   call print_option('-printmapping', 'Print the 1-based atom permutation after each RMSD')
   call print_option('-printstats', 'Print optimisation statistics (needs -align -remap)')
   call print_option('-random', 'Seed the random-number generator from the system clock')
end subroutine

subroutine report_result(irec)
! Report result number irec. Without -coords, print the RMSD (and the atom mapping with -printmapping).
! With -coords, the coordinates of molecule 2 (as aligned, if at all, and
! reordered) are written in the -coords format with the RMSD in the title:
! to standard output instead of the RMSD line or, with -savecoords, preceded
! by the coordinates of molecule 1 in the file <prefix>_<irec>.<format>,
! and the RMSD line is printed as usual.
   integer(IK), intent(in) :: irec
   ! Local variables
   integer(IK) :: pair_unit

   if (coords_flag) call set_coords( atoms2, full_coords2r)

   if (savecoords_flag) then
      call open2write_indexed( pairs_prefix, irec, out_format, pair_unit)
      call write_file( pair_unit, out_format, title1, atoms1, bonds1, ref_mapping)
      call write_file( pair_unit, out_format, 'rmsd=' // str( rmsd), &
            atoms2, bonds2, full_mapping1)
      close (pair_unit)
   end if

   if (coords_flag .and. .not. savecoords_flag) then
      call write_file( int(stdout, IK), out_format, 'rmsd=' // str( rmsd), &
            atoms2, bonds2, full_mapping1)
   else
      write (stdout,'(A)',advance='no') str( rmsd)
      if (printmapping_flag) then
         write (stdout,'(1X)',advance='no')
         call print_permutation( full_mapping1)
      end if
      write (stdout,*)
   end if
end subroutine

subroutine print_option(head, desc)
! One line of the option list: the option and its argument, then the
! description.
   character(*), intent(in) :: head, desc
   character(18) :: column

   column = head
   write (stderr, '(2X,A,1X,A)') column, desc
end subroutine

end program
!> @}
