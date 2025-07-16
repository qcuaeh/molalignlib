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

!> @defgroup program_rmsd_conformer Program RMSD Conformer
!> @brief Program to calculate RMSD between conformers
!> @{
program rmsd_conformer
use parameters
use options
use molecule
use spatial_transforms
use strutils
use chemistry
use adjacency
use permutation
use file_read
use file_write
use argparse
use biasing
use pruning
use registration
use assignment_cluster

implicit none

integer, allocatable :: atomperm(:)
character(:), allocatable :: title1, title2
character(:), allocatable :: arg, pathout, dummy
character(:), allocatable :: extin1, extin2, extout, extpipe
logical :: align_flag, remap_flag, write_flag, pipe_flag
type(strlist_type) :: posargs(2)
type(atom_t), dimension(:), allocatable :: atoms1, atoms2
type(partition_t) :: atomtypes
type(registry_t) :: registry
real(rk) :: rmsd
real(rk) :: center1(3), center2(3), rotquat(4)
real(rk), dimension(:), allocatable :: weights1, weights2
real(rk), dimension(:,:), allocatable :: coords1, coords2, wcoords1, wcoords2, rcoords2
integer :: unitin1, unitin2, unitout
integer :: i

! Set default options

adjacency_flag = .false.
iter_flag = .false.
test_flag = .false.
stats_flag = .false.
mirror_flag = .false.
align_flag = .false.
remap_flag = .false.
write_flag = .false.
pipe_flag = .false.

max_records = 1
max_count = 10
max_trials = huge( max_trials)

atomic_weights => ones
pathout = 'aligned.xyz'

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
!      iter_flag = .true.
      prune_procedure => prune_none
   case ('-prune')
      iter_flag = .true.
      prune_procedure => prune_rd
   case ('-tol')
      call read_optarg(arg, prune_tol)
   case ('-mass')
      atomic_weights => atomic_masses
   case ('-mirror')
      mirror_flag = .true.
   case ('-count')
      call read_optarg( arg, max_count)
   case ('-trials')
      call read_optarg( arg, max_trials)
   case ('-N')
      call read_optarg( arg, max_records)
   case ('-O')
      write_flag = .true.
      call read_optarg( arg, pathout)
   case ('-pipe')
      pipe_flag = .true.
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
      stop
   case (1)
      write (stderr, '(A)') 'Error: Too few file paths'
      stop
   case (2)
      call split_path( posargs(1)%arg, dummy, dummy, extin1)
      call split_path( posargs(2)%arg, dummy, dummy, extin2)
      call open2read( posargs(1)%arg, unitin1)
      call open2read( posargs(2)%arg, unitin2)
   case default
      write (stderr, '(A)') 'Error: Too many file paths'
      stop
   end select
   if (write_flag) then
      call open2write( pathout, unitout)
   end if
end if

! Read coordinates
call readmol( unitin1, extin1, title1, atoms1)
call readmol( unitin2, extin2, title2, atoms2)

! Abort if molecules have different number of atoms
if (size(atoms1) /= size(atoms2)) then
   write (stderr, '(A)') 'Error: These molecules are not isomers'
   stop
end if

! Abort if molecules are not isomers
if (any(sorted(atoms1%elnum) /= sorted(atoms2%elnum))) then
   write (stderr, '(A)') 'Error: These molecules are not isomers'
   stop
end if

! Compute atomic types
call collect_atomtypes( atoms1, atoms2, atomtypes)

! Abort if there are conflicting atomic types
if (any(atomtypes%parts%num_items1 /= atomtypes%parts%num_items2)) then
   write (stderr, '(A)') 'Error: There are conflicting atomic types'
   stop
end if

! Get standard coordinates
coords1 = get_coords( atoms1)
coords2 = get_coords( atoms2)
weights1 = atomic_weights(atoms1%elnum)
weights2 = atomic_weights(atoms2%elnum)
center1 = centroid( coords1, weights1)
center2 = centroid( coords2, weights2)
call translate_coords( coords2, center1 - center2)

! Get normalized coordinates
wcoords1 = get_coords( atoms1)
wcoords2 = get_coords( atoms2)
call translate_coords( wcoords1, -center1)
call translate_coords( wcoords2, -center2)
call weight_coords( wcoords1, weights1)
call weight_coords( wcoords2, weights2)

if (align_flag .and. remap_flag) then

   ! Remap atoms to minimize the MSD
   call optimize_atomperm_cluster( atoms1, atoms2, atomtypes, registry)

   ! Print optimization stats
   if (stats_flag) then
      call print_records( registry)
   end if

   do i = 1, registry%num_records
      atomperm = registry%records(i)%atomperm
!      rotquat = registry%records(i)%rotquat
      rotquat = least_rotquat( atomperm, wcoords1, wcoords2)
      rcoords2 = rotated_coords( coords2, rotquat, center1)
      rmsd = sqrt( total_sqdist( atomperm, weights1, coords1, rcoords2))

      write (stderr,'(A)') str( rmsd, 4)

      if (write_flag .or. pipe_flag) then
         title2 = 'RMSD=' // str( rmsd, 4)
         call set_coords( atoms2, rcoords2)
         call writemol( unitout, extout, title2, atoms2, atomperm)
      end if
   end do

else if (align_flag) then

   rotquat = least_rotquat( wcoords1, wcoords2)
   rcoords2 = rotated_coords( coords2, rotquat, center1)
   rmsd = sqrt( total_sqdist( weights1, coords1, rcoords2))

   write (stderr,'(A)') str( rmsd, 4)

   if (write_flag .or. pipe_flag) then
      title2 = 'RMSD=' // str( rmsd, 4)
      call set_coords( atoms2, rcoords2)
      call writemol( unitout, extout, title2, atoms2)
   end if

else if (remap_flag) then

   ! Not implemented

else

   ! Not implemented

end if

end program
!> @}
