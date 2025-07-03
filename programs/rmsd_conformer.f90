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
use globals
use molecule
use spatial_transforms
use strutils
use chemutils
use adjacency
use permutation
use fileio
use argparse
use biasing
use pruning
use assignment_conformer
use molalignlib
use registration

implicit none

integer, allocatable :: atomperm(:)
character(:), allocatable :: arg
character(:), allocatable :: pathout
character(:), allocatable :: fmtin1, fmtin2, fmtout, fmtpipe
logical :: align_flag, remap_flag, write_flag, pipe_flag, nrec_flag
real(rk) :: center1(3)
integer :: adjd
real(rk) :: rmsd
type(strlist_type) :: posargs(2)
type(mol_type) :: mol1, mol2, auxmol
type(registry_t) :: registry
real(rk), dimension(:), allocatable :: weights1, weights2
real(rk), dimension(:,:), allocatable :: coords1, coords2
logical, dimension(:,:), allocatable :: adjmat1, adjmat2
integer :: unitin1, unitin2, unitout
integer :: i

! Set default options

iter_flag = .true.
test_flag = .false.
stats_flag = .false.
mirror_flag = .false.
align_flag = .false.
remap_flag = .false.
write_flag = .false.
pipe_flag = .false.
nrec_flag = .false.

max_records = 1
max_count = 10
max_trials = huge( max_trials)

atomic_weights => ones
pathout = 'aligned.xyz'

! Get user options

call init_args()

do while (get_arg(arg))
   select case (arg)
   case ('-align')
      align_flag = .true.
   case ('-remap')
      remap_flag = .true.
   case ('-mass')
      atomic_weights => atomic_masses
   case ('-mirror')
      mirror_flag = .true.
   case ('-count')
      call read_optarg( arg, max_count)
   case ('-trials')
      call read_optarg( arg, max_trials)
   case ('-N')
      nrec_flag = .true.
      call read_optarg( arg, max_records)
   case ('-out')
      write_flag = .true.
      call read_optarg( arg, pathout)
   case ('-pipe')
      pipe_flag = .true.
      call read_optarg( arg, fmtpipe)
   case ('-stats')
      stats_flag = .true.
   case ('-test')
      test_flag = .true.
   case default
      call read_posarg( arg, posargs)
   end select
end do

if (pipe_flag) then
   unitin1 = stdin
   unitin2 = stdin
   unitout = stdout
   fmtin1 = fmtpipe
   fmtin2 = fmtpipe
   fmtout = fmtpipe
else
   select case (ipos)
   case (0)
      write (stderr, '(a)') 'Error: Missing file paths'
      stop
   case (1)
      write (stderr, '(a)') 'Error: Too few file paths'
      stop
   case (2)
      call open2read( posargs(1)%arg, unitin1, fmtin1)
      call open2read( posargs(2)%arg, unitin2, fmtin2)
   case default
      write (stderr, '(a)') 'Error: Too many file paths'
      stop
   end select
   if (write_flag) then
      call open2write( pathout, unitout, fmtout)
   end if
end if

! Read coordinates
call read_file( unitin1, fmtin1, mol1)
call read_file( unitin2, fmtin2, mol2)

! Initialization
call set_bonds( mol1)
call set_bonds( mol2)
adjmat1 = get_adjmat( mol1)
coords1 = get_coords( mol1)
weights1 = atomic_weights(mol1%atoms%elnum)
weights2 = atomic_weights(mol2%atoms%elnum)
call weight_coords( coords1, weights1)
center1 = centroid( coords1)
allocate (auxmol%atoms(size(mol2%atoms)))

if (align_flag .and. remap_flag) then

   ! Remap atoms to minimize the MSD
   call optimize_atomperm_conformer( mol1, mol2, registry)

   ! Print optimization stats
   if (stats_flag) then
      call print_records( registry)
   end if

   if (write_flag .and. .not. nrec_flag) then
      call write_file( unitout, fmtout, mol1)
   end if

   do i = 1, registry%num_records
      atomperm = registry%records(i)%atomperm
      adjmat2 = registry%records(i)%adjmat2
      adjd = adjacencydiff( atomperm, adjmat1, adjmat2)
      coords2 = registry%records(i)%coords2
      rmsd = sqrt( total_sqdist( atomperm, coords1, coords2))
      call unweight_coords( coords2, weights2)

      write (stderr, "(a,',',a)") str( adjd), str( rmsd, 4)
      if (write_flag) then
         auxmol%title = 'RMSD=' // str( rmsd, 4)
         auxmol%atoms%elnum = mol2%atoms(atomperm)%elnum
         auxmol%atoms%label = mol2%atoms(atomperm)%label
         call set_coords( auxmol, coords2(:,atomperm))
         call write_file( unitout, fmtout, auxmol)
      end if
   end do

else if (align_flag) then

   adjmat2 = get_adjmat( mol2)
   adjd = adjacencydiff( atomperm, adjmat1, adjmat2)
   call align_atoms( mol1, mol2, coords2)
   call unweight_coords( coords2, weights2)
   call translate_coords( coords2, center1)
   rmsd = sqrt( total_sqdist( atomperm, coords1, coords2))

   write (stderr, "(a,',',a)") str( adjd), str( rmsd, 4)
   if (write_flag) then
      auxmol%title = 'RMSD=' // str( rmsd, 4)
      auxmol%atoms%elnum = mol2%atoms%elnum
      auxmol%atoms%label = mol2%atoms%label
      call set_coords( auxmol, coords2)
      call write_file( unitout, fmtout, mol2)
   end if

else if (remap_flag) then

   ! Not implemented

else

   ! Not implemented

end if

end program
!> @}
