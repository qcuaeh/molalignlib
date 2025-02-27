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

program molalign
use parameters
use globals
use molecule
use printing
use rotation
use rigid_body
use strutils
use chemutils
use alignment
use adjacency
use permutation
use fileio
use argparse
use molalignlib
use biasing
use pruning

implicit none

integer :: i
integer :: read_unit1, read_unit2, write_unit
integer, allocatable :: atomperm(:)
character(:), allocatable :: arg
character(:), allocatable :: fmtin1, fmtin2, fmtout
character(:), allocatable :: optfmtin, optfmtout
character(:), allocatable :: pathout
logical :: fmtin_flag, fmtout_flag
logical :: remap_flag, pipe_flag, nrec_flag
real(rk) :: travec1(3), travec2(3), rotquat(4)
integer :: adjd
real(rk) :: rmsd
type(strlist_type) :: posargs(2)
type(mol_type) :: mol1, mol2, auxmol
type(registry_type) :: results

integer :: num_atoms1
real(rk), dimension(:, :), allocatable :: coords1, coords2

! Set default options

iter_flag = .true.
test_flag = .false.
reac_flag = .false.
stats_flag = .false.
mirror_flag = .false.
remap_flag = .false.
pipe_flag = .false.
fmtin_flag = .false.
fmtout_flag = .false.
nrec_flag = .false.

max_records = 1
max_count = 10
max_trials = huge(max_trials)

atomic_weights => ones
pathout = 'aligned.xyz'

!prune_tol = 0.5
!prune_procedure => prune_none

! Get user options

call init_args()

do while (get_arg(arg))

   select case (arg)
   case ('-stats')
      stats_flag = .true.
   case ('-test')
      test_flag = .true.
   case ('-remap')
      remap_flag = .true.
!   case ('-near')
!      iter_flag = .false.
!      prune_procedure => prune_none
!   case ('-prune')
!      iter_flag = .true.
!      prune_procedure => prune_rd
   case ('-reac')
      reac_flag = .true.
   case ('-mass')
      atomic_weights => atomic_masses
   case ('-mirror')
      mirror_flag = .true.
   case ('-count')
      call read_optarg(arg, max_count)
   case ('-trials')
      call read_optarg(arg, max_trials)
   case ('-tol')
      call read_optarg(arg, prune_tol)
   case ('-N')
      nrec_flag = .true.
      call read_optarg(arg, max_records)
   case ('-out')
      call read_optarg(arg, pathout)
   case ('-fmtin')
      fmtin_flag = .true.
      call read_optarg(arg, optfmtin)
   case ('-fmtout')
      fmtout_flag = .true.
      call read_optarg(arg, optfmtout)
   case ('-pipe')
      pipe_flag = .true.
   case default
      call read_posarg(arg, posargs)
   end select

end do

if (pipe_flag) then
   read_unit1 = stdin
   read_unit2 = stdin
   fmtin1 = 'xyz'
   fmtin2 = 'xyz'
else
   select case (ipos)
   case (0)
      write (stderr, '(a)') 'Error: Missing file paths'
      stop
   case (1)
      write (stderr, '(a)') 'Error: Too few file paths'
      stop
   case (2)
      call open2read(posargs(1)%arg, read_unit1, fmtin1)
      call open2read(posargs(2)%arg, read_unit2, fmtin2)
   case default
      write (stderr, '(a)') 'Error: Too many file paths'
      stop
   end select
end if

if (fmtin_flag) then
   fmtin1 = optfmtin
   fmtin2 = optfmtin
end if

! Read coordinates
call readfile( read_unit1, fmtin1, mol1)
call readfile( read_unit2, fmtin2, mol2)

! Allocate arrays
if (pipe_flag) then
   write_unit = stdout
   fmtout = 'xyz'
else
   call open2write( pathout, write_unit, fmtout)
end if

if (fmtout_flag) then
   fmtout = optfmtout
end if

if (remap_flag) then

   ! Remap atoms to minimize the MSD
   call molecule_remap( &
      mol1, &
      mol2, &
      results)

   ! Print optimization stats
   if (stats_flag) then
      call print_stats_adjd( results)
!      call print_stats_rmsd( results)
      call print_final_stats( results)
   end if

   if (.not. nrec_flag) then
!      call writefile( write_unit, fmtout, mol1)
   end if

   num_atoms1 = size(mol1%atoms)
   coords1 = mol1%get_weighted_coords()
   coords2 = mol2%get_weighted_coords()

   ! Calculate centroids
   travec1 = -centroid( coords1)
   travec2 = -centroid( coords2)

   allocate (auxmol%atoms(size(mol2%atoms)))
   allocate (auxmol%adjmat(size(mol2%adjmat, 1), size(mol2%adjmat, 2)))

   do i = 1, results%num_records

      atomperm = results%records(i)%atomperm
      coords2 = mol2%get_weighted_coords()

      ! Calculate optimal rotation matrix
      rotquat = leasteigquat( &
         atomperm, &
         translated_coords( coords1, travec1), &
         translated_coords( coords2, travec2) &
      )

      coords2 = mol2%get_weighted_coords()
      call translate_coords( coords2, travec2)
      call rotate_coords( coords2, rotquat)
      call translate_coords( coords2, -travec1)

      adjd = adjacencydiff( atomperm, mol1%adjmat, mol2%adjmat)
      rmsd = sqrt(sqdistsum( atomperm, coords1, coords2))

      write (stderr, "(a,',',a)") str(adjd), str(rmsd, 4)
!      write (stderr, "(a)") str(rmsd, 4)

      auxmol%title = 'RMSD='//str(rmsd, 4)
      auxmol%atoms%elnum = mol2%atoms(atomperm)%elnum
      auxmol%atoms%label = mol2%atoms(atomperm)%label
      auxmol%atoms%weight = mol2%atoms(atomperm)%weight
      auxmol%adjmat = mol2%adjmat(atomperm, atomperm)
      call auxmol%set_unweighted_coords(coords2(:, atomperm))
!      call writefile( write_unit, fmtout, auxmol)

   end do

else

   ! Align atoms to minimize RMSD
   call molecule_align( &
      mol1, &
      mol2, &
      travec1, &
      travec2, &
      rotquat)

   coords2 = mol2%get_weighted_coords()
   call translate_coords( coords2, travec2)
   call rotate_coords( coords2, rotquat)
   call translate_coords( coords2, -travec1)

   adjd = adjacencydiff( identity_perm(num_atoms1), mol1%adjmat, mol2%adjmat)
   rmsd = sqrt( sqdistsum( identity_perm(num_atoms1), coords1, coords2))

   write (stderr, "(a,',',a)") str(adjd), str(rmsd, 4)
!   write (stderr, "(a)") str(rmsd, 4)

   mol2%title = 'RMSD='//str(rmsd, 4)
   auxmol%atoms%elnum = mol2%atoms%elnum
   auxmol%atoms%label = mol2%atoms%label
   auxmol%atoms%weight = mol2%atoms%weight
   auxmol%adjmat = mol2%adjmat
   call auxmol%set_unweighted_coords(coords2)
   call writefile( write_unit, fmtout, mol2)

end if

end program
