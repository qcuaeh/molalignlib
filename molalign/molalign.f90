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
use kinds
use types
use flags
use bounds
use printing
use rotation
use translation
use strutils
use chemutils
use alignment
use adjacency
use discrete
use fileio
use argparse
use molalignlib

implicit none

integer :: i
integer :: nrec, error
integer :: read_unit0, read_unit1, write_unit
integer, allocatable, dimension(:) :: countlist
integer, allocatable, dimension(:, :) :: maplist
character(:), allocatable :: arg
character(:), allocatable :: pathin1, pathin2, pathout
character(:), allocatable :: fmtin0, fmtin1, fmtout, optfmtin, optfmtout
character(ll) :: posargs(2)
real(wp) :: travec(3), rotmat(3, 3)
logical :: remap_flag
logical :: stdin_flag, stdout_flag
logical :: fmtin_flag, fmtout_flag
type(Molecule) :: mol0, mol1, auxmol0, auxmol1
integer :: adjd, minadjd
real(wp) :: rmsd, minrmsd

procedure(f_realint), pointer :: weight_func

! Set default options

iter_flag = .false.
bias_flag = .false.
bond_flag = .false.
back_flag = .false.
test_flag = .false.
reac_flag = .false.
stats_flag = .false.
trial_flag = .false.
mirror_flag = .false.
remap_flag = .false.
stdin_flag = .false.
stdout_flag = .false.
fmtin_flag = .false.
fmtout_flag = .false.

maxrec = 1
maxcount = 10
maxcoord = 16
maxlevel = 16

bias_tol = 0.35
bias_scale = 1.e3

pathout = 'aligned.xyz'

weight_func => unity
print_stats => print_stats_dist

! Get user options

call initarg()

do while (getarg(arg))

   select case (arg)
   case ('-stats')
      stats_flag = .true.
   case ('-test')
      test_flag = .true.
   case ('-remap')
      remap_flag = .true.
   case ('-back')
      back_flag = .true.
   case ('-fast')
      iter_flag = .true.
      bias_flag = .true.
   case ('-bond')
      bond_flag = .true.
      print_stats => print_stats_diff
   case ('-reac')
      reac_flag = .true.
   case ('-mass')
      weight_func => stdmass
   case ('-mirror')
      mirror_flag = .true.
   case ('-count')
      call readoptarg(arg, maxcount)
   case ('-trials')
      trial_flag = .true.
      call readoptarg(arg, maxtrials)
   case ('-tol')
      call readoptarg(arg, bias_tol)
!      case ('-scale')
!         call readoptarg(arg, bias_scale)
   case ('-rec')
      call readoptarg(arg, maxrec)
   case ('-out')
      call readoptarg(arg, pathout)
   case ('-fmtin')
      fmtin_flag = .true.
      call readoptarg(arg, optfmtin)
   case ('-fmtout')
      fmtout_flag = .true.
      call readoptarg(arg, optfmtout)
   case ('-stdin')
      stdin_flag = .true.
   case ('-stdout')
      stdout_flag = .true.
   case default
      call readposarg(arg, posargs)
   end select

end do

if (stdin_flag) then
   read_unit0 = stdin
   read_unit1 = stdin
   fmtin0 = 'xyz'
   fmtin1 = 'xyz'
else
   select case (ipos)
   case (0)
      write (stderr, '(a)') 'Error: Missing arguments'
      stop
   case (1)
      write (stderr, '(a)') 'Error: Too few arguments'
      stop
   case (2)
      pathin1 = trim(posargs(1))
      pathin2 = trim(posargs(2))
      call open2read(pathin1, read_unit0, fmtin0)
      call open2read(pathin2, read_unit1, fmtin1)
   case default
      write (stderr, '(a)') 'Error: Too many arguments'
      stop
   end select
end if

if (fmtin_flag) then
   fmtin0 = optfmtin
   fmtin1 = optfmtin
end if

! Read coordinates

call readfile(read_unit0, fmtin0, mol0)
call readfile(read_unit1, fmtin1, mol1)

! Allocate arrays

allocate(maplist(mol0%natom, maxrec))
allocate(countlist(maxrec))

! Get atomic numbers, types and weights

do i = 1, mol0%natom
   call readlabel(mol0%atoms(i)%label, mol0%atoms(i)%znum, mol0%atoms(i)%ztype)
   mol0%atoms(i)%weight = weight_func(mol0%atoms(i)%znum)
end do

do i = 1, mol1%natom
   call readlabel(mol1%atoms(i)%label, mol1%atoms(i)%znum, mol1%atoms(i)%ztype)
   mol1%atoms(i)%weight = weight_func(mol1%atoms(i)%znum)
end do

if (stdout_flag) then
   write_unit = stdout
   fmtout = 'xyz'
else
   call open2write(pathout, write_unit, fmtout)
end if

if (fmtout_flag) then
   fmtout = optfmtout
end if

! Get adjacency matrices and lists

if (bond_flag) then
   call getadjmat(mol0%natom, mol0%get_coords(), mol0%get_znums(), mol0%adjmat)
   call getadjmat(mol1%natom, mol1%get_coords(), mol1%get_znums(), mol1%adjmat)
else
   mol0%adjmat = .false.
   mol1%adjmat = .false.
end if

! Sort atoms to minimize MSD

if (remap_flag) then

   auxmol0 = mol0
   auxmol1 = mol1

   call remap_atoms( &
      auxmol0, &
      auxmol1, &
      maplist, &
      countlist, &
      nrec, &
      error)

   if (error /= 0) stop

!      mol0%title = 'Reference'
!      call writefile(write_unit, fmtout, mol0)

   auxmol1 = mol1
   minrmsd = huge(rmsd)
   minadjd = huge(adjd)

   do i = 1, nrec

      mol1 = auxmol1
      call mol1%permutate_atoms(maplist(:, i))

      call align_atoms( &
         mol0, &
         mol1, &
         travec, &
         rotmat, &
         error)

      call mol1%rotate_coords(rotmat)
      call mol1%translate_coords(travec)

      rmsd = get_rmsd(mol0, mol1)
      adjd = get_adjd(mol0, mol1)
      minrmsd = min(minrmsd, rmsd)
      minadjd = min(minadjd, adjd)
      mol1%title = 'Map='//istr(i)//' RMSD='//rstr(rmsd, 4)
      call writefile(write_unit, fmtout, mol1)

   end do

   if (.not. stats_flag) then
      if (bond_flag) then
         write (stdout, '(a,1x,a)') istr(minadjd), rstr(minrmsd, 4)
      else
         write (stdout, '(a,1x,a)') rstr(minrmsd, 4)
      end if
   end if

else

   ! Align atoms to minimize RMSD

   call align_atoms( &
      mol0, &
      mol1, &
      travec, &
      rotmat, &
      error)

   if (error /= 0) stop

   call mol1%rotate_coords(rotmat)
   call mol1%translate_coords(travec)

   rmsd = get_rmsd(mol0, mol1)
   adjd = get_adjd(mol0, mol1)

!      write (stdout, '(a)') istr(adjd)
   write (stdout, '(a)') rstr(rmsd, 4)

!      mol0%title = 'Reference'
!      call writefile(write_unit, fmtout, mol0)

   mol1%title = 'RMSD='//rstr(rmsd, 4)
   call writefile(write_unit, fmtout, mol1)

end if

end program
