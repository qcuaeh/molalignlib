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
use moltypes
use molalignlib

implicit none

integer :: i
integer :: nrec, error
integer :: read_unit0, read_unit1, write_unit
integer, allocatable, dimension(:) :: countlist
integer, allocatable, dimension(:, :) :: maplist
character(:), allocatable :: arg
character(:), allocatable :: pathin1, pathin2, pathout
character(:), allocatable :: fmtin0, fmtin1, fmtstdin, fmtout
character(ll) :: posargs(2)
real(wp) :: travec(3), rotmat(3, 3)
logical :: remap_flag, stdin_flag, stdout_flag
type(Molecule) :: mol0, mol1, auxmol
integer :: adjd, minadjd
real(wp) :: rmsd, minrmsd

procedure(f_realint), pointer :: weight_func

! Set default options

iter_flag = .false.
bias_flag = .false.
bond_flag = .false.
back_flag = .false.
test_flag = .false.
react_flag = .false.
stats_flag = .false.
trial_flag = .false.
stdin_flag = .false.
stdout_flag = .false.
mirror_flag = .false.
remap_flag = .false.

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
   case ('-sort')
      remap_flag = .true.
   case ('-back')
      back_flag = .true.
   case ('-fast')
      iter_flag = .true.
      bias_flag = .true.
   case ('-bond')
      bond_flag = .true.
!         react_flag = .true.
      print_stats => print_stats_diff
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
   case ('-stdin')
      stdin_flag = .true.
      call readoptarg(arg, fmtstdin)
   case ('-stdout')
      stdout_flag = .true.
      call readoptarg(arg, fmtout)
   case default
      call readposarg(arg, posargs)
   end select

end do

if (stdin_flag) then
   fmtin0 = fmtstdin
   fmtin1 = fmtstdin
   read_unit0 = input_unit
   read_unit1 = input_unit
else
   select case (ipos)
   case (0)
      write (error_unit, '(a)') 'Error: Missing arguments'
      stop
   case (1)
      write (error_unit, '(a)') 'Error: Too few arguments'
      stop
   case (2)
      pathin1 = trim(posargs(1))
      pathin2 = trim(posargs(2))
      call open2read(pathin1, read_unit0, fmtin0)
      call open2read(pathin2, read_unit1, fmtin1)
   case default
      write (error_unit, '(a)') 'Error: Too many arguments'
      stop
   end select
end if

! Read coordinates

call readfile(read_unit0, fmtin0, mol0)
call readfile(read_unit1, fmtin1, mol1)

! Allocate arrays

allocate(maplist(mol0%natom, maxrec))
allocate(countlist(maxrec))

! Get atomic numbers, types and weights

do i = 1, mol0%natom
   call readlabel(mol0%atoms(i)%label, mol0%atoms(i)%znum, mol0%atoms(i)%type)
   mol0%atoms(i)%weight = weight_func(mol0%atoms(i)%znum)
end do

do i = 1, mol1%natom
   call readlabel(mol1%atoms(i)%label, mol1%atoms(i)%znum, mol1%atoms(i)%type)
   mol1%atoms(i)%weight = weight_func(mol1%atoms(i)%znum)
end do

if (stdout_flag) then
   write_unit = output_unit
else
   fmtout = baseext(pathout)
   if (len(fmtout) > 0) then
      call open2write(pathout, write_unit)
   else
      write (error_unit, '(a)') 'Error: Output file must have an extension'
      stop
   end if
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

   call assign_atoms( &
      mol0%natom, &
      mol0%get_znums(), &
      mol0%get_types(), &
      mol0%get_weights(), &
      mol0%get_coords(), &
      mol0%adjmat, &
      mol1%natom, &
      mol1%get_znums(), &
      mol1%get_types(), &
      mol1%get_weights(), &
      mol1%get_coords(), &
      mol1%adjmat, &
      maplist, &
      countlist, &
      nrec, &
      error)

   if (error /= 0) stop

!      mol0%title = 'Reference'
!      call writefile(write_unit, fmtout, mol0)

   auxmol = mol1
   minrmsd = huge(rmsd)
   minadjd = huge(adjd)

   do i = 1, nrec

      mol1 = auxmol
      call mol1%permutate(maplist(:, i))

      call align_atoms( &
         mol0%natom, &
         mol0%get_znums(), &
         mol0%get_types(), &
         mol0%get_weights(), &
         mol0%get_coords(), &
         mol1%natom, &
         mol1%get_znums(), &
         mol1%get_types(), &
         mol1%get_weights(), &
         mol1%get_coords(), &
         travec, &
         rotmat, &
         error)

      call mol1%rotate(rotmat)
      call mol1%translate(travec)

      rmsd = get_rmsd(mol0, mol1)
      adjd = get_adjd(mol0, mol1)
      minrmsd = min(minrmsd, rmsd)
      minadjd = min(minadjd, adjd)
      mol1%title = 'Map='//istr(i)//' RMSD='//rstr(rmsd, 4)
      call writefile(write_unit, fmtout, mol1)

   end do

   if (.not. stats_flag) then
      write (output_unit, '(a,1x,a)') istr(minadjd), rstr(minrmsd, 4)
   end if

else

   ! Align atoms to minimize RMSD

   call align_atoms( &
      mol0%natom, &
      mol0%get_znums(), &
      mol0%get_types(), &
      mol0%get_weights(), &
      mol0%get_coords(), &
      mol1%natom, &
      mol1%get_znums(), &
      mol1%get_types(), &
      mol1%get_weights(), &
      mol1%get_coords(), &
      travec, &
      rotmat, &
      error)

   if (error /= 0) stop

   call mol1%rotate(rotmat)
   call mol1%translate(travec)

   rmsd = get_rmsd(mol0, mol1)
   adjd = get_adjd(mol0, mol1)

!      write (output_unit, '(a)') istr(adjd)
   write (output_unit, '(a)') rstr(rmsd, 4)

!      mol0%title = 'Reference'
!      call writefile(write_unit, fmtout, mol0)

   mol1%title = 'RMSD='//rstr(rmsd, 4)
   call writefile(write_unit, fmtout, mol1)

end if

end program
