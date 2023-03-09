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
   use library
   use fileio
   use argparse

   implicit none

   integer :: i
   integer :: nrec, error
   integer :: natom0, natom1
   integer :: read_unit0, read_unit1, write_unit
   integer, allocatable, dimension(:) :: znums0, znums1
   integer, allocatable, dimension(:) :: types0, types1
   integer, allocatable, dimension(:) :: countlist
   integer, allocatable, dimension(:, :) :: permlist
   character(:), allocatable :: arg
   character(:), allocatable :: pathin1, pathin2, pathout
   character(:), allocatable :: fmtin, fmtin0, fmtin1, fmtout
   character(:), allocatable :: title0, title1
   character(ll) :: posargs(2)
   character(wl), allocatable, dimension(:) :: labels0, labels1
   integer :: adjd
   real(wp) :: rmsd, travec(3), rotmat(3, 3)
   real(wp), allocatable, dimension(:) :: weights0, weights1
   real(wp), allocatable, dimension(:, :) :: coords0, coords1, aligned1
   logical :: sort_flag, mirror_flag, stdin_flag, stdout_flag
   logical, dimension(:, :), allocatable :: adjmat0, adjmat1

   procedure(f_realint), pointer :: weight_func

   ! Set default options

   sort_flag = .false.
   iter_flag = .false.
   bias_flag = .false.
   bond_flag = .false.
   trial_flag = .false.
   stdin_flag = .false.
   stdout_flag = .false.
   test_flag = .false.
   stats_flag = .false.
   mirror_flag = .false.
   maxlvl_flag = .false.
   back_flag = .false.

   maxrec = 1
   maxcount = 10
   maxlevel = 10
   maxcoord = 16
   bias_tol = 0.35
   bias_scale = 1.e3
   bias_ratio = 0.5
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
         sort_flag = .true.
      case ('-back')
         back_flag = .true.
      case ('-fast')
         iter_flag = .true.
         bias_flag = .true.
      case ('-bond')
         bond_flag = .true.
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
      case ('-scale')
         call readoptarg(arg, bias_scale)
      case ('-level')
         maxlvl_flag = .true.
         call readoptarg(arg, maxlevel)
      case ('-rec')
         call readoptarg(arg, maxrec)
      case ('-out')
         call readoptarg(arg, pathout)
      case ('-stdin')
         stdin_flag = .true.
         call readoptarg(arg, fmtin0)
         call readoptarg(arg, fmtin1)
      case ('-stdout')
         stdout_flag = .true.
         call readoptarg(arg, fmtout)
      case default
         call readposarg(arg, posargs)
      end select

   end do

   if (stdin_flag) then
      read_unit0 = input_unit
      read_unit1 = input_unit
   else
      select case (ipos)
      case (0)
         write (error_unit, '(a)') 'Error: No file paths were specified'
         stop
      case (1)
         pathin1 = trim(posargs(1))
         call open2read(pathin1, read_unit0, fmtin0)
         read_unit1 = read_unit0
         fmtin1 = fmtin0
      case (2)
         pathin1 = trim(posargs(1))
         pathin2 = trim(posargs(2))
         call open2read(pathin1, read_unit0, fmtin0)
         call open2read(pathin2, read_unit1, fmtin1)
      end select
   end if

   ! Read coordinates

   call readfile(read_unit0, fmtin0, title0, natom0, labels0, coords0, adjmat0)
   call readfile(read_unit1, fmtin1, title1, natom1, labels1, coords1, adjmat1)

   if (mirror_flag) then
      coords1(1, :) = -coords1(1, :)
   end if

   ! Allocate arrays

   allocate(znums0(natom0), znums1(natom1))
   allocate(types0(natom0), types1(natom1))
   allocate(weights0(natom0), weights1(natom1))
   allocate(aligned1(3, natom1))
   allocate(permlist(natom0, maxrec))
   allocate(countlist(maxrec))

   ! Get atomic numbers, types and weights

   do i = 1, natom0
      call readlabel(labels0(i), znums0(i), types0(i))
      weights0(i) = weight_func(znums0(i))
   end do

   do i = 1, natom1
      call readlabel(labels1(i), znums1(i), types1(i))
      weights1(i) = weight_func(znums1(i))
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
      call getadjmat(natom0, coords0, znums0, adjmat0)
      call getadjmat(natom1, coords1, znums1, adjmat1)
   else
      adjmat0(:, :) = .false.
      adjmat1(:, :) = .false.
   end if

   ! Sort atoms to minimize MSD

   if (sort_flag) then

      call assign_atoms( &
         natom0, &
         znums0, &
         types0, &
         weights0, &
         coords0, &
         adjmat0, &
         natom1, &
         znums1, &
         types1, &
         weights1, &
         coords1, &
         adjmat1, &
         permlist, &
         countlist, &
         nrec, &
         error)

      if (error /= 0) stop

      do i = 1, nrec

         call align_atoms( &
            natom0, &
            znums0, &
            types0, &
            weights0, &
            coords0, &
            natom1, &
            znums1(permlist(:, i)), &
            types1(permlist(:, i)), &
            weights1(permlist(:, i)), &
            coords1(:, permlist(:, i)), &
            travec, &
            rotmat, &
            error)

         if (error /= 0) stop

         adjd = adjacencydiff(natom0, adjmat0, adjmat1, permlist(:, i))
         aligned1 = translated(natom1, rotated(natom1, coords1, rotmat), travec)
         rmsd = sqrt(squaredist(natom0, weights0, coords0, aligned1, permlist(:, i))/sum(weights0))

         if (i == 1) then
            if (bond_flag) write (output_unit, '(a)') 'Optimized AdjD = ' // intstr(adjd)
            write (output_unit, '(a)') 'Optimized RMSD = ' // realstr(rmsd, 4)
            call writefile(write_unit, fmtout, 'Reference', natom0, znums0, coords0, adjmat0)
         end if

         call writefile(write_unit, fmtout, 'RMSD ' // realstr(rmsd, 4), natom1, znums1(permlist(:, i)), &
            aligned1(:, permlist(:, i)), adjmat1(permlist(:, i), permlist(:, i)))

      end do

   else

      ! Align atoms to minimize RMSD

      call align_atoms( &
         natom0, &
         znums0, &
         types0, &
         weights0, &
         coords0, &
         natom1, &
         znums1, &
         types1, &
         weights1, &
         coords1, &
         travec, &
         rotmat, &
         error)

      if (error /= 0) stop

      adjd = adjacencydiff(natom0, adjmat0, adjmat1, identityperm(natom0))
      aligned1 = translated(natom1, rotated(natom1, coords1, rotmat), travec)
      rmsd = sqrt(squaredist(natom0, weights0, coords0, aligned1, identityperm(natom0))/sum(weights0))

      if (bond_flag) write (output_unit, '(a)') 'AdjD = ' // intstr(adjd)
      write (output_unit, '(a)') 'RMSD = ' // realstr(rmsd, 4)
      call writefile(write_unit, fmtout, 'Reference', natom0, znums0, coords0, adjmat0)
      call writefile(write_unit, fmtout, 'RMSD ' // realstr(rmsd, 4), natom1, znums1, aligned1, adjmat1)

   end if

end program
