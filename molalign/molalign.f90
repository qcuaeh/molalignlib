! MolAlign
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
   use settings
   use strutils
   use optparse
   use chemutils
   use readmol
   use writemol
   use translation
   use rotation
   use alignment
   use molalignlib

   implicit none

   integer :: error
   integer :: i, nmap, nrec
   integer :: natom0, natom1
   integer :: unit, unit0, unit1
   integer, allocatable :: mapind(:, :)
   integer, allocatable :: mapcount(:)
   integer, allocatable, dimension(:) :: znums0, znums1, types0, types1
   character(arg_len) :: arg, files(2), fmtstdin, fmtin0, fmtin1, fmtout
   character(title_len) :: title0, title1
   character(label_len), allocatable, dimension(:) :: labels0, labels1
   real(wp) :: dist
   real(wp), allocatable :: mapdist(:)
   real(wp), allocatable :: weights0(:)
   real(wp), dimension(:, :), allocatable :: coords0, coords1, aligned1
   logical :: sort_flag, mirror_flag, stdin_flag

   procedure(f_realint), pointer :: weight_function

   ! Set default options

   bias_flag = .false.
   iter_flag = .false.
   sort_flag = .false.
   mirror_flag = .false.
   free_flag = .true.
   test_flag = .false.
   live_flag = .false.
   debug_flag = .false.
   stdin_flag = .false.

   nrec = 1
   max_count = 10
   bias_tol = 0.35
   bias_scale = 1.e3
   fmtstdin = 'xyz'
   fmtout = 'xyz'

   weight_function => unity

   ! Get user options

   call initarg()

   do while (getarg(arg))

      select case (arg)
      case ('-live')
         live_flag = .true.
      case ('-test')
         test_flag = .true.
      case ('-sort')
         sort_flag = .true.
      case ('-debug')
         debug_flag = .true.
      case ('-fast')
         bias_flag = .true.
         iter_flag = .true.
      case ('-mass')
         weight_function => stdmass
      case ('-mirror')
         mirror_flag = .true.
      case ('-count')
         call readoptarg(arg, max_count)
      case ('-trials')
         free_flag = .false.
         call readoptarg(arg, max_trials)
      case ('-tol')
         call readoptarg(arg, bias_tol)
      case ('-rec')
         call readoptarg(arg, nrec)
      case ('-out')
         call readoptarg(arg, fmtout)
      case ('-stdin')
         stdin_flag = .true.
         call readoptarg(arg, fmtstdin)
!      case ('-bias')
!         bias_flag = .true.
!         iter_flag = .false.
!      case ('-iter')
!         bias_flag = .false.
!         iter_flag = .true.
!      case ('-none')
!         bias_flag = .false.
!         iter_flag = .false.
!      case ('-scale')
!         call readoptarg(arg, bias_scale)
      case default
         call readarg(arg, files)
      end select

   end do

   if (stdin_flag) then
      unit0 = input_unit
      unit1 = input_unit
      fmtin0 = fmtstdin
      fmtin1 = fmtstdin
   else
      select case (ipos)
      case (0)
         write (error_unit, '(a)') 'Error: No file paths were specified'
         stop
      case (1)
         call open2read(files(1), unit0, fmtin0)
         unit1 = unit0
         fmtin1 = fmtin0
      case (2)
         call open2read(files(1), unit0, fmtin0)
         call open2read(files(2), unit1, fmtin1)
      end select
   end if

   ! Read coordinates

   call readfile(unit0, fmtin0, natom0, title0, labels0, coords0)
   call readfile(unit1, fmtin1, natom1, title1, labels1, coords1)

   if (mirror_flag) then
      coords1(1, :) = -coords1(1, :)
   end if

   ! Allocate arrays

   allocate(znums0(natom0), znums1(natom1))
   allocate(types0(natom0), types1(natom1))
   allocate(weights0(natom0))
   allocate(aligned1(3, natom1))
   allocate(mapind(natom0, nrec))
   allocate(mapcount(nrec))
   allocate(mapdist(nrec))

   ! Get atomic numbers and types

   do i = 1, natom0
      call readlabel(labels0(i), znums0(i), types0(i))
   end do

   do i = 1, natom1
      call readlabel(labels1(i), znums1(i), types1(i))
   end do

   ! Set weights

   do i = 1, natom0
      weights0(i) = weight_function(znums0(i))
   end do

   ! Normalize weights

   weights0 = weights0/sum(weights0)

   ! Sort atoms to minimize MSD

   if (sort_flag) then

      call assign_atoms(natom0, natom1, znums0, znums1, types0, types1, coords0, coords1, &
         weights0, nrec, nmap, mapind, mapcount, mapdist, error)

      if (error /= 0) stop

      ! Write aligned coordinates

      do i = 1, nmap

         call align_atoms(natom0, natom1, znums0, znums1(mapind(:, i)), &
            types0, types1(mapind(:, i)), coords0, coords1(:, mapind(:, i)), &
            weights0, aligned1, dist, error)

         if (error /= 0) stop

         call open2write('aligned_'//str(i)//'.'//trim(fmtout), unit)
         call writefile(unit, fmtout, natom0, title0, znums0, coords0)
         call writefile(unit, fmtout, natom1, title1, znums1(mapind(:, i)), aligned1)
         close(unit)

      end do

   else

      ! Align atoms to minimize RMSD

      call align_atoms(natom0, natom1, znums0, znums1, types0, types1, coords0, coords1, &
         weights0, aligned1, dist, error)

      if (error /= 0) stop

      ! Write aligned coordinates

      write (output_unit, '(a,1x,f0.4,1x,a)') 'RMSD:', sqrt(dist), '(only alignment performed)'

      call open2write('aligned_0.'//trim(fmtout), unit)
      call writefile(unit, fmtout, natom0, title0, znums0, coords0)
      call writefile(unit, fmtout, natom1, title1, znums1, aligned1) 
      close(unit)

   end if

end program
