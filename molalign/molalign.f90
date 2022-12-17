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
   use library

   implicit none

   integer :: error
   integer :: i, nmap, nrec
   integer :: natom0, natom1
   integer :: unit, unit0, unit1
   integer, allocatable :: mapping(:, :)
   integer, allocatable :: mapcount(:)
   integer, allocatable, dimension(:) :: znums0, znums1, types0, types1
   character(title_len) :: title0, title1
   character(arg_len) :: arg, files(2), fmtstdin, fmtin0, fmtin1, fmtout
   character(label_len), allocatable, dimension(:) :: labels0, labels1
   real(wp) :: travec(3), rotmat(3, 3), dist2
   real(wp), allocatable :: mapdist2(:)
   real(wp), allocatable :: weights0(:), weights1(:)
   real(wp), dimension(:, :), allocatable :: coords0, coords1, aligned1
   logical :: sort_flag, enan_flag, stdin_flag, test_flag

   procedure(f_realint), pointer :: weight_function

   ! Set default options

   bias_flag = .false.
   iter_flag = .false.
   sort_flag = .false.
   trial_flag = .false.
   stdin_flag = .false.
   repro_flag = .false.
   enan_flag = .false.
   test_flag = .false.
   live_flag = .false.

   nrec = 1
   maxcount = 10
   biastol = 0.35
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
         repro_flag = .true.
      case ('-sort')
         sort_flag = .true.
      case ('-fast')
         bias_flag = .true.
         iter_flag = .true.
      case ('-mass')
         weight_function => stdmass
      case ('-enan')
         enan_flag = .true.
      case ('-count')
         call readoptarg(arg, maxcount)
      case ('-trials')
         trial_flag = .false.
         call readoptarg(arg, maxtrials)
      case ('-tol')
         call readoptarg(arg, biastol)
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

   if (enan_flag) then
      coords1(1, :) = -coords1(1, :)
   end if

   ! Allocate arrays

   allocate(znums0(natom0), znums1(natom1))
   allocate(types0(natom0), types1(natom1))
   allocate(weights0(natom0), weights1(natom1))
   allocate(aligned1(3, natom1))
   allocate(mapping(natom0, nrec))
   allocate(mapcount(nrec))
   allocate(mapdist2(nrec))

   ! Get atomic numbers, types and weights

   do i = 1, natom0
      call readlabel(labels0(i), znums0(i), types0(i))
      weights0(i) = weight_function(znums0(i))
   end do

   do i = 1, natom1
      call readlabel(labels1(i), znums1(i), types1(i))
      weights1(i) = weight_function(znums1(i))
   end do

   ! Sort atoms to minimize MSD

   if (sort_flag) then

      call assign_atoms( &
         natom0, &
         znums0, &
         types0, &
         coords0, &
         weights0, &
         natom1, &
         znums1, &
         types1, &
         coords1, &
         weights1, &
         nrec, &
         nmap, &
         mapping, &
         mapcount, &
         mapdist2, &
         error)

      if (error /= 0) stop

      do i = 1, nmap

         call align_atoms( &
            natom0, &
            znums0, &
            types0, &
            coords0, &
            weights0, &
            natom1, &
            znums1(mapping(:, i)), &
            types1(mapping(:, i)), &
            coords1(:, mapping(:, i)), &
            weights1, &
            travec, &
            rotmat, &
            dist2, &
            error)

         if (error /= 0) stop

         aligned1 = translated(natom1, rotated(natom1, coords1, rotmat), travec)

         if (.not. test_flag) then

            call open2write('aligned_'//str(i)//'.'//trim(fmtout), unit)
            call writefile(unit, fmtout, natom0, title0, znums0, coords0)
            call writefile(unit, fmtout, natom1, title1, znums1(mapping(:, i)), aligned1(:, mapping(:, i)))
            close(unit)

         end if

      end do

   else

      ! Align atoms to minimize RMSD

      call align_atoms( &
         natom0, &
         znums0, &
         types0, &
         coords0, &
         weights0, &
         natom1, &
         znums1, &
         types1, &
         coords1, &
         weights1, &
         travec, &
         rotmat, &
         dist2, &
         error)

      if (error /= 0) stop

      aligned1 = translated(natom1, rotated(natom1, coords1, rotmat), travec)
      write (output_unit, '(a,1x,a,1x,a)') 'RMSD =', str(sqrt(dist2)), '(only alignment performed)'

      if (.not. test_flag) then

         call open2write('aligned.'//trim(fmtout), unit)
         call writefile(unit, fmtout, natom0, title0, znums0, coords0)
         call writefile(unit, fmtout, natom1, title1, znums1, aligned1) 
         close(unit)

      end if

   end if

end program
