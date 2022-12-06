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

   integer :: i, nmap, nrec, stat
   integer :: natom0, natom1
   integer :: unit, unit0, unit1
   integer, allocatable :: maplist(:, :)
   integer, allocatable :: countlist(:)
   integer, allocatable, dimension(:) :: znums0, znums1, types0, types1
   character(arg_len) :: arg, files(2), fmtstdin, fmtin0, fmtin1, fmtout
   character(title_len) :: title0, title1
   character(label_len), allocatable, dimension(:) :: labels0, labels1
   real(wp) :: rmsd
   real(wp), allocatable :: rmsdlist(:)
   real(wp), allocatable :: weights0(:)
   real(wp), dimension(:, :), allocatable :: coords0, coords1, aligned1
   logical :: sort_flag, stdin_flag

   procedure(f_realint), pointer :: weight_function

   ! Set default options

   bias_flag = .false.
   iter_flag = .false.
   sort_flag = .false.
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
!        case ('-biasonly')
!            bias_flag = .true.
!            iter_flag = .false.
!        case ('-iteronly')
!            bias_flag = .false.
!            iter_flag = .true.
      case ('-mass')
         weight_function => stdmass
      case ('-count')
         call readoptarg(arg, max_count)
      case ('-trials')
         free_flag = .false.
         call readoptarg(arg, max_trials)
      case ('-tol')
         call readoptarg(arg, bias_tol)
!        case ('-scale')
!            call readoptarg(arg, bias_scale)
      case ('-rec')
         call readoptarg(arg, nrec)
      case ('-out')
         call readoptarg(arg, fmtout)
      case ('-stdin')
         stdin_flag = .true.
         call readoptarg(arg, fmtstdin)
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

   ! Allocate arrays

   allocate(znums0(natom0), znums1(natom1))
   allocate(types0(natom0), types1(natom1))
   allocate(weights0(natom0))
   allocate(aligned1(3, natom1))
   allocate(maplist(natom0, nrec))
   allocate(countlist(nrec))
   allocate(rmsdlist(nrec))

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
         weights0, nrec, nmap, maplist, countlist, rmsdlist)

      ! Write aligned coordinates

      do i = 1, nmap

         call align_atoms(natom0, natom1, znums0, znums1(maplist(:, i)), &
            types0, types1(maplist(:, i)), coords0, coords1(:, maplist(:, i)), &
            weights0, rmsd, aligned1)

         call open2write('aligned_'//str(i)//'.'//trim(fmtout), unit)
         call writefile(unit, fmtout, natom0, title0, znums0, coords0)
         call writefile(unit, fmtout, natom1, title1, znums1(maplist(:, i)), aligned1)
         close(unit)

      end do

   else

      ! Align atoms to minimize RMSD

      call align_atoms(natom0, natom1, znums0, znums1, types0, types1, coords0, coords1, &
         weights0, rmsd, aligned1)

      ! Write aligned coordinates

      write (output_unit, '(a,1x,f0.4,1x,a)') 'RMSD:', rmsd, '(only alignment performed)'

      call open2write('aligned_0.'//trim(fmtout), unit)
      call writefile(unit, fmtout, natom0, title0, znums0, coords0)
      call writefile(unit, fmtout, natom1, title1, znums1, aligned1) 
      close(unit)

   end if

end program
