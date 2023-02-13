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

module argparse
use stdio
use kinds
use bounds

implicit none
integer :: iarg, ipos

private
public ipos
public getarg
public readposarg
public readoptarg
public initarg

interface readoptarg
   module procedure readstroptarg
   module procedure readintoptarg
   module procedure readrealoptarg
end interface

contains

subroutine initarg()

   iarg = 0
   ipos = 0

end subroutine

subroutine readposarg(arg, posargs)
   character(*), intent(in) :: arg
   character(ll), intent(out) :: posargs(:)

   if (arg(1:1) == '-') then
      write (error_unit, '(a,1x,a)') 'Unknown option:', arg
      stop
   end if

   ipos = ipos + 1

   if (ipos > size(posargs)) then
      write (error_unit, '(a)') 'Too many positional arguments'
      stop
   end if

   posargs(ipos) = arg

end subroutine

logical function getarg(arg)
   character(:), allocatable, intent(out) :: arg
   integer arglen

   iarg = iarg + 1

   if (iarg <= command_argument_count()) then
      call get_command_argument(iarg, length=arglen)
      allocate(character(arglen) :: arg)
      call get_command_argument(iarg, arg)
      getarg = .true.
   else
      getarg = .false.
   end if

end function

subroutine getoptarg(option, arg)
   character(*), intent(in) :: option
   character(:), allocatable, intent(out) :: arg
   integer arglen

   iarg = iarg + 1

   if (iarg <= command_argument_count()) then
      call get_command_argument(iarg, length=arglen)
      allocate(character(arglen) :: arg)
      call get_command_argument(iarg, arg)
      if (arg(1:1) /= '-') then
         return
      end if
   end if

   write (error_unit, '(a,1x,a,1x,a)') 'Option', option, 'requires an argument'
   stop

end subroutine

subroutine readstroptarg(option, optval)
   character(*), intent(in) :: option
   character(:), allocatable, intent(out) :: optval
   character(:), allocatable :: arg 

   call getoptarg(option, arg)
   optval = arg

end subroutine

subroutine readintoptarg(option, optval)
   character(*), intent(in) :: option
   integer, intent(out) :: optval
   character(:), allocatable :: arg 
   integer :: stat

   call getoptarg(option, arg)
   read (arg, *, iostat=stat) optval
   if (stat /= 0) then
      write (error_unit, '(a,1x,a,1x,a)') 'Option', option, 'requires an integer argument'
      stop
   end if

end subroutine

subroutine readrealoptarg(option, optval)
   character(*), intent(in) :: option
   real(wp), intent(out) :: optval
   character(:), allocatable :: arg 
   integer :: stat

   call getoptarg(option, arg)
   read (arg, *, iostat=stat) optval
   if (stat /= 0) then
      write (error_unit, '(a,1x,a,1x,a)') 'Option', option, 'requires a numeric argument'
      stop
   end if

end subroutine

end module
