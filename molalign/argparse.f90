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
public get_arg
public read_posarg
public read_optarg
public init_args

type, public :: p_char
   character(:), pointer :: var
end type

interface read_optarg
   module procedure int_read_optarg
   module procedure real_read_optarg
   module procedure char_read_optarg
end interface

contains

subroutine init_args()

   iarg = 0
   ipos = 0

end subroutine

subroutine read_posarg(arg, posargs)
   character(*), intent(in) :: arg
   type(p_char), intent(out) :: posargs(:)

   if (arg(1:1) == '-') then
      write (stderr, '(a,1x,a)') 'Unknown option:', arg
      stop
   end if

   ipos = ipos + 1

   if (ipos > size(posargs)) then
      write (stderr, '(a)') 'Too many positional arguments'
      stop
   end if

   posargs(ipos)%var = arg

end subroutine

logical function get_arg(arg) result(success)
   character(:), allocatable, intent(out) :: arg
   integer :: arglen

   iarg = iarg + 1

   if (iarg <= command_argument_count()) then
      call get_command_argument(iarg, length=arglen)
      allocate (character(arglen) :: arg)
      call get_command_argument(iarg, arg)
      success = .true.
   else
      success = .false.
   end if

end function

subroutine get_optarg(option, optarg)
   character(*), intent(in) :: option
   character(:), allocatable, intent(out) :: optarg
   integer :: arglen

   iarg = iarg + 1

   if (iarg > command_argument_count()) then
      write (stderr, '(a,1x,a,1x,a)') 'Option', option, 'requires an argument'
      stop
   else
      call get_command_argument(iarg, length=arglen)
      allocate (character(arglen) :: optarg)
      call get_command_argument(iarg, optarg)
      if (optarg(1:1) == '-') then
         write (stderr, '(a,1x,a,1x,a)') 'Option', option, 'requires an argument'
         stop
      end if
   end if

end subroutine

subroutine char_read_optarg(option, optarg)
   character(*), intent(in) :: option
   character(:), allocatable, intent(out) :: optarg

   call get_optarg(option, optarg)

end subroutine

subroutine int_read_optarg(option, optarg)
   character(*), intent(in) :: option
   integer, intent(out) :: optarg
   character(:), allocatable :: chararg 
   integer :: stat

   call get_optarg(option, chararg)
   read (chararg, *, iostat=stat) optarg
   if (stat /= 0) then
      write (stderr, '(a,1x,a,1x,a)') 'Option', option, 'requires an integer argument'
      stop
   end if

end subroutine

subroutine real_read_optarg(option, optarg)
   character(*), intent(in) :: option
   real(rk), intent(out) :: optarg
   character(:), allocatable :: chararg 
   integer :: stat

   call get_optarg(option, chararg)
   read (chararg, *, iostat=stat) optarg
   if (stat /= 0) then
      write (stderr, '(a,1x,a,1x,a)') 'Option', option, 'requires a numeric argument'
      stop
   end if

end subroutine

end module
