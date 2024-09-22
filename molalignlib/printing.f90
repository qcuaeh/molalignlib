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

module printing
use stdio
use kinds
use flags

implicit none

abstract interface
   subroutine print_stats_proc(nrec, matches, avgsteps, avgtotalrot, avgrealrot, recadjd, recrmsd)
      use kinds
      integer, intent(in) :: nrec
      integer, dimension(:), intent(in) :: matches, recadjd
      real(rk), dimension(:), intent(in) :: avgsteps, avgtotalrot, avgrealrot, recrmsd
   end subroutine
end interface

procedure(print_stats_proc), pointer :: print_stats

character(*), parameter :: line1 = '-------------------------------------------'
character(*), parameter :: line2 = '--------------------------------------------------'

contains

subroutine print_stats_dist(nrec, matches, avgsteps, avgtotalrot, avgrealrot, recadjd, recrmsd)
   integer, intent(in) :: nrec
   integer, dimension(:), intent(in) :: matches, recadjd
   real(rk), dimension(:), intent(in) :: avgsteps, avgtotalrot, avgrealrot, recrmsd
   integer :: irec
   write (stderr, '(1x,a,4x,a,4x,a,5x,a,7x,a)') 'Map', 'Count', 'Steps', 'Angle', 'RMSD'
   write (stderr, '(a)') line1
   do irec = 1, nrec
      write (stderr, '(i4,3x,i6,5x,f4.1,5x,f5.1,3x,f8.4)') &
         irec, matches(irec), avgsteps(irec), 90./asin(1.)*avgrealrot(irec), recrmsd(irec)
   end do
   write (stderr, '(a)') line1
   flush(stderr)
end subroutine

subroutine print_stats_diff(nrec, matches, avgsteps, avgtotalrot, avgrealrot, recadjd, recrmsd)
   integer, intent(in) :: nrec
   integer, dimension(:), intent(in) :: matches, recadjd
   real(rk), dimension(:), intent(in) :: avgsteps, avgtotalrot, avgrealrot, recrmsd
   integer :: irec
   write (stderr, '(1x,a,4x,a,4x,a,5x,a,3x,a,7x,a)') 'Map', 'Count', 'Steps', 'Angle', 'AdjD', 'RMSD'
   write (stderr, '(a)') line2
   do irec = 1, nrec
      write (stderr, '(i4,3x,i6,5x,f4.1,5x,f5.1,3x,i4,3x,f8.4)') &
         irec, matches(irec), avgsteps(irec), 90./asin(1.)*avgrealrot(irec), recadjd(irec), recrmsd(irec)
   end do
   write (stderr, '(a)') line2
   flush(stderr)
end subroutine

subroutine print_final_stats(overflow, maxrec, nrec, ntrial, nstep)
   logical, intent(in) :: overflow
   integer, intent(in) :: maxrec, nrec, ntrial, nstep
   write (stderr, '(a,1x,i0)') 'Random trials =', ntrial
   write (stderr, '(a,1x,i0)') 'Minimization steps =', nstep
   if (overflow) then
      write (stderr, '(a,1x,i0)') 'Visited local minima >', maxrec
   else
      write (stderr, '(a,1x,i0)') 'Visited local minima =', nrec
   end if
   flush(stderr)
end subroutine

end module
