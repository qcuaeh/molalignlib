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
use parameters
use settings

implicit none

contains

subroutine print_header()
   write (output_unit, '(1x,a,4x,a,4x,a,5x,a,6x,a,7x,a)') 'Map', 'Count', 'Steps', 'Total', 'Real', 'RMSD'
   write (output_unit, '(a)') '-----------------------------------------------------'
end subroutine

subroutine print_body(imap, matches, avgsteps, avgtotalrot, avgrealrot, dist2)
   integer, intent(in) :: imap, matches
   real(wp), intent(in) :: avgsteps, avgtotalrot, avgrealrot, dist2
   write (output_unit, '(i4,3x,i6,5x,f4.1,5x,f5.1,5x,f5.1,3x,f8.4)') &
      imap, matches, avgsteps, 90./asin(1.)*avgtotalrot, 90./asin(1.)*avgrealrot, sqrt(dist2)
end subroutine

subroutine print_footer()
   if (live_flag) write (output_unit, '(a)', advance='no') achar(27)//'[K'
   write (output_unit, '(a)') '-----------------------------------------------------'
end subroutine

subroutine print_stats(overflow, nrec, nmap, ntrial, nstep)
   logical, intent(in) :: overflow
   integer, intent(in) :: nrec, nmap, ntrial, nstep
   write (output_unit, '(a,1x,i0)') 'Random trials =', ntrial
   write (output_unit, '(a,1x,i0)') 'Minimization steps =', nstep
   if (overflow) then
      write (output_unit, '(a,1x,i0)') 'Visited local minima >', nrec
   else
      write (output_unit, '(a,1x,i0)') 'Visited local minima =', nmap
   end if
end subroutine

end module
