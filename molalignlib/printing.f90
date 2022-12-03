module printing
use parameters
use settings

implicit none

contains

subroutine print_header()
    if (live_flag) write (output_unit, '(a)', advance='no') achar(27)//'[K'
    write (output_unit, '(1x,a,4x,a,4x,a,5x,a,6x,a,7x,a)') 'Map', 'Count', 'Steps', 'Total', 'Real', 'RMSD'
    if (live_flag) write (output_unit, '(a)', advance='no') achar(27)//'[K'
    write (output_unit, '(a)') '-----------------------------------------------------'
end subroutine

subroutine print_stats(imap, matches, avgiter, avgtotalrot, avgrealrot, rmsd)
    integer, intent(in) :: imap, matches
    real(wp), intent(in) :: avgiter, avgtotalrot, avgrealrot, rmsd
    if (live_flag) write (output_unit, '(a)', advance='no') achar(27)//'[K'
    write (output_unit, '(i4,3x,i6,5x,f4.1,5x,f5.1,5x,f5.1,3x,f8.4)') &
        imap, matches, avgiter, 90./asin(1.)*avgtotalrot, 90./asin(1.)*avgrealrot, rmsd
end subroutine

subroutine print_footer(overflow, nmap, itrial)
    logical, intent(in) :: overflow
    integer, intent(in) :: nmap, itrial
    if (live_flag) write (output_unit, '(a)', advance='no') achar(27)//'[K'
    write (output_unit, '(a)') '-----------------------------------------------------'
    if (live_flag) write (output_unit, '(a)', advance='no') achar(27)//'[K'
    if (overflow) then
        write (output_unit, '(a,1x,i0,1x,a,1x,i0,1x,a)') 'Found more than', nmap, 'mapping(s) in', &
            itrial, 'random trial(s)'
    else
        write (output_unit, '(a,1x,i0,1x,a,1x,i0,1x,a)') 'Found', nmap, 'mapping(s) in', itrial, 'random trial(s)'
    end if
end subroutine

end module