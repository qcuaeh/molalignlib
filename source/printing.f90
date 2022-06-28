module printing

use iso_fortran_env, only: output_unit

use options

implicit none

contains

subroutine print_header()
    write (output_unit, '(a)', advance='no') achar(27)//'[K'
    write (output_unit, '(1x, a, 3x, a, 3x, a, 3x, a, 3x, a, 4x, a, 5x, a)') &
        'Map', 'Trial', 'Count', 'Cycles', 'Meanrot', 'Totalrot', 'wRMSD'
    write (output_unit, '(a)', advance='no') achar(27)//'[K'
    write (output_unit, '(a)') '-------------------------------------------------------------'
end subroutine

subroutine print_stats(imap, earliest, matches, avgiter, avgmeanrot, avgangle, mindist)
    integer, intent(in) :: imap, earliest, matches
    real(wp), intent(in) :: avgiter, avgmeanrot, avgangle, mindist
    write (output_unit, '(a)', advance='no') achar(27)//'[K'
    write (output_unit, '(i4, 4x, i4, 2x, i6, 5x, f4.1, 5x, f5.1, 5x, f5.1, 3x, f9.4)') &
        imap, earliest, matches, avgiter, 90./asin(1.)*avgmeanrot, 90./asin(1.)*avgangle, sqrt(mindist)
end subroutine

subroutine print_footer(overflow, nrecord, itrial)
    logical, intent(in) :: overflow
    integer, intent(in) :: nrecord, itrial
    write (output_unit, '(a)', advance='no') achar(27)//'[K'
    write (output_unit, '(a)') '-------------------------------------------------------------'
    write (output_unit, '(a)', advance='no') achar(27)//'[K'
    if (overflow) then
        write (output_unit, '(a,x,i0,x,a,x,i0,x,a)') 'Found more than', nrecord, 'mapping(s) in', &
            itrial, 'random trial(s)'
    else
        write (output_unit, '(a,x,i0,x,a,x,i0,x,a)') 'Found', nrecord, 'mapping(s) in', itrial, 'random trial(s)'
    end if
end subroutine

end module
