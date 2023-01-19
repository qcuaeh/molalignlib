module settings
use parameters

implicit none
integer, parameter :: maxcoord = 16
logical :: bias_flag, iter_flag, trial_flag, repro_flag, stats_flag, live_flag
integer :: maxcount, maxtrials
real(wp) :: bias_tol, bias_scale

end module
