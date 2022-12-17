module settings
use parameters
implicit none

logical :: bias_flag, iter_flag, trial_flag, repro_flag, live_flag
integer :: maxcount, maxtrials
real(wp) :: biastol

integer, parameter :: maxcoord = 16
real(wp), parameter :: biascale = 1.E3

end module
