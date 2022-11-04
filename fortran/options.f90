module options

integer, parameter :: arg_len = 256
integer, parameter :: title_len = 256
integer, parameter :: label_len = 16
integer, parameter :: maxcoord = 16
logical, parameter :: conv_flag = .true.

integer maxcount, maxtrials
logical bias_flag, halt_flag, test_flag, debug_flag, live_flag
real bias_scale, bias_tol

end module
