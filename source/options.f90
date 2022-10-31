module options

integer, parameter :: arg_len = 256
integer, parameter :: title_len = 256
integer, parameter :: label_len = 16
integer, parameter :: maxcoord = 16

integer maxcount, maxtrials
logical sort_flag, bias_flag, conv_flag, halt_flag, test_flag, debug_flag, live_flag
real bias_scale, bias_tol

end module
