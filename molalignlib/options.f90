module options

integer, parameter :: arg_len = 256
integer, parameter :: title_len = 256
integer, parameter :: label_len = 16
integer, parameter :: maxcoord = 16

logical bias_flag, iter_flag, test_flag, live_flag, stop_flag, debug_flag
integer max_count, max_trials
real bias_tol, bias_scale

end module
