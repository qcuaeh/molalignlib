module options

integer, parameter :: arg_len = 256
integer, parameter :: title_len = 256
integer, parameter :: label_len = 16
integer, parameter :: maxcoord = 16

integer convcount, maxtrials
logical biased, iterated, terminate, converge, testing, debug, live
real bias_scale, bias_tol

end module
