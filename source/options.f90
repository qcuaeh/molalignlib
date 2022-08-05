module options

integer, parameter :: arg_len = 256
integer, parameter :: title_len = 256
integer, parameter :: label_len = 16
integer, parameter :: temp_file_unit = 1000
integer, parameter :: first_file_unit = 1001
integer, parameter :: second_file_unit = 1002
integer, parameter :: maxcoord = 16

integer maxtrial, mincount
logical live, iterative, biased, bounded, converged, testing
character(arg_len) formatout
real lenscale, tolerance

end module
