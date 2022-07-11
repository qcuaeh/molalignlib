module options

integer, parameter :: lbllen = 16
integer, parameter :: optlen = 256
integer, parameter :: ttllen = 256
integer, parameter :: maxcoord = 16
integer, parameter :: temp_file_unit = 1000
integer, parameter :: first_file_unit = 1001
integer, parameter :: second_file_unit = 1002

integer maxtrials, maxcount
logical live, iterative, biased, aborting, counting, testing
character(optlen) outformat
real lenscale, tolerance

end module
