module options

integer, parameter :: wp = 8
integer, parameter :: lbllen = 16
integer, parameter :: optlen = 128
integer, parameter :: ttllen = 256
integer, parameter :: maxcoord = 16
integer, parameter :: first_file_unit = 1001
integer, parameter :: second_file_unit = 1002
integer, parameter :: third_file_unit = 1003

integer maxtrial, maxmatch
logical live, iterative, biased, trialing, matching, testing
character(optlen) weighter, outformat
real(wp) biasscale, tolerance

end module
