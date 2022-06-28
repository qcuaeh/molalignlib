module options

integer, parameter :: wp = 8
integer, parameter :: lbllen = 16
integer, parameter :: optlen = 256
integer, parameter :: ttllen = 256
integer, parameter :: maxcoord = 16
integer, parameter :: temp_file_unit = 1000
integer, parameter :: first_file_unit = 1001
integer, parameter :: second_file_unit = 1002

integer maxtrial, maxcount
logical live, iterative, biased, trialing, counting, testing
character(optlen) weighter, outformat
real(wp) biasscale, tolerance

end module
