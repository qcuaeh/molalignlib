module options

integer, parameter :: arg_len = 256
integer, parameter :: title_len = 256
integer, parameter :: label_len = 16
integer, parameter :: maxcoord = 16

integer convcount, maxtrial
logical biased, iterated, complete, converge, testing, live
real biasscale, biastol

end module
