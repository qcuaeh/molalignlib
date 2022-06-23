module options

real scaling, tolerance
integer maxcoord, maxcount, maxlevel, maxrecord
logical live, iterative, biased, converge, testing
character(32) weighting, output_format

procedure (generic_test), pointer :: stop_test => null()

abstract interface
    logical function generic_test(x, y)
        integer, intent(in) :: x, y
    end function
end interface

end module
