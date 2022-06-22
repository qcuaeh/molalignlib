module globals

real scaling, levelratio, tolerance
integer maxcoord, maxcount, maxlevel, maxrecord
logical topology, iterative, backtrack, permutate, debug, live
character(32) weighting, output_format
character(128) output_file

procedure (p_biasing_test), pointer :: bias_test => null()
procedure (p_counting_test), pointer :: stop_test => null()

abstract interface
    logical function p_counting_test(x, y)
        integer, intent(in) :: x, y
    end function
end interface

abstract interface
    logical function p_biasing_test(x)
        real, dimension(:), intent(in) :: x
    end function
end interface

end module
