module pointers
use kinds

abstract interface
   real(wp) function f_realint(z)
      use kinds
      integer, intent(in) :: z
   end function
end interface

procedure(f_realint), pointer :: weight_func

end module
