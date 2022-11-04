module sorting

use iso_fortran_env, only: error_unit

use options
use maputils

implicit none

private
public sort
public order

interface sort
    module procedure intquicksort
    module procedure realquicksort
end interface

interface order
    module procedure intorder
    module procedure charorder
end interface

contains

function intorder(x, n) result(o)
    integer, intent(in) :: n
    integer, dimension(:), intent(in) :: x
    integer o(n), t((n+1)/2)
    call initializemap(o)
    call intmergesort(x, o, n, t)
end function

function charorder(x, n) result(o)
    integer, intent(in) :: n
    character(*), dimension(:), intent(in) :: x
    integer o(n), t((n+1)/2)
    call initializemap(o)
    call charmergesort(x, o, n, t)
end function

recursive subroutine intquicksort(x, m, n)
    integer, intent(in) :: m, n
    integer, dimension(:), intent(inout) :: x

    integer i, j
    integer xi, xp

    xp = x((m+n)/2)
    i = m
    j = n

    do
        do while (x(i) < xp)
            i = i + 1
        end do
        do while (xp < x(j))
            j = j - 1
        end do
        if (i >= j) exit
        xi = x(i); x(i) = x(j); x(j) = xi
        i = i + 1
        j = j - 1
    end do

    if (m < i-1) call intquicksort(x, m, i-1)
    if (j+1 < n) call intquicksort(x, j+1, n)
end subroutine

recursive subroutine realquicksort(x, m, n)
    integer, intent(in) :: m, n
    real, dimension(:), intent(inout) :: x

    integer i, j
    real xi, xp

    xp = x((m+n)/2)
    i = m
    j = n

    do
        do while (x(i) < xp)
            i = i + 1
        end do
        do while (xp < x(j))
            j = j - 1
        end do
        if (i >= j) exit
        xi = x(i); x(i) = x(j); x(j) = xi
        i = i + 1
        j = j - 1
    end do

    if (m < i-1) call realquicksort(x, m, i-1)
    if (j+1 < n) call realquicksort(x, j+1, n)
end subroutine

subroutine intmerge(x, a, na, b, nb, c, nc)
    integer, intent(in) :: na, nb, nc
    integer, intent(in) :: x(:)
    integer, intent(in) :: b(nb)
    integer, intent(inout) :: a(na), c(nc)
     
    integer :: i, j, k
     
    i = 1; j = 1; k = 1;
    do while(i <= na .and. j <= nb)
        if (x(a(i)) <= x(b(j))) then
            c(k) = a(i)
            i = i+1
        else
            c(k) = b(j)
            j = j+1
        end if
        k = k + 1
    end do
    do while (i <= na)
        c(k) = a(i)
        i = i + 1
        k = k + 1
    end do
end subroutine
 
recursive subroutine intmergesort(x, o, n, t)
    integer, intent(in) :: n
    integer, dimension(:), intent(in) :: x
    integer, dimension(n), intent(inout) :: o
    integer, dimension((n+1)/2), intent(out) :: t
     
    integer :: i, o1
     
    if (n < 2) return
    if (n == 2) then
        if (x(o(1)) > x(o(2))) then
            o1 = o(1); o(1) = o(2); o(2) = o1
        end if
        return
    end if      
    i = (n+1)/2
     
    call intmergesort(x, o, i, t)
    call intmergesort(x, o(i+1), n-i, t)
     
    if (x(o(i)) > x(o(i+1))) then
        t(1:i) = o(1:i)
        call intmerge(x, t, i, o(i+1), n-i, o, n)
    end if
end subroutine

subroutine charmerge(x, a, na, b, nb, c, nc)
    integer, intent(in) :: na, nb, nc
    character(*), intent(in) :: x(:)
    integer, intent(in) :: b(nb)
    integer, intent(inout) :: a(na), c(nc)
     
    integer :: i, j, k
     
    i = 1; j = 1; k = 1;
    do while(i <= na .and. j <= nb)
        if (x(a(i)) <= x(b(j))) then
            c(k) = a(i)
            i = i+1
        else
            c(k) = b(j)
            j = j+1
        end if
        k = k + 1
    end do
    do while (i <= na)
        c(k) = a(i)
        i = i + 1
        k = k + 1
    end do
end subroutine
 
recursive subroutine charmergesort(x, o, n, t)
    integer, intent(in) :: n
    character(*), dimension(:), intent(in) :: x
    integer, dimension(n), intent(inout) :: o
    integer, dimension((n+1)/2), intent(out) :: t
     
    integer :: i, o1
     
    if (n < 2) return
    if (n == 2) then
        if (x(o(1)) > x(o(2))) then
            o1 = o(1); o(1) = o(2); o(2) = o1
        end if
        return
    end if      
    i = (n+1)/2
     
    call charmergesort(x, o, i, t)
    call charmergesort(x, o(i+1), n-i, t)
     
    if (x(o(i)) > x(o(i+1))) then
        t(1:i) = o(1:i)
        call charmerge(x, t, i, o(i+1), n-i, o, n)
    end if
end subroutine

end module
