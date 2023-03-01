module subset
!https://people.sc.fsu.edu/~jburkardt/m_src/subset/subset.html
contains

subroutine perm1_next ( n, p, more, even )

!*****************************************************************************80
!
!! PERM1_NEXT computes permutations of (1,...,N), one at a time.
!
!  Discussion:
!
!    The routine is initialized by calling with MORE = TRUE, in which case
!    it returns the identity permutation.
!
!    If the routine is called with MORE = FALSE, then the successor of the
!    input permutation is computed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2001
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer P(N), the permutation, in standard
!    index form.  On the first call, the input value is unimportant.
!    On subsequent calls, the input value should be the same
!    as the output value from the previous call.  In other words, the
!    user should just leave P alone.
!    On output, contains the "next" permutation.
!
!    Input/output, logical MORE.
!    Set MORE = FALSE before the first call.
!    MORE will be reset to TRUE and a permutation will be returned.
!    Each new call produces a new permutation until
!    MORE is returned FALSE.
!
!    Input/output, logical EVEN.
!    The input value of EVEN should simply be its output value from the
!    previous call; (the input value on the first call doesn't matter.)
!    On output, EVEN is TRUE if the output permutation is even, that is,
!    involves an even number of transpositions.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  logical even
  integer i
  integer i1
  integer ia
  integer id
  integer is
  integer j
  integer l
  integer m
  logical more
  integer p(n)

  if ( .not. more ) then

    call i4vec_indicator1 ( n, p )
    more = .true.
    even = .true.

    if ( n == 1 ) then
      more = .false.
      return
    end if

    if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
      return
    end if

    do i = 1, n - 3
      if ( p(i+1) /= p(i)+1 ) then
        return
      end if
    end do

    more = .false.

  else

    if ( n == 1 ) then
      p(1) = 0
      more = .false.
      return
    end if

    if ( even ) then

      ia = p(1)
      p(1) = p(2)
      p(2) = ia
      even = .false.

      if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
        return
      end if

      do i = 1, n - 3
        if ( p(i+1) /= p(i)+1 ) then
          return
        end if
      end do

      more = .false.
      return

    else

      more = .false.

      is = 0

      do i1 = 2, n

        ia = p(i1)
        i = i1 - 1
        id = 0

        do j = 1, i
          if ( ia < p(j) ) then
            id = id + 1
          end if
        end do

        is = id + is
        if ( id /= i * mod ( is, 2 ) ) then
          more = .true.
          exit
        end if

      end do

      if ( .not. more ) then
        p(1) = 0
        return
      end if

    end if

    m = mod ( is + 1, 2 ) * ( n + 1 )

    do j = 1, i

      if ( sign ( 1, p(j)-ia ) /= sign ( 1, p(j)-m ) ) then
        m = p(j)
        l = j
      end if

    end do

    p(l) = ia
    p(i1) = m
    even = .true.

  end if

  return
end subroutine

subroutine i4vec_indicator1 ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR1 sets an I4VEC to the indicator vector (1,...,N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, integer A(N), the array to be initialized.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)
  integer i

  do i = 1, n
    a(i) = i
  end do

  return
end subroutine

end module
