! MolAlignLib
! Copyright (C) 2022 José M. Vásquez

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module permutation

implicit none

contains

! Get an identity permutation
function identity_permutation(n) result(mapping)
   integer, intent(in) :: n
   ! Local variables
   integer :: i
   integer :: mapping(n)

   do i = 1, n
      mapping(i) = i
   end do

end function

! Get the inverse mapping of mapping
function inverse_permutation(mapping)
   integer, intent(in) :: mapping(:)
   ! Local variables
   integer :: i
   integer, allocatable :: inverse_permutation(:)

   allocate (inverse_permutation(size(mapping)))

   do i = 1, size(mapping)
      inverse_permutation(mapping(i)) = i
   end do

end function

logical function is_permutation(arr)
   implicit none
   integer, intent(in) :: arr(:)
   integer :: i, N
   logical :: seen(size(arr))
   
   N = size(arr)
   seen = .false.
   is_permutation = .true.
   
   do i = 1, N
      if (arr(i) < 1 .or. arr(i) > N .or. seen(arr(i))) then
         is_permutation = .false.
         return
      end if
      seen(arr(i)) = .true.
   end do

end function

subroutine perm1_next3 ( n, p, more, rank )
!This subroutine was obtained from:
!https://people.sc.fsu.edu/~jburkardt/f_src/subset/subset.f90

!*****************************************************************************80
!
!! perm1_next3() computes permutations of (1,...,N).
!
!  Discussion:
!
!    The routine is initialized by calling with MORE = TRUE, in which case
!    it returns the identity permutation.
!
!    If the routine is called with MORE = FALSE, then the successor of the
!    input permutation is computed.
!
!    Trotter's algorithm is used.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 November 2018
!
!  Author:
!
!    Original Fortran77 version by Hale Trotter,
!    This version by John Burkardt
!
!  Reference:
!
!    Hale Trotter,
!    Algorithm 115:
!    PERM,
!    Communications of the Association for Computing Machinery,
!    Volume 5, 1962, pages 434-435.
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer P(N), the permutation, in standard
!    index form.  If MORE is TRUE, then P is assumed to contain the
!    "previous" permutation, and on P(I) is the value
!    of the I-th object under the next permutation.
!    Otherwise, P will be set to the "first" permutation.
!
!    Input/output, logical MORE.
!    Set MORE = FALSE before first calling this routine.
!    MORE will be reset to TRUE and a permutation will be returned.
!    Each new call produces a new permutation until MORE is returned FALSE.
!
!    Input/output, integer RANK, the rank of the current permutation.
!
  integer :: n
  integer :: p(n)
  logical :: more
  integer :: rank

  integer :: i
  integer :: m2
  integer :: n2
  integer :: q
  integer :: s
  integer :: t

  if ( .not. more ) then

    do i = 1, n
      p(i) = i
    end do

    more = .true.
    rank = 1

  else

    n2 = n
    m2 = rank
    s = n

    do

      q = mod ( m2, n2 )
      t = mod ( m2, 2 * n2 )

      if ( q /= 0 ) then
        exit
      end if

      if ( t == 0 ) then
        s = s - 1
      end if

      m2 = m2 / n2
      n2 = n2 - 1

      if ( n2 == 0 ) then
        do i = 1, n
          p(i) = i
        end do
        more = .false.
        rank = 1
        exit
      end if

    end do

    if ( n2 /= 0 ) then

      if ( q == t ) then
        s = s - q
      else
        s = s + q - n2
      end if
!
!  Swap.
!
      t      = p(s)
      p(s)   = p(s+1)
      p(s+1) = t

      rank = rank + 1

    end if

  end if

  return
end subroutine

end module
