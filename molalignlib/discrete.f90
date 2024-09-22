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

module discrete

implicit none

contains

! Get an identity permutation
function identity_mapping(n) result(mapping)
   integer, intent(in) :: n
   ! Local variables
   integer :: i
   integer :: mapping(n)

   do i = 1, n
      mapping(i) = i
   end do

end function

! Get the inverse mapping of mapping
function inverse_mapping(mapping)
   integer, intent(in) :: mapping(:)
   ! Local variables
   integer :: i
   integer, allocatable :: inverse_mapping(:)

   allocate(inverse_mapping(size(mapping)))

   do i = 1, size(mapping)
      inverse_mapping(mapping(i)) = i
   end do

end function

function intersection(list1, list2, hash_size)
   integer, intent(in) :: hash_size
   integer, intent(in) :: list1(:), list2(:)
   ! Local variables
   integer :: i, n
   logical, allocatable :: hash_table(:)
   integer, allocatable :: intersection(:)

   allocate(hash_table(hash_size))

   hash_table = .false.

   ! Add elements in list1 to hash table
   do i = 1, size(list1)
      hash_table(list1(i)) = .true.
   end do

   ! Count intersection elements
   n = 0
   do i = 1, size(list2)
      if (hash_table(list2(i))) then
         n = n + 1
      end if
   end do

   allocate(intersection(n))

   ! Store intersection elements
   n = 0
   do i = 1, size(list2)
      if (hash_table(list2(i))) then
         n = n + 1
         intersection(n) = list2(i)
      end if
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
