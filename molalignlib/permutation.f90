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
use parameters
implicit none

! Derived type to store permutation subsets
type :: subperm_t
   integer :: size                    ! number of assigned entries
   integer, allocatable :: subset(:)  ! indices of assigned entries in permutation
   integer, allocatable :: forward(:)  ! forward(i) = j means atom i -> atom j
   integer, allocatable :: backward(:)  ! backward(j) = i means atom i -> atom j
end type

interface operator(==)
   module procedure subperm_equality
end interface

contains

elemental function subperm_equality(left, right) result(equality)
   type(subperm_t), intent(in) :: left, right
   logical :: equality
   equality = all(left%forward == right%forward)
end function

subroutine subperm_init(subperm, perm_size)
   type(subperm_t), intent(out) :: subperm
   integer, intent(in) :: perm_size

   allocate(subperm%subset(perm_size))
   allocate(subperm%forward(perm_size))
   subperm%size = 0
   subperm%forward = identity_permutation(perm_size)
   subperm%backward = identity_permutation(perm_size)
end subroutine

subroutine subperm_add(subperm, i1, j1)
   ! Merge source assignment into target assignment
   type(subperm_t), intent(inout) :: subperm
   integer, intent(in) :: i1, j1
   ! Local variables
   integer :: n, i2, j2

   ! Check for conflicts in existing assignments
   if (any(subperm%subset(1:subperm%size) == i1)) then
      write(stderr, '(A,I0,A)') "ERROR: Attempting to overwrite permutation at position ", i1
      error stop "Assignment merge conflict"
   end if

   n = subperm%size + 1
   i2 = subperm%backward(j1)
   j2 = subperm%forward(i1)

   subperm%size = n
   subperm%subset(n) = i1
   subperm%forward(i1) = j1
   subperm%backward(j1) = i1
   subperm%forward(i2) = j2
   subperm%backward(j2) = i2
end subroutine

subroutine subperm_merge(subperm, other_subperm)
   ! Merge source assignment into target assignment
   type(subperm_t), intent(inout) :: subperm
   type(subperm_t), intent(in) :: other_subperm
   integer :: i, i1, j1

   do i = 1, other_subperm%size
      i1 = other_subperm%subset(i)
      j1 = other_subperm%forward(i1)
      call subperm_add( subperm, i1, j1)
   end do
end subroutine

! Get an identity permutation
function identity_permutation(n) result(perm)
   integer, intent(in) :: n
   integer, allocatable :: perm(:)
   ! Local variables
   integer :: i

   allocate (perm(n))

   do i = 1, n
      perm(i) = i
   end do
end function

! Get the inverse permutation of perm
function inverse_permutation(perm) result(invperm)
   integer, dimension(:), intent(in) :: perm
   integer, dimension(:), allocatable :: invperm
   ! Local variables
   integer :: n, i

   n = size(perm)
   allocate (invperm(n))

   do i = 1, n
      invperm(perm(i)) = i
   end do
end function

logical function is_permutation(perm) result(isperm)
   implicit none
   integer, dimension(:), intent(in) :: perm
   logical :: seen(size(perm))
   integer :: n, i

   n = size(perm)
   seen = .false.
   isperm = .true.

   do i = 1, n
      if (perm(i) < 1 .or. perm(i) > n) then
         isperm = .false.
         return
      end if
      if (seen(perm(i))) then
         isperm = .false.
         return
      end if
      seen(perm(i)) = .true.
   end do
end function

subroutine check_permutation(perm)
   implicit none
   integer, dimension(:), intent(in) :: perm
   logical :: seen(size(perm))
   logical :: has_errors
   integer :: n, i
   
   n = size(perm)
   seen = .false.
   has_errors = .false.
   
   do i = 1, n
      ! Check if out of bounds
      if (perm(i) < 1 .or. perm(i) > n) then
         write(stderr, '(A,I0,A,I0,A)') 'Index ', perm(i), ' at position ', i, ' is out of range'
         has_errors = .true.
      else
         ! Only check for repetition if within bounds
         if (seen(perm(i))) then
            write(stderr, '(A,I0,A)') 'Index ', perm(i), ' appears multiple times in permutation'
            has_errors = .true.
         else
            seen(perm(i)) = .true.
         end if
      end if
   end do
   
   if (has_errors) then
      error stop 'Permutation is not valid'
   else
      write(stderr, '(A)') 'Permutation is valid'
   end if
end subroutine

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
