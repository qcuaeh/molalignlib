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
   integer :: perm_size  ! total permutation size 
   integer :: count  ! partial permutation count
   integer, allocatable :: subset(:)  ! indices of assigned entries in permutation
   integer, allocatable :: permut(:)  ! permut(i) = j means atom i -> atom j
   integer, allocatable :: backward(:)  ! backward(j) = i means atom i -> atom j
end type

interface operator(==)
   module procedure subperm_equality
end interface

contains

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

elemental function subperm_equality(left, right) result(equality)
   type(subperm_t), intent(in) :: left, right
   logical :: equality
   integer :: i

   if (left%count /= right%count) then
      equality = .false.
      return
   end if

!   do i = 1, left%count
!      if (left%subset(i) /= right%subset(i)) then
!         equality = .false.
!         return
!      end if
!   end do

   do i = 1, left%count
      if (left%permut(i) /= right%permut(i)) then
         equality = .false.
         return
      end if
   end do

   equality = .true.
end function

subroutine subperm_init(subperm, perm_size)
   type(subperm_t), intent(out) :: subperm
   integer, intent(in) :: perm_size

   allocate(subperm%subset(perm_size))
   allocate(subperm%permut(perm_size))
   subperm%perm_size = perm_size
   subperm%count = 0
end subroutine

subroutine subperm_add(subperm, i1, i2)
   ! Merge source assignment into target assignment
   type(subperm_t), intent(inout) :: subperm
   integer, intent(in) :: i1, i2
   ! Local variables
   integer :: n

block
   ! Check for merge conflicts
   integer :: i
   do i = 1, subperm%count
      if (subperm%subset(i) == i1) then
         error stop "Merge conflict in subset"
      end if
      if (subperm%permut(i) == i2) then
         error stop "Merge conflict in permut"
      end if
   end do
end block

   n = subperm%count + 1
   subperm%subset(n) = i1
   subperm%permut(n) = i2
   subperm%count = n
end subroutine

subroutine subperm_merge(subperm, other_subperm)
   ! Merge source assignment into target assignment
   type(subperm_t), intent(inout) :: subperm
   type(subperm_t), intent(in) :: other_subperm
   integer :: i, i1, i2

   do i = 1, other_subperm%count
      i1 = other_subperm%subset(i)
      i2 = other_subperm%permut(i)
      call subperm_add( subperm, i1, i2)
   end do
end subroutine

subroutine subperm_to_perm(subperm, perm, invperm)
   ! Merge source assignment into target assignment
   type(subperm_t), intent(in) :: subperm
   integer, dimension(:), allocatable, intent(out) :: perm, invperm
   ! Local variables
   integer :: i, i1, i2, j1, j2

   perm = identity_permutation(subperm%perm_size)
   invperm = identity_permutation(subperm%perm_size)

   do i = 1, subperm%count
      i1 = subperm%subset(i)
      i2 = subperm%permut(i)
      j1 = invperm(i2)
      j2 = perm(i1)
      perm(i1) = i2
      invperm(i2) = i1
      perm(j1) = j2
      invperm(j2) = j1
   end do
end subroutine

subroutine check_subperm(subperm)
   implicit none
   type(subperm_t), intent(in) :: subperm
   logical :: subset_seen(subperm%perm_size)
   logical :: permut_seen(subperm%perm_size)
   logical :: has_errors
   integer :: i, src_idx, tgt_idx
   
   has_errors = .false.
   subset_seen = .false.
   permut_seen = .false.
   
   ! Check if count is within valid bounds
   if (subperm%count < 0 .or. subperm%count > subperm%perm_size) then
      write(stderr, '(A,I0,A,I0,A)') 'Count ', subperm%count, &
         ' is out of range [0,', subperm%perm_size, ']'
      has_errors = .true.
   end if
   
   ! Check each pair in the partial permutation
   do i = 1, subperm%count
      src_idx = subperm%subset(i)
      tgt_idx = subperm%permut(i)
      
      ! Check if source index is out of bounds
      if (src_idx < 1 .or. src_idx > subperm%perm_size) then
         write(stderr, '(A,I0,A,I0,A,I0,A)') 'Source index ', src_idx, &
            ' at position ', i, ' is out of range [1,', subperm%perm_size, ']'
         has_errors = .true.
      else
         ! Only check for repetition if within bounds
         if (subset_seen(src_idx)) then
            write(stderr, '(A,I0,A)') 'Source index ', src_idx, &
               ' appears multiple times in subset'
            has_errors = .true.
         else
            subset_seen(src_idx) = .true.
         end if
      end if
      
      ! Check if target index is out of bounds
      if (tgt_idx < 1 .or. tgt_idx > subperm%perm_size) then
         write(stderr, '(A,I0,A,I0,A,I0,A)') 'Target index ', tgt_idx, &
            ' at position ', i, ' is out of range [1,', subperm%perm_size, ']'
         has_errors = .true.
      else
         ! Only check for repetition if within bounds
         if (permut_seen(tgt_idx)) then
            write(stderr, '(A,I0,A)') 'Target index ', tgt_idx, &
               ' appears multiple times in permutation'
            has_errors = .true.
         else
            permut_seen(tgt_idx) = .true.
         end if
      end if
   end do
   
   if (has_errors) then
      error stop 'Sub-permutation is not valid'
   else
      write(stderr, '(A,I0,A,I0,A)') 'Sub-permutation is valid (', &
         subperm%count, '/', subperm%perm_size, ' assignments)'
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
