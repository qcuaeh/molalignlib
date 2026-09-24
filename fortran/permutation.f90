! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

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

! Partial atom permutation, built by appending assigned pairs.
!   atomset(1:atomset_size) lists the assigned atoms of molecule 1, in the
!   order they were added; atomperm(i1) is the partner of i1 in molecule 2.
!   Only atomperm(atomset(1:atomset_size)) is defined. Other entries are
!   stale or uninitialised and must never be read; to export a permutation,
!   scatter the defined entries into a zeroed array.
! Both arrays have fixed capacity n_atoms. Pairs are only ever appended and
! an atom is never re-added, so a subperm can be rolled back to an earlier
! state by restoring atomset_size, and copied cheaply with subperm_merge
! (O(assigned pairs)) instead of intrinsic assignment (O(n_atoms)).
type :: subperm_t
   integer(ik) :: atomset_size = 0
   integer(ik), allocatable :: atomset(:)
   integer(ik), allocatable :: atomperm(:)
end type

abstract interface
   integer function int_f(i)
      use parameters
      integer(ik), intent(in) :: i
   end function
end interface

interface init_array
   module procedure init_array_integer
end interface

contains

integer function identity(i)
   integer(ik), intent(in) :: i
   identity = i
end function

subroutine init_array_integer(array, n, f)
   integer(ik), dimension(:), allocatable, intent(out) :: array
   integer, intent(in) :: n
   procedure(int_f) :: f
   integer(ik) :: i

   allocate (array(n))

   do i = 1, n
      array(i) = f(i)
   end do
end subroutine

subroutine print_permutation(permutation)
   integer(ik), dimension(:), intent(in) :: permutation
   integer(ik) :: i

   write (stdout,'(I0)',advance='no') permutation(1)
   do i = 2, size(permutation)
      write (stdout,'(",",I0)',advance='no') permutation(i)
   end do
end subroutine

function inverse_permutation(permutation)
   integer(ik), dimension(:), intent(in) :: permutation
   integer(ik), dimension(:), allocatable :: inverse_permutation
   ! Local variables
   integer(ik) :: i

   allocate (inverse_permutation, mold=permutation)

   do i = 1, size(permutation)
      inverse_permutation(permutation(i)) = i
   end do
end function

logical(lk) function is_permutation(perm) result(isperm)
   implicit none
   integer(ik), dimension(:), intent(in) :: perm
   logical(lk) :: seen(size(perm))
   integer(ik) :: n, i

   n = size(perm)
   seen = .FALSE.
   isperm = .TRUE.

   do i = 1, n
      if (perm(i) < 1 .or. perm(i) > n) then
         isperm = .FALSE.
         return
      end if
      if (seen(perm(i))) then
         isperm = .FALSE.
         return
      end if
      seen(perm(i)) = .TRUE.
   end do
end function

subroutine subperm_init(subperm, perm_size)
! Allocate an empty subperm. Call it once per (unallocated) object; to reuse
! an initialised subperm just set atomset_size = 0. No array is filled:
! only entries listed in atomset are ever read.
   type(subperm_t), intent(inout) :: subperm
   integer(ik), intent(in) :: perm_size

   allocate (subperm%atomset(perm_size), subperm%atomperm(perm_size))
   subperm%atomset_size = 0
end subroutine

subroutine subperm_add(subperm, i1, i2)
   type(subperm_t), intent(inout) :: subperm
   integer(ik), intent(in) :: i1, i2
   integer(ik) :: n

   if (DEBUG_TESTS) then
   block
      integer(ik) :: i
      ! Check if i1 is already in atomset
      do i = 1, subperm%atomset_size
         if (subperm%atomset(i) == i1) then
            error stop 'Index i1 is already in atomset'
         end if
      end do
      ! Check if i2 is already assigned to something in atomset
      do i = 1, subperm%atomset_size
         if (subperm%atomperm(subperm%atomset(i)) == i2) then
            error stop 'Index i2 is already assigned'
         end if
      end do
   end block
   end if

   n = subperm%atomset_size + 1
   subperm%atomset(n) = i1
   subperm%atomperm(i1) = i2
   subperm%atomset_size = n
end subroutine

subroutine subperm_merge(subperm, other_subperm)
   type(subperm_t), intent(inout) :: subperm
   type(subperm_t), intent(in) :: other_subperm
   integer(ik) :: i, i1, i2

   do i = 1, other_subperm%atomset_size
      i1 = other_subperm%atomset(i)
      i2 = other_subperm%atomperm(i1)
      call subperm_add(subperm, i1, i2)
   end do
end subroutine

! Original at https://people.sc.fsu.edu/~jburkardt/f_src/atomset/atomset.f90
subroutine perm1_next3 ( n, p, more, rank )

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

  integer(ik) n
  integer(ik) i
  integer(ik) m2
  logical(lk) more
  integer(ik) n2
  integer(ik) p(n)
  integer(ik) q
  integer(ik) rank
  integer(ik) s
  integer(ik) t

  if ( .not. more ) then
    do i = 1, n
      p(i) = i
    end do
    more = .TRUE.
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
        more = .FALSE.
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
