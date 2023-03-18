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

module biasing
use kinds
use flags
use bounds
use sorting

implicit none

abstract interface
   subroutine bias_func_proc(natom, nblk, blklen, nadj0, adjlist0, nadj1, adjlist1, coords0, coords1, biasmat)
      use kinds
      integer, intent(in) :: natom, nblk
      integer, dimension(:), intent(in) :: blklen
      integer, dimension(:), intent(in) :: nadj0, nadj1
      integer, dimension(:, :), intent(in) :: adjlist0, adjlist1
      real(wp), dimension(:, :), intent(in) :: coords0, coords1
      real(wp), dimension(:, :), intent(out) :: biasmat
   end subroutine
end interface

real(wp) :: bias_tol
real(wp) :: bias_scale
real(wp) :: bias_ratio
procedure(bias_func_proc), pointer :: bias_func

contains

subroutine nocrossbias(natom, nblk, blklen, nadj0, adjlist0, nadj1, adjlist1, coords0, coords1, biasmat)
! Purpose: Set biases from sorted neighbors' distances equivalence

   integer, intent(in) :: natom, nblk
   integer, dimension(:), intent(in) :: blklen
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   integer, dimension(:), intent(in) :: nadj0, nadj1
   integer, dimension(:, :), intent(in) :: adjlist0, adjlist1
   real(wp), dimension(:, :), intent(out) :: biasmat

   integer :: h, i, j, offset

   offset = 0

   do h = 1, nblk
      do i = offset + 1, offset + blklen(h)
         do j = offset + 1, offset + blklen(h)
            biasmat(i, j) = 0
         end do
      end do
      offset = offset + blklen(h)
   end do

end subroutine

subroutine sndcrossbias(natom, nblk, blklen, nadj0, adjlist0, nadj1, adjlist1, coords0, coords1, biasmat)
! Purpose: Set biases from sorted neighbors' distances equivalence

   integer, intent(in) :: natom, nblk
   integer, dimension(:), intent(in) :: blklen
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   integer, dimension(:), intent(in) :: nadj0, nadj1
   integer, dimension(:, :), intent(in) :: adjlist0, adjlist1
   real(wp), dimension(:, :), intent(out) :: biasmat

   integer :: h, i, j, offset
   real(wp), allocatable :: d0(:, :), d1(:, :)

   allocate(d0(natom, natom), d1(natom, natom))

   do i = 1, natom
      offset = 0
      do h = 1, nblk
         do j = offset + 1, offset + blklen(h)
            d0(j, i) = sqrt(sum((coords0(:, j) - coords0(:, i))**2))
         end do
         call sort(d0(:, i), offset + 1, offset + blklen(h))
         offset = offset + blklen(h)
      end do
   end do

   do i = 1, natom
      offset = 0
      do h = 1, nblk
         do j = offset + 1, offset + blklen(h)
            d1(j, i) = sqrt(sum((coords1(:, j) - coords1(:, i))**2))
         end do
         call sort(d1(:, i), offset + 1, offset + blklen(h))
         offset = offset + blklen(h)
      end do
   end do

   offset = 0

   do h = 1, nblk
      do i = offset + 1, offset + blklen(h)
         do j = offset + 1, offset + blklen(h)
            if (all(abs(d1(:, j) - d0(:, i)) < bias_tol)) then
               biasmat(i, j) = 0
            else
               biasmat(i, j) = bias_scale**2
            end if
         end do
      end do
      offset = offset + blklen(h)
   end do

end subroutine

subroutine mnacrossbias(natom, nblk, blklen, nadj0, adjlist0, nadj1, adjlist1, coords0, coords1, biasmat)
! Purpose: Set biases from sorted distances to neighbors equivalence

   integer, intent(in) :: natom, nblk
   integer, dimension(:), intent(in) :: blklen
   integer, dimension(:), intent(in) :: nadj0, nadj1
   integer, dimension(:, :), intent(in) :: adjlist0, adjlist1
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   real(wp), dimension(:, :), intent(out) :: biasmat

   integer :: h, i, j, offset, lev, nin, nout
   integer, dimension(natom) :: intype0, intype1, outype0, outype1

   lev = 0
   nin = nblk

   offset = 0
   do h = 1, nblk
      intype0(offset+1:offset+blklen(h)) = h
      intype1(offset+1:offset+blklen(h)) = h
      offset = offset + blklen(h)
   end do

   offset = 0
   do h = 1, nblk
      do i = offset + 1, offset + blklen(h)
         do j = offset + 1, offset + blklen(h)
            biasmat(i, j) = 0
         end do
      end do
      offset = offset + blklen(h)
   end do

   do

      lev = lev + 1

      if (maxlvl_flag) then
         if (lev > maxlevel) exit
      end if

      call getmnacrosstypes(natom, nin, intype0, nadj0, adjlist0, intype1, nadj1, adjlist1, &
         nout, outype0, outype1)

      offset = 0

      do h = 1, nblk
         do i = offset + 1, offset + blklen(h)
            do j = offset + 1, offset + blklen(h)
               if (lev >= 1 .and. outype0(i) /= outype1(j)) then
                  biasmat(i, j) = biasmat(i, j) + bias_scale**2*bias_ratio**(lev - 1)
               end if
            end do
         end do
         offset = offset + blklen(h)
      end do

      if (all(outype0 == intype0) .and. all((outype1 == intype1))) exit

      nin = nout
      intype0 = outype0
      intype1 = outype1

   end do

!    print *, natom, nout

end subroutine

subroutine getmnacrosstypes(natom, nin, intype0, nadj0, adjlist0, intype1, nadj1, adjlist1, &
      nout, outype0, outype1)
   integer, intent(in) :: natom, nin
   integer, dimension(:), intent(in) :: intype0, intype1, nadj0, nadj1
   integer, dimension(:, :), intent(in) :: adjlist0, adjlist1
   integer, intent(out) :: nout
   integer, dimension(:), intent(out) :: outype0, outype1

   integer i, j
   integer archetype(natom)
   logical untyped(natom)

   nout = 0
   untyped(:) = .true.

   do i = 1, natom
      if (untyped(i)) then
         nout = nout + 1
         outype0(i) = nout
         archetype(nout) = i
         do j = i + 1, natom
            if (untyped(j)) then
               if (intype0(j) == intype0(i)) then
                  if (sameadjacency(nin, intype0, nadj0(i), adjlist0(:, i), intype0, nadj0(j), adjlist0(:, j))) then
                     outype0(j) = nout
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

   untyped(:) = .true.

   do i = 1, nout
      do j = 1, natom
         if (untyped(j)) then
            if (intype1(j) == intype0(archetype(i))) then
               if (sameadjacency(nin, intype0, nadj0(archetype(i)), adjlist0(:, archetype(i)), intype1, nadj1(j), &
                            adjlist1(:, j))) then
                  outype1(j) = i
                  untyped(j) = .false.
               end if
            end if
         end if
      end do
   end do

   do i = 1, natom
      if (untyped(i)) then
         nout = nout + 1
         outype1(i) = nout
         do j = i + 1, natom
            if (untyped(j)) then
               if (intype1(j) == intype1(i)) then
                  if (sameadjacency(nin, intype1, nadj1(i), adjlist1(:, i), intype1, nadj1(j), adjlist1(:, j))) then
                     outype1(j) = nout
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

end subroutine

function sameadjacency(ntype, atomtype0, nadj0, adjlist0, atomtype1, nadj1, adjlist1)
   integer, intent(in) :: ntype, nadj0, nadj1
   integer, dimension(:), intent(in) :: adjlist0, adjlist1
   integer, dimension(:) :: atomtype0, atomtype1
   logical :: sameadjacency

   integer i0, i1
   integer n0(ntype), n1(ntype)
!   real atoms0(3, maxcoord), atoms1(3, maxcoord)
!   integer typelist0(maxcoord, nin), typelist1(maxcoord, nin)

   sameadjacency = .true.

   if (nadj0 /= nadj1) then
      sameadjacency = .false.
      return
   end if

!   If coordination number is the same check if coordinated atoms are the same

   n0(:) = 0
   n1(:) = 0

   do i0 = 1, nadj0
      n0(atomtype0(adjlist0(i0))) = n0(atomtype0(adjlist0(i0))) + 1
!       typelist0(n0(atomtype0(adjlist0(i0))), atomtype0(adjlist0(i0))) = i0
   end do

   do i1 = 1, nadj1
      n1(atomtype1(adjlist1(i1))) = n1(atomtype1(adjlist1(i1))) + 1
!       typelist1(n1(atomtype1(adjlist1(i1))), atomtype1(adjlist1(i1))) = i1
   end do

   if (any(n0 /= n1)) then
      sameadjacency = .false.
      return
   end if

!   print *, typelist0(:nadj0), '/', typelist1(:nadj1)

end function

end module
