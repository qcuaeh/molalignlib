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

module assorting
use stdio
use kinds
use flags
use bounds
use discrete
use sorting
use chemdata

implicit none

contains

subroutine groupatoms(natom, znums, ztypes, weights, nblk, blklen, blkidx)
! Purpose: Group atoms by atomic numbers and types
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: znums
   integer, dimension(:), intent(in) :: ztypes
   real(wp), dimension(:), intent(in) :: weights
   integer, intent(out) :: nblk
   integer, dimension(:), intent(out) :: blklen
   integer, dimension(:), intent(out) :: blkidx

   integer :: i, j
   logical :: remaining(natom)
   integer :: blkznum(natom)
   integer :: blktype(natom)
   integer :: blkorder(natom)

   ! Initialization

   nblk = 0
   remaining = .true.

   ! Create block list

   do i = 1, natom
      if (remaining(i)) then
         nblk = nblk + 1
         blklen(nblk) = 1
         blkznum(nblk) = znums(i)
         blktype(nblk) = ztypes(i)
         blkidx(i) = nblk
         do j = i + 1, natom
            if (remaining(i)) then
               if (znums(j) == znums(i) .and. ztypes(j) == ztypes(i)) then
                  if (weights(j) == weights(i)) then
                     blkidx(j) = nblk
                     blklen(nblk) = blklen(nblk) + 1
                     remaining(j) = .false.
                  else
                     ! Abort if there are inconsistent weights
                     write (stderr, '(a)') 'Error: There are incosistent weights'
                     stop
                  end if
               end if
            end if
         end do
      end if
   end do

   ! Order blocks by atomic type

   blkorder(:nblk) = sortorder(blktype, nblk)
   blklen(:nblk) = blklen(blkorder(:nblk))
   blkorder(:nblk) = inverseperm(blkorder(:nblk))
   blkidx = blkorder(blkidx)

   ! Order blocks by atomic number

   blkorder(:nblk) = sortorder(blkznum, nblk)
   blklen(:nblk) = blklen(blkorder(:nblk))
   blkorder(:nblk) = inverseperm(blkorder(:nblk))
   blkidx = blkorder(blkidx)

   ! Order blocks by block size

!   blkorder(:nblk) = sortorder(blklen, nblk)
!   blklen(:nblk) = blklen(blkorder(:nblk))
!   blkorder(:nblk) = inverseperm(blkorder(:nblk))
!   blkidx = blkorder(blkidx)

end subroutine

subroutine groupbytype(nelem, elements, types, groupid, ngroup, groupsize)
! Purpose: Categorize atoms by types
    integer, intent(in) :: nelem
    integer, intent(out) :: ngroup
    integer, dimension(:), intent(in) :: types, elements
    integer, dimension(:), intent(out) :: groupid, groupsize

    integer i, j
    integer, dimension(maxcoord) :: grouptype, grouporder
    logical, dimension(maxcoord) :: remaining

! Initialization

    ngroup = 0
    remaining = .true.

! Create group lists

    do i = 1, nelem
        if (remaining(i)) then
            ngroup = ngroup + 1
            groupsize(ngroup) = 1
            grouptype(ngroup) = types(elements(i))
            groupid(i) = ngroup
            do j = i + 1, nelem
                if (remaining(i)) then
                    if (types(elements(j)) == types(elements(i))) then
                        groupid(j) = ngroup
                        groupsize(ngroup) = groupsize(ngroup) + 1
                        remaining(j) = .false.
                    end if
                end if
            end do
        end if
    end do

! Order groups by category type 

    grouporder(:ngroup) = sortorder(grouptype, ngroup)
    grouptype(:ngroup) = grouptype(grouporder(:ngroup))
    groupsize(:ngroup) = groupsize(grouporder(:ngroup))
    grouporder(:ngroup) = inverseperm(grouporder(:ngroup))
    groupid(:nelem) = grouporder(groupid(:nelem))

!    print *, groupid(:ngroup)
!    print *, typess(grouporder)

end subroutine

subroutine groupneighbors(natom, neqv, eqvlen, nadj, adjlist, nadjeqv, adjeqvlen)
! Purpose: Categorize atoms by eqtypes
   integer, intent(in) :: natom, neqv
   integer, dimension(:), intent(in) :: nadj, eqvlen
   integer, dimension(:, :), intent(inout) :: adjlist
   integer, dimension(:), intent(out) :: nadjeqv
   integer, dimension(:, :), intent(out) :: adjeqvlen

   integer i, h, offset
   integer :: eqvidx(natom)
   integer, dimension(maxcoord) :: adjeqvid, atomorder

   ! assing equivalence group

   offset = 0
   do h = 1, neqv
      eqvidx(offset+1:offset+eqvlen(h)) = h
      offset = offset + eqvlen(h)
   end do

   do i = 1, natom
      call groupbytype(nadj(i), adjlist(:, i), eqvidx, adjeqvid, nadjeqv(i), adjeqvlen(:, i))
      atomorder(:nadj(i)) = sortorder(adjeqvid, nadj(i))
      adjlist(:nadj(i), i) = adjlist(atomorder(:nadj(i)), i)
   end do

!    do i = 1, natom
!        print '(a, 1x, i0, 1x, a)', '--------------', i, '--------------'
!        offset = 0
!        do j = 1, nadjeqv(i)
!            print *, adjlist(offset+1:offset+adjeqvlen(j, i), i)
!            offset = offset + adjeqvlen(j, i)
!         end do
!    end do

end subroutine

subroutine getmnatypes(natom, nin, intype, nadj, adjlist, nout, outype, outsize, uptype)
   integer, intent(in) :: natom, nin
   integer, dimension(:), intent(in) :: intype, nadj
   integer, dimension(:, :), intent(in) :: adjlist
   integer, intent(out) :: nout
   integer, dimension(:), intent(out) :: outype, outsize, uptype

   integer :: i, j
   logical :: untyped(natom)

   nout = 0
   untyped(:) = .true.

   do i = 1, natom
      if (untyped(i)) then
         nout = nout + 1
         outype(i) = nout
         outsize(nout) = 1
         uptype(nout) = intype(i)
         do j = i + 1, natom
!               print '(a, x, i0, x, i0)', trim(elsym(intype0(i))), i, j
            if (untyped(j)) then
               if (intype(j) == intype(i)) then
                  if (sameadjacency(nin, intype, nadj(i), adjlist(:, i), intype, nadj(j), adjlist(:, j))) then
                     outype(j) = nout
                     outsize(nout) = outsize(nout) + 1
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

end subroutine

subroutine groupequivatoms(natom, nblk, blkidx, nadj, adjlist, neqv, eqvlen, eqvidx)
! Group atoms by MNA type at infinite level
   integer, intent(in) :: natom, nblk
   integer, dimension(:), intent(in) :: nadj
   integer, dimension(:), intent(in) :: blkidx
   integer, dimension(:, :), intent(in) :: adjlist
   integer, dimension(:), intent(out) :: eqvidx
   integer, dimension(:), intent(out) :: eqvlen

   integer i, nin, neqv
   integer, dimension(natom) :: intype, grouporder, uptype, basetype

   ! Determine MNA types iteratively

   nin = nblk
   intype = blkidx
   basetype = [(i, i=1, natom)]

   do

      call getmnatypes(natom, nin, intype, nadj, adjlist, neqv, eqvidx, eqvlen, uptype)
      basetype(:neqv) = basetype(uptype(:neqv))

      if (all(eqvidx == intype)) exit

      nin = neqv
      intype = eqvidx

   end do

   grouporder(:neqv) = sortorder(basetype, neqv)
   eqvlen(:neqv) = eqvlen(grouporder(:neqv))
   grouporder(:neqv) = inverseperm(grouporder(:neqv))
   eqvidx = grouporder(eqvidx)

!    do i = 1, natom
!        print *, i, elsym(znum(i)), eqvidx(i)
!    end do

end subroutine

subroutine calcequivmat(natom, nblk, blklen, nadj0, adjlist0, nadjmna0, adjmnalen0, adjmnalist0, &
   nadj1, adjlist1, nadjmna1, adjmnalen1, adjmnalist1, equivmat)
! Purpose: Calculate the maximum common MNA level for all atom cross assignments
   integer, intent(in) :: natom, nblk
   integer, dimension(:), intent(in) :: blklen
   integer, dimension(:), intent(in) :: nadj0, nadj1
   integer, dimension(:, :), intent(inout) :: adjlist0, adjlist1
   integer, dimension(:, :), intent(out) :: nadjmna0, nadjmna1
   integer, dimension(:, :, :), intent(out) :: adjmnalen0, adjmnalen1
   integer, dimension(:, :, :), intent(out) :: adjmnalist0, adjmnalist1
   integer, dimension(:, :), intent(out) :: equivmat

   integer :: h, i, j, offset, level, nin, nout
   integer, dimension(natom) :: intype0, intype1, outype0, outype1
   integer, dimension(maxcoord) :: indices, atomorder

   level = 1
   nin = nblk

   offset = 0
   do h = 1, nblk
      intype0(offset+1:offset+blklen(h)) = h
      intype1(offset+1:offset+blklen(h)) = h
      offset = offset + blklen(h)
   end do

   do

      offset = 0
      do h = 1, nblk
         do i = offset + 1, offset + blklen(h)
            do j = offset + 1, offset + blklen(h)
               if (intype0(i) == intype1(j)) then
                  equivmat(j, i) = level
               end if
            end do
         end do
         offset = offset + blklen(h)
      end do

      do i = 1, natom
         call groupbytype(nadj0(i), adjlist0(:, i), intype0, indices, nadjmna0(i, level), adjmnalen0(:, i, level))
         atomorder(:nadj0(i)) = sortorder(indices, nadj0(i))
         adjmnalist0(:nadj0(i), i, level) = adjlist0(atomorder(:nadj0(i)), i)
      end do

      do i = 1, natom
         call groupbytype(nadj1(i), adjlist1(:, i), intype1, indices, nadjmna1(i, level), adjmnalen1(:, i, level))
         atomorder(:nadj1(i)) = sortorder(indices, nadj1(i))
         adjmnalist1(:nadj1(i), i, level) = adjlist1(atomorder(:nadj1(i)), i)
      end do

      call getmnacrosstypes(natom, nin, intype0, nadj0, adjlist0, intype1, nadj1, adjlist1, &
         nout, outype0, outype1)

      if (all(outype0 == intype0) .and. all((outype1 == intype1))) exit

      nin = nout
      intype0 = outype0
      intype1 = outype1
      level = level + 1

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

   integer :: i0, i1
   integer, dimension(ntype) :: n0, n1
!   real :: atoms0(3, maxcoord), atoms1(3, maxcoord)
!   integer :: typelist0(maxcoord, nin), typelist1(maxcoord, nin)

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
