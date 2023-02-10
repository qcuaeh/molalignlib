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
use kinds
use discrete
use sorting
use chemdata

implicit none

contains

subroutine grouptypes(natom, znums, types, nblk, blksz, blkid)
! Purpose: Group atoms by atomic numbers and types
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: znums, types
   integer, intent(out) :: nblk
   integer, dimension(:), intent(out) :: blksz, blkid

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
         blksz(nblk) = 1
         blkznum(nblk) = znums(i)
         blktype(nblk) = types(i)
         blkid(i) = nblk
         do j = i + 1, natom
            if (remaining(i)) then
               if (znums(j) == znums(i) .and. types(j) == types(i)) then
                  blkid(j) = nblk
                  blksz(nblk) = blksz(nblk) + 1
                  remaining(j) = .false.
               end if
            end if
         end do
      end if
   end do

! Order blocks by atomic type

   blkorder(:nblk) = order(blktype, nblk)
   blksz(:nblk) = blksz(blkorder(:nblk))
   blkorder(:nblk) = inverseperm(blkorder(:nblk))
   blkid = blkorder(blkid)

! Order blocks by atomic number

   blkorder(:nblk) = order(blkznum, nblk)
   blksz(:nblk) = blksz(blkorder(:nblk))
   blkorder(:nblk) = inverseperm(blkorder(:nblk))
   blkid = blkorder(blkid)

! Order blocks by block size

!   blkorder(:nblk) = order(blksz, nblk)
!   blksz(:nblk) = blksz(blkorder(:nblk))
!   blkorder(:nblk) = inverseperm(blkorder(:nblk))
!   blkid = blkorder(blkid)

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

subroutine getmnatypes(natom, nin, intype, nadj, adjlist, nout, outype, outsize, parentype)
    integer, intent(in) :: natom, nin
    integer, dimension(:), intent(in) :: intype, nadj
    integer, dimension(:, :), intent(in) :: adjlist
    integer, intent(out) :: nout
    integer, dimension(:), intent(out) :: outype, outsize, parentype

    integer :: i, j
    logical :: untyped(natom)

    nout = 0
    untyped(:) = .true.

    do i = 1, natom
        if (untyped(i)) then
            nout = nout + 1
            outype(i) = nout
            outsize(nout) = 1
            parentype(nout) = intype(i)
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

subroutine groupequiv(natom, nblk, blkid, nadj, adjlist, neqv, eqvsz, eqvid)
! Group atoms by MNA at infinite lavel
    integer, intent(in) :: natom, nblk
    integer, dimension(:), intent(in) :: nadj, blkid
    integer, dimension(:, :), intent(in) :: adjlist
    integer, dimension(:), intent(out) :: eqvid, eqvsz

    integer i, nin, neqv
    integer, dimension(natom) :: intype, grouporder, parentype, basetype

    nin = nblk
    intype = blkid
    basetype = [(i, i=1, natom)]

    do

        call getmnatypes(natom, nin, intype, nadj, adjlist, neqv, eqvid, eqvsz, parentype)
        basetype(:neqv) = basetype(parentype(:neqv))

        if (all(eqvid == intype)) exit

        nin = neqv
        intype = eqvid

    end do

    grouporder(:neqv) = order(basetype, neqv)
    eqvsz(:neqv) = eqvsz(grouporder(:neqv))
    grouporder(:neqv) = inverseperm(grouporder(:neqv))
    eqvid = grouporder(eqvid)

!    do i = 1, natom
!        print *, i, elsym(znum(i)), eqvid(i)
!    end do

end subroutine

end module
