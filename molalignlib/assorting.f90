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
use types
use flags
use bounds
use discrete
use sorting
use chemdata

implicit none

contains

subroutine assort_atoms(mol)
! Purpose: Group atoms by atomic numbers and types
   type(Molecule), intent(inout) :: mol

   integer :: i, j
   integer :: natomtype
   integer :: atomtypeidcs(mol%natom)
   logical :: remaining(mol%natom)
   integer, dimension(mol%natom) :: blkatomnums, blkatomtags
   integer, dimension(mol%natom) :: foreorder, backorder

   integer, allocatable :: atomnums(:), atomtags(:)
   real(wp), allocatable :: weights(:)

   ! Initialization

   natomtype = 0
   remaining = .true.
   atomnums = mol%get_atomnums()
   atomtags = mol%get_atomtags()
   weights = mol%get_weights()

   ! Create block list

   do i = 1, mol%natom

      if (remaining(i)) then

         natomtype = natomtype + 1
         atomtypeidcs(i) = natomtype
         blkatomnums(natomtype) = atomnums(i)
         blkatomtags(natomtype) = atomtags(i)
         remaining(i) = .false.

         do j = 1, mol%natom
            if (remaining(j)) then
               if (atomnums(i) == atomnums(j) .and. atomtags(i) == atomtags(j)) then
!                  if (weights(i) == weights(j)) then
                     atomtypeidcs(j) = natomtype
                     atomnums(j) = atomnums(i)
                     weights(j) = weights(i)
                     remaining(j) = .false.
!                  else
!                     ! Abort if there are inconsistent weights
!                     write (stderr, '(a)') 'Error: There are incosistent weights'
!                     stop
!                  end if
               end if
            end if
         end do

      end if

   end do

   ! Order parts by atom tag

   foreorder(:natomtype) = sorted_order(blkatomtags, natomtype)
   backorder(:natomtype) = inverse_mapping(foreorder(:natomtype))
   atomtypeidcs = backorder(atomtypeidcs)

   ! Order parts by atomic number

   foreorder(:natomtype) = sorted_order(blkatomnums, natomtype)
   backorder(:natomtype) = inverse_mapping(foreorder(:natomtype))
   atomtypeidcs = backorder(atomtypeidcs)

   call mol%set_atomtypeidcs(natomtype, atomtypeidcs) 

end subroutine

subroutine set_equiv_atoms(mol)
! Group atoms by MNA type at infinite level
   type(Molecule), intent(inout) :: mol

   integer i, nin, natomequiv
   integer, dimension(mol%natom) :: atomequividcs
   integer, dimension(mol%natom) :: intype, atomtypemap
   integer, dimension(mol%natom) :: foreorder, backorder

   ! Determine MNA types iteratively

   nin = mol%get_natomtype()
   intype = mol%get_atomtypeidcs()
   atomtypemap = [(i, i=1, nin)]

   do

      call getmnatypes(mol, nin, intype, natomequiv, atomequividcs, atomtypemap)

      if (all(atomequividcs == intype)) exit

      nin = natomequiv
      intype = atomequividcs

   end do

   foreorder(:natomequiv) = sorted_order(atomtypemap, natomequiv)
   backorder(:natomequiv) = inverse_mapping(foreorder(:natomequiv))
   atomequividcs = backorder(atomequividcs)

   call mol%set_atomequividcs(natomequiv, atomequividcs)

end subroutine

subroutine groupbytype(nelem, elements, types, ngroup, groupsize, groupid)
! Purpose: Categorize atoms by types
    integer, intent(in) :: nelem
    integer, intent(out) :: ngroup
    integer, dimension(:), intent(in) :: types, elements
    integer, dimension(:), intent(out) :: groupid, groupsize

    integer i, j
    integer, dimension(maxcoord) :: grouptype
    logical, dimension(maxcoord) :: remaining
    integer, dimension(maxcoord) :: foreorder, backorder

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

    foreorder(:ngroup) = sorted_order(grouptype, ngroup)
    backorder(:ngroup) = inverse_mapping(foreorder(:ngroup))
    grouptype(:ngroup) = grouptype(foreorder(:ngroup))
    groupsize(:ngroup) = groupsize(foreorder(:ngroup))
    groupid(:nelem) = backorder(groupid(:nelem))

!    print *, groupid(:ngroup)
!    print *, types(foreorder)

end subroutine

subroutine assort_neighbors(mol)
! Purpose: Categorize atoms by eqtypes
   type(Molecule), intent(inout) :: mol

   integer :: i, h
   integer :: nadjs(mol%natom), adjlists(maxcoord, mol%natom)
   integer :: nadjequivs(mol%natom), adjequivlenlists(maxcoord, mol%natom)
   integer :: adjeqvid(maxcoord), atomorder(maxcoord), atomequividcs(mol%natom)

   nadjs = mol%get_nadjs()
   adjlists = mol%get_adjlists()
   atomequividcs = mol%get_atomequividcs()

   do i = 1, mol%natom
      call groupbytype(nadjs(i), adjlists(:, i), atomequividcs, &
            nadjequivs(i), adjequivlenlists(:, i), adjeqvid)
      atomorder(:nadjs(i)) = sorted_order(adjeqvid, nadjs(i))
      adjlists(:nadjs(i), i) = adjlists(atomorder(:nadjs(i)), i)
   end do

   call mol%set_adjlists(nadjs, adjlists)
   call mol%set_adjequivlenlists(nadjequivs, adjequivlenlists)

end subroutine

subroutine getmnatypes(mol, nin, intype, nout, outype, atomtypemap)
   type(Molecule), intent(inout) :: mol
   integer, intent(in) :: nin
   integer, dimension(:), intent(in) :: intype
   integer, intent(out) :: nout
   integer, dimension(:), intent(out) :: outype, atomtypemap

   integer :: i, j
   integer :: nadjs(mol%natom)
   integer :: adjlists(maxcoord, mol%natom)
   integer :: parentype(mol%natom)
   logical :: untyped(mol%natom)

   nadjs = mol%get_nadjs()
   adjlists = mol%get_adjlists()

   nout = 0
   untyped(:) = .true.

   do i = 1, mol%natom
      if (untyped(i)) then
         nout = nout + 1
         outype(i) = nout
         parentype(nout) = intype(i)
         do j = i + 1, mol%natom
!               print '(a, x, i0, x, i0)', trim(elsym(intype0(i))), i, j
            if (untyped(j)) then
               if (intype(j) == intype(i)) then
                  if (same_adjacency(nin, intype, nadjs(i), adjlists(:, i), intype, &
                        nadjs(j), adjlists(:, j))) then
                     outype(j) = nout
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

   atomtypemap(:nout) = atomtypemap(parentype(:nout))

end subroutine

subroutine calcequivmat(mol0, mol1, nadjmna0, adjmnalen0, adjmnalist0, &
   nadjmna1, adjmnalen1, adjmnalist1, equivmat)
! Purpose: Calculate the maximum common MNA level for all atom cross assignments
   type(Molecule), intent(in) :: mol0, mol1

   integer, dimension(:, :), intent(out) :: nadjmna0, nadjmna1
   integer, dimension(:, :, :), intent(out) :: adjmnalen0, adjmnalen1
   integer, dimension(:, :, :), intent(out) :: adjmnalist0, adjmnalist1
   integer, dimension(:, :), intent(out) :: equivmat

   integer :: natom, natomtype
   integer, allocatable :: atomtypelenlist(:)
   integer, dimension(:), allocatable :: nadjs0, nadjs1
   integer, dimension(:, :), allocatable :: adjlists0, adjlists1

   integer :: h, i, j, offset, level, nin, nout
   integer, dimension(mol0%natom) :: intype0, intype1, outype0, outype1
   integer, dimension(maxcoord) :: indices, atomorder

   natom = mol0%get_natom()
   natomtype = mol0%get_natomtype()
   atomtypelenlist = mol0%get_atomtypelenlist()
   nadjs0 = mol0%get_sorted_nadjs()
   nadjs1 = mol1%get_sorted_nadjs()
   adjlists0 = mol0%get_sorted_adjlists()
   adjlists1 = mol1%get_sorted_adjlists()

   nin = natomtype
   intype0 = mol0%get_sorted_atomtypeidcs()
   intype1 = mol1%get_sorted_atomtypeidcs()
   level = 1

   do

      offset = 0
      do h = 1, natomtype
         do i = offset + 1, offset + atomtypelenlist(h)
            do j = offset + 1, offset + atomtypelenlist(h)
               if (intype0(i) == intype1(j)) then
                  equivmat(j, i) = level
               end if
            end do
         end do
         offset = offset + atomtypelenlist(h)
      end do

      do i = 1, natom
         call groupbytype(nadjs0(i), adjlists0(:, i), intype0, nadjmna0(i, level), adjmnalen0(:, i, level), indices)
         atomorder(:nadjs0(i)) = sorted_order(indices, nadjs0(i))
         adjmnalist0(:nadjs0(i), i, level) = adjlists0(atomorder(:nadjs0(i)), i)
      end do

      do i = 1, natom
         call groupbytype(nadjs1(i), adjlists1(:, i), intype1, nadjmna1(i, level), adjmnalen1(:, i, level), indices)
         atomorder(:nadjs1(i)) = sorted_order(indices, nadjs1(i))
         adjmnalist1(:nadjs1(i), i, level) = adjlists1(atomorder(:nadjs1(i)), i)
      end do

      call getmnacrosstypes(mol0%natom, nin, intype0, nadjs0, adjlists0, intype1, nadjs1, adjlists1, &
         nout, outype0, outype1)

      if (all(outype0 == intype0) .and. all((outype1 == intype1))) exit

      nin = nout
      intype0 = outype0
      intype1 = outype1
      level = level + 1

   end do

!    print *, natom, nout

end subroutine

subroutine getmnacrosstypes(natom, nin, intype0, nadjs0, adjlists0, intype1, nadjs1, adjlists1, &
      nout, outype0, outype1)
   integer, intent(in) :: natom, nin
   integer, dimension(:), intent(in) :: intype0, intype1, nadjs0, nadjs1
   integer, dimension(:, :), intent(in) :: adjlists0, adjlists1
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
                  if (same_adjacency(nin, intype0, nadjs0(i), adjlists0(:, i), intype0, nadjs0(j), adjlists0(:, j))) then
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
               if (same_adjacency(nin, intype0, nadjs0(archetype(i)), adjlists0(:, archetype(i)), intype1, nadjs1(j), &
                            adjlists1(:, j))) then
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
                  if (same_adjacency(nin, intype1, nadjs1(i), adjlists1(:, i), intype1, nadjs1(j), adjlists1(:, j))) then
                     outype1(j) = nout
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

end subroutine

function same_adjacency(natomtype, atomtype0, nadjs0, adjlists0, atomtype1, nadjs1, adjlists1) result(sameadj)
   integer, intent(in) :: natomtype, nadjs0, nadjs1
   integer, dimension(:), intent(in) :: adjlists0, adjlists1
   integer, dimension(:) :: atomtype0, atomtype1
   logical :: sameadj

   integer :: i0, i1
   integer, dimension(natomtype) :: n0, n1
!   real :: atoms0(3, maxcoord), atoms1(3, maxcoord)
!   integer :: typelist0(maxcoord, nin), typelist1(maxcoord, nin)

   sameadj = .true.

   if (nadjs0 /= nadjs1) then
      sameadj = .false.
      return
   end if

!   If coordination number is the same check if coordinated atoms are the same

   n0(:) = 0
   n1(:) = 0

   do i0 = 1, nadjs0
      n0(atomtype0(adjlists0(i0))) = n0(atomtype0(adjlists0(i0))) + 1
!       typelist0(n0(atomtype0(adjlists0(i0))), atomtype0(adjlists0(i0))) = i0
   end do

   do i1 = 1, nadjs1
      n1(atomtype1(adjlists1(i1))) = n1(atomtype1(adjlists1(i1))) + 1
!       typelist1(n1(atomtype1(adjlists1(i1))), atomtype1(adjlists1(i1))) = i1
   end do

   if (any(n0 /= n1)) then
      sameadj = .false.
      return
   end if

!   print *, typelist0(:nadjs0), '/', typelist1(:nadjs1)

end function

end module
