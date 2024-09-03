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
   type(cMol), intent(inout) :: mol

   integer :: i, j
   integer :: neltype
   integer :: atomeltypes(mol%natom)
   logical :: remaining(mol%natom)
   integer, dimension(mol%natom) :: foreorder, backorder
   integer, dimension(mol%natom) :: blockelnums, blocklabels

   integer, allocatable :: atomelnums(:), atomlabels(:)
!   real(wp), allocatable :: weights(:)

   ! Initialization

   neltype = 0
   remaining = .true.
   atomelnums = mol%get_atomelnums()
   atomlabels = mol%get_atomlabels()
!   weights = mol%get_atomweights()

   ! Create block list

   do i = 1, mol%natom

      if (remaining(i)) then

         neltype = neltype + 1
         atomeltypes(i) = neltype
         blockelnums(neltype) = atomelnums(i)
         blocklabels(neltype) = atomlabels(i)
         remaining(i) = .false.

         do j = 1, mol%natom
            if (remaining(j)) then
               if (atomelnums(i) == atomelnums(j) .and. atomlabels(i) == atomlabels(j)) then
!                  if (weights(i) == weights(j)) then
                     remaining(j) = .false.
                     atomeltypes(j) = neltype
                     atomelnums(j) = atomelnums(i)
!                     weights(j) = weights(i)
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

   foreorder(:neltype) = sorted_order(blocklabels, neltype)
   backorder(:neltype) = inverse_mapping(foreorder(:neltype))
   atomeltypes = backorder(atomeltypes)

   ! Order parts by atomic number

   foreorder(:neltype) = sorted_order(blockelnums, neltype)
   backorder(:neltype) = inverse_mapping(foreorder(:neltype))
   atomeltypes = backorder(atomeltypes)

   call mol%set_eltypes(neltype, atomeltypes) 

end subroutine

subroutine set_equiv_atoms(mol)
! Group atoms by MNA type at infinite level
   type(cMol), intent(inout) :: mol

   integer :: i, natom, nin, nout
   integer, allocatable, dimension(:) :: intype, outype, typemap
   integer, allocatable, dimension(:) :: foreorder, backorder

   allocate(outype(size(mol%atoms)))
   allocate(typemap(size(mol%atoms)))
   allocate(foreorder(size(mol%atoms)))
   allocate(backorder(size(mol%atoms)))

   ! Initialization

   nin = mol%get_neltype()
   intype = mol%get_atomeltypes()
   typemap(:nin) = [(i, i=1, nin)]

   ! Determine MNA types iteratively

   do

      call getmnatypes(mol, nin, nout, intype, outype, typemap)

      if (all(outype == intype)) exit

      nin = nout
      intype = outype

   end do

   foreorder(:nout) = sorted_order(typemap, nout)
   backorder(:nout) = inverse_mapping(foreorder(:nout))
   outype = backorder(outype)

   call mol%set_mnatypes(nout, outype)

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
   type(cMol), intent(inout) :: mol

   integer :: i, h
   integer :: nadjs(mol%natom), adjlists(maxcoord, mol%natom)
   integer :: natomneimnatypes(mol%natom), atomneimnatypepartlens(maxcoord, mol%natom)
   integer :: adjeqvid(maxcoord), atomorder(maxcoord), atommnatypes(mol%natom)

   nadjs = mol%get_nadjs()
   adjlists = mol%get_adjlists()
   atommnatypes = mol%get_atommnatypes()

   do i = 1, mol%natom
      call groupbytype(nadjs(i), adjlists(:, i), atommnatypes, natomneimnatypes(i), atomneimnatypepartlens(:, i), adjeqvid)
      atomorder(:nadjs(i)) = sorted_order(adjeqvid, nadjs(i))
      adjlists(:nadjs(i), i) = adjlists(atomorder(:nadjs(i)), i)
   end do

   call mol%set_adjlists(nadjs, adjlists)
   call mol%set_atomneimnatypepartlens(natomneimnatypes, atomneimnatypepartlens)

end subroutine

subroutine getmnatypes(mol, nin, nout, intype, outype, typemap)
   type(cMol), intent(in) :: mol
   integer, intent(in) :: nin
   integer, dimension(:), intent(in) :: intype
   integer, intent(out) :: nout
   integer, dimension(:), intent(out) :: outype, typemap

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
                  if (same_adjacency(nin, intype, nadjs(i), adjlists(:, i), intype, nadjs(j), adjlists(:, j))) then
                     outype(j) = nout
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

   typemap(:nout) = typemap(parentype(:nout))

end subroutine

subroutine calcequivmat(mol0, mol1, mnaord0, mnaord1, nadjmna0, adjmnalen0, adjmnalist0, &
   nadjmna1, adjmnalen1, adjmnalist1, equivmat)
! Purpose: Calculate the maximum common MNA level for all atom cross assignments
   type(cMol), intent(in) :: mol0, mol1
   integer, intent(in), dimension(:) :: mnaord0, mnaord1

   integer, dimension(:, :), intent(out) :: nadjmna0, nadjmna1
   integer, dimension(:, :, :), intent(out) :: adjmnalen0, adjmnalen1
   integer, dimension(:, :, :), intent(out) :: adjmnalist0, adjmnalist1
   integer, dimension(:, :), intent(out) :: equivmat

   integer :: natom, neltype
   integer, allocatable :: eltypepartlens(:)
   integer, dimension(:), allocatable :: nadjs0, nadjs1
   integer, dimension(:, :), allocatable :: adjlists0, adjlists1

   integer :: h, i, j, offset, level, nin, nout
   integer, dimension(mol0%natom) :: intype0, intype1, outype0, outype1
   integer, dimension(maxcoord) :: indices, atomorder

   natom = mol0%get_natom()
   neltype = mol0%get_neltype()
   eltypepartlens = mol0%get_eltypepartlens()
   nadjs0 = mol0%get_nadjs(mnaord0)
   nadjs1 = mol1%get_nadjs(mnaord1)
   adjlists0 = mol0%get_adjlists(mnaord0)
   adjlists1 = mol1%get_adjlists(mnaord1)

   nin = neltype
   intype0 = mol0%get_atomeltypes(mnaord0)
   intype1 = mol1%get_atomeltypes(mnaord1)
   level = 1

   do

      offset = 0
      do h = 1, neltype
         do i = offset + 1, offset + eltypepartlens(h)
            do j = offset + 1, offset + eltypepartlens(h)
               if (intype0(i) == intype1(j)) then
                  equivmat(j, i) = level
               end if
            end do
         end do
         offset = offset + eltypepartlens(h)
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

function same_adjacency(neltype, atomtype0, nadjs0, adjlists0, atomtype1, nadjs1, adjlists1) result(sameadj)
   integer, intent(in) :: neltype, nadjs0, nadjs1
   integer, dimension(:), intent(in) :: adjlists0, adjlists1
   integer, dimension(:) :: atomtype0, atomtype1
   logical :: sameadj

   integer :: i0, i1
   integer, dimension(neltype) :: n0, n1
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
