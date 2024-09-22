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
   type(t_mol), intent(inout) :: mol

   integer :: i, j
   integer :: neltype
   integer :: eltypes(mol%natom)
   logical :: remaining(mol%natom)
   integer, dimension(mol%natom) :: elnums, labels
   integer, dimension(mol%natom) :: blockelnums, blocklabels
   integer, dimension(mol%natom) :: foreorder, backorder

!   real(rk), allocatable :: weights(:)

   ! Initialization

   neltype = 0
   remaining = .true.
   elnums = mol%get_elnums()
   labels = mol%get_labels()
!   weights = mol%get_weights()

   ! Create block list

   do i = 1, mol%natom

      if (remaining(i)) then

         neltype = neltype + 1
         eltypes(i) = neltype
         blockelnums(neltype) = elnums(i)
         blocklabels(neltype) = labels(i)
         remaining(i) = .false.

         do j = 1, mol%natom
            if (remaining(j)) then
               if (elnums(i) == elnums(j) .and. labels(i) == labels(j)) then
!                  if (weights(i) == weights(j)) then
                     remaining(j) = .false.
                     eltypes(j) = neltype
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
   eltypes = backorder(eltypes)

   ! Order parts by atomic number

   foreorder(:neltype) = sorted_order(blockelnums, neltype)
   backorder(:neltype) = inverse_mapping(foreorder(:neltype))
   eltypes = backorder(eltypes)

   call mol%set_eltypes(neltype, eltypes) 

end subroutine

subroutine set_equiv_atoms(mol)
! Group atoms by MNA type at infinite level
   type(t_mol), intent(inout) :: mol

   integer :: i, natom, nin, nout
   integer, allocatable, dimension(:) :: intypes, outtypes, typemap
   integer, allocatable, dimension(:) :: foreorder, backorder

   allocate(outtypes(size(mol%atoms)))
   allocate(typemap(size(mol%atoms)))
   allocate(foreorder(size(mol%atoms)))
   allocate(backorder(size(mol%atoms)))

   ! Initialization

   nin = mol%get_neltype()
   intypes = mol%get_eltypes()
   typemap(:nin) = [(i, i=1, nin)]

   ! Determine MNA types iteratively

   do

      call getmnatypes(mol, nin, nout, intypes, outtypes, typemap)

      if (all(outtypes == intypes)) exit

      nin = nout
      intypes = outtypes

   end do

   foreorder(:nout) = sorted_order(typemap, nout)
   backorder(:nout) = inverse_mapping(foreorder(:nout))
   outtypes = backorder(outtypes)

   call mol%set_mnatypes(nout, outtypes)

end subroutine

subroutine assort_neighbors(mol)
! Purpose: Categorize atoms by eqtypes
   type(t_mol), intent(inout) :: mol

   integer :: i, h
   integer :: nadjs(mol%natom), adjlists(maxcoord, mol%natom)
   integer :: nadjmnatypes(mol%natom), adjmnatypepartlens(maxcoord, mol%natom)
   integer :: adjeqvid(maxcoord), atomorder(maxcoord), mnatypes(mol%natom)

   nadjs = mol%get_nadjs()
   adjlists = mol%get_adjlists()
   mnatypes = mol%get_mnatypes()

   do i = 1, mol%natom
      call groupbytype(adjlists(:nadjs(i), i), mnatypes, nadjmnatypes(i), adjmnatypepartlens(:, i), adjeqvid)
      atomorder(:nadjs(i)) = sorted_order(adjeqvid, nadjs(i))
      adjlists(:nadjs(i), i) = adjlists(atomorder(:nadjs(i)), i)
   end do

   call mol%set_adjlists(nadjs, adjlists)
   call mol%set_adjmnatypepartlens(nadjmnatypes, adjmnatypepartlens)

end subroutine

subroutine getmnatypes(mol, nin, nout, intypes, outtypes, typemap)
   type(t_mol), intent(in) :: mol
   integer, intent(in) :: nin
   integer, dimension(:), intent(in) :: intypes
   integer, intent(out) :: nout
   integer, dimension(:), intent(out) :: outtypes, typemap

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
         outtypes(i) = nout
         parentype(nout) = intypes(i)
         do j = i + 1, mol%natom
!               print '(a, x, i0, x, i0)', trim(elsym(intypes0(i))), i, j
            if (untyped(j)) then
               if (intypes(j) == intypes(i)) then
                  if (same_adjacency(nin, intypes, adjlists(:nadjs(i), i), intypes, adjlists(:nadjs(j), j))) then
                     outtypes(j) = nout
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

   typemap(:nout) = typemap(parentype(:nout))

end subroutine

subroutine calcequivmat(mol0, mol1, nadjmna0, adjmnalen0, adjmnalist0, &
   nadjmna1, adjmnalen1, adjmnalist1, equivmat)
! Purpose: Calculate the maximum common MNA level for all atom cross assignments
   type(t_mol), intent(in) :: mol0, mol1

   integer, dimension(:, :), intent(out) :: nadjmna0, nadjmna1
   integer, dimension(:, :, :), intent(out) :: adjmnalen0, adjmnalen1
   integer, dimension(:, :, :), intent(out) :: adjmnalist0, adjmnalist1
   integer, dimension(:, :), intent(out) :: equivmat

   integer :: natom, neltype
   integer, allocatable :: eltypepartlens(:)
!   integer, dimension(:), allocatable :: nadjs0, nadjs1
!   integer, dimension(:, :), allocatable :: adjlists0, adjlists1

   integer :: h, i, j, offset, level, nin, nout
   integer, dimension(mol0%natom) :: intypes0, intypes1, outtypes0, outtypes1
   integer, dimension(maxcoord) :: indices, atomorder
   type(t_atomlist), allocatable, dimension(:) :: adjlists0, adjlists1

   natom = mol0%natom
   neltype = mol0%get_neltype()
   eltypepartlens = mol0%get_eltypepartlens()
   adjlists0 = mol0%get_sorted_newadjlists()
   adjlists1 = mol1%get_sorted_newadjlists()

   nin = neltype
   intypes0 = mol0%get_sorted_eltypes()
   intypes1 = mol1%get_sorted_eltypes()
   level = 1

   do

      offset = 0
      do h = 1, neltype
         do i = offset + 1, offset + eltypepartlens(h)
            do j = offset + 1, offset + eltypepartlens(h)
               if (intypes0(i) == intypes1(j)) then
                  equivmat(j, i) = level
               end if
            end do
         end do
         offset = offset + eltypepartlens(h)
      end do

      do i = 1, natom
         call groupbytype(adjlists0(i)%atomidcs, intypes0, nadjmna0(i, level), adjmnalen0(:, i, level), indices)
         atomorder(:size(adjlists0(i)%atomidcs)) = sorted_order(indices, size(adjlists0(i)%atomidcs))
         adjmnalist0(:size(adjlists0(i)%atomidcs), i, level) = adjlists0(i)%atomidcs(atomorder(:size(adjlists0(i)%atomidcs)))
      end do

      do i = 1, natom
         call groupbytype(adjlists1(i)%atomidcs, intypes1, nadjmna1(i, level), adjmnalen1(:, i, level), indices)
         atomorder(:size(adjlists1(i)%atomidcs)) = sorted_order(indices, size(adjlists1(i)%atomidcs))
         adjmnalist1(:size(adjlists1(i)%atomidcs), i, level) = adjlists1(i)%atomidcs(atomorder(:size(adjlists1(i)%atomidcs)))
      end do

      call getmnacrosstypes(adjlists0, adjlists1, nin, intypes0, intypes1, nout, outtypes0, outtypes1)

      if (all(outtypes0 == intypes0) .and. all((outtypes1 == intypes1))) exit

      nin = nout
      intypes0 = outtypes0
      intypes1 = outtypes1
      level = level + 1

   end do

end subroutine

subroutine getmnacrosstypes(adjlists0, adjlists1, nin, intypes0, intypes1, &
      nout, outtypes0, outtypes1)
   type(t_atomlist), dimension(:), intent(in) :: adjlists0, adjlists1
   integer, intent(in) :: nin
   integer, dimension(:), intent(in) :: intypes0, intypes1
!   integer, dimension(:), intent(in) :: nadjs0, nadjs1
!   integer, dimension(:, :), intent(in) :: adjlists0, adjlists1
   integer, intent(out) :: nout
   integer, dimension(:), intent(out) :: outtypes0, outtypes1

   integer :: i, j
   integer :: archatom(size(adjlists0))
   logical :: untyped(size(adjlists0))

   nout = 0
   untyped(:) = .true.

   do i = 1, size(adjlists0)
      if (untyped(i)) then
         nout = nout + 1
         outtypes0(i) = nout
         archatom(nout) = i
         do j = i + 1, size(adjlists0)
            if (untyped(j)) then
               if (intypes0(j) == intypes0(i)) then
                  if (same_adjacency(nin, intypes0, adjlists0(i)%atomidcs, intypes0, adjlists0(j)%atomidcs)) then
                     outtypes0(j) = nout
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

   untyped(:) = .true.

   do i = 1, nout
      do j = 1, size(adjlists0)
         if (untyped(j)) then
            if (intypes1(j) == intypes0(archatom(i))) then
               if (same_adjacency(nin, intypes0, adjlists0(archatom(i))%atomidcs, intypes1, adjlists1(j)%atomidcs)) then
                  outtypes1(j) = i
                  untyped(j) = .false.
               end if
            end if
         end if
      end do
   end do

   do i = 1, size(adjlists0)
      if (untyped(i)) then
         nout = nout + 1
         outtypes1(i) = nout
         do j = i + 1, size(adjlists0)
            if (untyped(j)) then
               if (intypes1(j) == intypes1(i)) then
                  if (same_adjacency(nin, intypes1, adjlists1(i)%atomidcs, intypes1, adjlists1(j)%atomidcs)) then
                     outtypes1(j) = nout
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

end subroutine

subroutine groupbytype(items, itemtypes, ngroup, grouplens, groupidcs)
! Purpose: Categorize atoms by types
    integer, intent(out) :: ngroup
    integer, dimension(:), intent(in) :: itemtypes, items
    integer, dimension(:), intent(out) :: groupidcs, grouplens

    integer :: i, j
    integer, dimension(maxcoord) :: grouptype
    logical, dimension(maxcoord) :: remaining
    integer, dimension(maxcoord) :: foreorder, backorder

! Initialization

    ngroup = 0
    remaining = .true.

! Create group lists

    do i = 1, size(items)
        if (remaining(i)) then
            ngroup = ngroup + 1
            grouplens(ngroup) = 1
            grouptype(ngroup) = itemtypes(items(i))
            groupidcs(i) = ngroup
            do j = i + 1, size(items)
                if (remaining(i)) then
                    if (itemtypes(items(j)) == itemtypes(items(i))) then
                        groupidcs(j) = ngroup
                        grouplens(ngroup) = grouplens(ngroup) + 1
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
    grouplens(:ngroup) = grouplens(foreorder(:ngroup))
    groupidcs(:size(items)) = backorder(groupidcs(:size(items)))

!    print *, groupidcs(:ngroup)
!    print *, itemtypes(foreorder)

end subroutine

function same_adjacency(neltype, atomtype0, adjlist0, atomtype1, adjlist1) result(sameadj)
   integer, intent(in) :: neltype
   integer, dimension(:), intent(in) :: adjlist0, adjlist1
   integer, dimension(:) :: atomtype0, atomtype1
   logical :: sameadj

   integer :: i0, i1
   integer, dimension(neltype) :: n0, n1
!   real :: adjlists0(3, maxcoord), adjlists1(3, maxcoord)
!   integer :: typelist0(maxcoord, nin), typelist1(maxcoord, nin)

   sameadj = .true.

   if (size(adjlist0) /= size(adjlist1)) then
      sameadj = .false.
      return
   end if

!   If coordination number is the same check if coordinated atoms are the same

   n0(:) = 0
   n1(:) = 0

   do i0 = 1, size(adjlist0)
      n0(atomtype0(adjlist0(i0))) = n0(atomtype0(adjlist0(i0))) + 1
!       typelist0(n0(atomtype0(adjlist0(i0))), atomtype0(adjlist0(i0))) = i0
   end do

   do i1 = 1, size(adjlist1)
      n1(atomtype1(adjlist1(i1))) = n1(atomtype1(adjlist1(i1))) + 1
!       typelist1(n1(atomtype1(adjlist1(i1))), atomtype1(adjlist1(i1))) = i1
   end do

   if (any(n0 /= n1)) then
      sameadj = .false.
      return
   end if

!   print *, typelist0(:size(adjlist0)), '/', typelist1(:size(adjlist1))

end function

end module
