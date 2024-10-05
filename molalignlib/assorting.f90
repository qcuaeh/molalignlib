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
use partition
use permutation
use molecule
use flags
use bounds
use sorting
use chemdata

implicit none

contains

! Group atoms by atomic number and label
subroutine resolve_eltypes(mol)
   type(molecule_type), intent(inout) :: mol
   ! Local variables
   integer :: i, j
   integer :: neltype
   integer :: atomeltypes(mol%natom)
   integer :: eltypepops(mol%natom)
   logical :: untyped(mol%natom)
   integer, dimension(mol%natom) :: elnums, labels
   integer, dimension(mol%natom) :: groupelnums, grouplabels
   integer, dimension(mol%natom) :: foreorder, backorder
!   real(rk), allocatable :: weights(:)

   neltype = 0
   untyped(:) = .true.
   elnums = mol%get_elnums()
   labels = mol%get_labels()
!   weights = mol%get_weights()

   do i = 1, mol%natom
      if (untyped(i)) then
         neltype = neltype + 1
         eltypepops(neltype) = 1
         groupelnums(neltype) = elnums(i)
         grouplabels(neltype) = labels(i)
         atomeltypes(i) = neltype
         untyped(i) = .false.
         do j = 1, mol%natom
            if (untyped(j)) then
               if (elnums(i) == elnums(j) .and. labels(i) == labels(j)) then
!                  if (weights(i) == weights(j)) then
                     untyped(j) = .false.
                     atomeltypes(j) = neltype
                     eltypepops(neltype) = eltypepops(neltype) + 1
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
   foreorder(:neltype) = sorted_order(grouplabels, neltype)
   backorder(:neltype) = inverse_permutation(foreorder(:neltype))
   eltypepops(:neltype) = eltypepops(foreorder(:neltype))
   atomeltypes = backorder(atomeltypes)

   ! Order parts by atomic number
   foreorder(:neltype) = sorted_order(groupelnums, neltype)
   backorder(:neltype) = inverse_permutation(foreorder(:neltype))
   eltypepops(:neltype) = eltypepops(foreorder(:neltype))
   atomeltypes = backorder(atomeltypes)

   call mol%set_eltypes(eltypepops(:neltype), atomeltypes)

end subroutine

! Assort atoms by MNA type at infinite level
subroutine set_equiv_atoms(mol)
   type(molecule_type), intent(inout) :: mol
   ! Local variables
   integer :: i, natom, nintype, ntype
   integer, allocatable, dimension(:) :: intypes, types, typemap
   integer, allocatable, dimension(:) :: foreorder, backorder
   type(atomlist_type), allocatable :: adjlists(:)

   allocate (types(size(mol%atoms)))
   allocate (typemap(size(mol%atoms)))
   allocate (foreorder(size(mol%atoms)))
   allocate (backorder(size(mol%atoms)))

   nintype = mol%get_neltype()
   intypes = mol%get_atomeltypes()
   typemap(:nintype) = [(i, i=1, nintype)]
   adjlists = mol%get_newadjlists()

   ! Determine MNA types iteratively

   do
      call resolve_mnatypes(adjlists, nintype, ntype, intypes, types, typemap)
      ! Exit the loop if types are unchanged
      if (all(types == intypes)) exit
      nintype = ntype
      intypes = types
   end do

   foreorder(:ntype) = sorted_order(typemap, ntype)
   backorder(:ntype) = inverse_permutation(foreorder(:ntype))
   types = backorder(types)

   call mol%set_mnatypes(ntype, types)

end subroutine

! Assort neighbors by MNA type
subroutine assort_neighbors(mol)
   type(molecule_type), intent(inout) :: mol
   ! Local variables
   integer :: i, h
   integer :: nadjs(mol%natom), adjlists(maxcoord, mol%natom)
   integer :: nadjmnatypes(mol%natom), adjmnatypepartlens(maxcoord, mol%natom)
   integer :: adjeqvid(maxcoord), atomorder(maxcoord), mnatypes(mol%natom)

   nadjs = mol%get_nadjs()
   adjlists = mol%get_oldadjlists()
   mnatypes = mol%get_atommnatypes()

   do i = 1, mol%natom
      call group_by(adjlists(:nadjs(i), i), mnatypes, nadjmnatypes(i), adjmnatypepartlens(:, i), adjeqvid)
      atomorder(:nadjs(i)) = sorted_order(adjeqvid, nadjs(i))
      adjlists(:nadjs(i), i) = adjlists(atomorder(:nadjs(i)), i)
   end do

   call mol%set_adjlists(nadjs, adjlists)
   call mol%set_adjmnatypepartlens(nadjmnatypes, adjmnatypepartlens)

end subroutine

! Find next level MNA types
subroutine resolve_mnatypes(adjlists, nintype, ntype, intypes, types, typemap)
   type(atomlist_type), intent(in) :: adjlists(:)
   integer, intent(in) :: nintype
   integer, dimension(:), intent(in) :: intypes
   integer, intent(out) :: ntype
   integer, dimension(:), intent(out) :: types, typemap
   ! Local variables
   integer :: i, j
   logical :: untyped(size(adjlists))
   integer :: parentype(size(adjlists))

   ntype = 0
   untyped(:) = .true.

   do i = 1, size(adjlists)
      if (untyped(i)) then
         ntype = ntype + 1
         types(i) = ntype
         parentype(ntype) = intypes(i)
         do j = i + 1, size(adjlists)
            if (untyped(j)) then
               if (intypes(j) == intypes(i)) then
                  if (same_adjacency(nintype, intypes, adjlists(i)%atomidcs, intypes, adjlists(j)%atomidcs)) then
                     types(j) = ntype
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

   typemap(:ntype) = typemap(parentype(:ntype))

end subroutine

! Find next level cross MNA types between mol0 and mol1
subroutine resolve_crossmnatypes(adjlists0, adjlists1, nintype, intypes0, intypes1, &
      ntype, types0, types1)
   type(atomlist_type), dimension(:), intent(in) :: adjlists0, adjlists1
   integer, intent(in) :: nintype
   integer, dimension(:), intent(in) :: intypes0, intypes1
   integer, intent(out) :: ntype
   integer, dimension(:), intent(out) :: types0, types1
   ! Local variables
   integer :: h, i, j
   logical :: untyped(size(adjlists0))
   integer :: archeatom(size(adjlists0))

   ntype = 0
   untyped(:) = .true.

   do i = 1, size(adjlists0)
      if (untyped(i)) then
         ntype = ntype + 1
         types0(i) = ntype
         archeatom(ntype) = i
         do j = i + 1, size(adjlists0)
            if (untyped(j)) then
               if (intypes0(j) == intypes0(i)) then
                  if (same_adjacency(nintype, intypes0, adjlists0(i)%atomidcs, intypes0, adjlists0(j)%atomidcs)) then
                     types0(j) = ntype
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

   untyped(:) = .true.

   do h = 1, ntype
      do i = 1, size(adjlists0)
         if (untyped(i)) then
            if (intypes1(i) == intypes0(archeatom(h))) then
               if (same_adjacency(nintype, intypes0, adjlists0(archeatom(h))%atomidcs, intypes1, adjlists1(i)%atomidcs)) then
                  types1(i) = h
                  untyped(i) = .false.
               end if
            end if
         end if
      end do
   end do

   do i = 1, size(adjlists0)
      if (untyped(i)) then
         ntype = ntype + 1
         types1(i) = ntype
         do j = i + 1, size(adjlists0)
            if (untyped(j)) then
               if (intypes1(j) == intypes1(i)) then
                  if (same_adjacency(nintype, intypes1, adjlists1(i)%atomidcs, intypes1, adjlists1(j)%atomidcs)) then
                     types1(j) = ntype
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

end subroutine

! Group atoms by types
subroutine group_by(items, itemtypes, ngroup, grouplens, groupidcs)
   integer, intent(out) :: ngroup
   integer, dimension(:), intent(in) :: itemtypes, items
   integer, dimension(:), intent(out) :: groupidcs, grouplens
   ! Local variables
   integer :: i, j
   integer, dimension(maxcoord) :: grouptype
   logical, dimension(maxcoord) :: untyped
   integer, dimension(maxcoord) :: foreorder, backorder

   ngroup = 0
   untyped = .true.

   do i = 1, size(items)
       if (untyped(i)) then
           ngroup = ngroup + 1
           grouplens(ngroup) = 1
           grouptype(ngroup) = itemtypes(items(i))
           groupidcs(i) = ngroup
           do j = i + 1, size(items)
               if (untyped(i)) then
                   if (itemtypes(items(j)) == itemtypes(items(i))) then
                       groupidcs(j) = ngroup
                       grouplens(ngroup) = grouplens(ngroup) + 1
                       untyped(j) = .false.
                   end if
               end if
           end do
       end if
   end do

! Order groups by type 

    foreorder(:ngroup) = sorted_order(grouptype, ngroup)
    backorder(:ngroup) = inverse_permutation(foreorder(:ngroup))
    grouptype(:ngroup) = grouptype(foreorder(:ngroup))
    grouplens(:ngroup) = grouplens(foreorder(:ngroup))
    groupidcs(:size(items)) = backorder(groupidcs(:size(items)))

end subroutine

! Test if adjacent atoms to a pair of atoms in mol0 and mol1 are the same
function same_adjacency(neltype, atomtype0, adjlist0, atomtype1, adjlist1) result(sameadj)
   integer, intent(in) :: neltype
   integer, dimension(:), intent(in) :: adjlist0, adjlist1
   integer, dimension(:) :: atomtype0, atomtype1
   ! Result variable
   logical :: sameadj
   ! Local variables
   integer :: i0, i1
   integer, dimension(neltype) :: n0, n1

   sameadj = .true.

   if (size(adjlist0) /= size(adjlist1)) then
      sameadj = .false.
      return
   end if

!  Find if adjacent atoms are the same

   n0(:) = 0
   n1(:) = 0

   do i0 = 1, size(adjlist0)
      n0(atomtype0(adjlist0(i0))) = n0(atomtype0(adjlist0(i0))) + 1
   end do

   do i1 = 1, size(adjlist1)
      n1(atomtype1(adjlist1(i1))) = n1(atomtype1(adjlist1(i1))) + 1
   end do

   if (any(n0 /= n1)) then
      sameadj = .false.
      return
   end if

end function

end module
