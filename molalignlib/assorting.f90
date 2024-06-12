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
use pointers
use discrete
use sorting
use chemdata
use chemutils

implicit none

contains

subroutine assort_atoms(mol)
! Purpose: Group atoms by atomic numbers and types
   type(Molecule), intent(inout) :: mol

   integer :: i, j
   integer :: ntype, ztype
   integer :: typeidcs(mol%natom)
   logical :: remaining(mol%natom)
   integer, dimension(mol%natom) :: typeaggs, blkznum, blkynum
   integer, dimension(mol%natom) :: order, idorder

   integer, allocatable :: znums(:)
   real(wp), allocatable :: weights(:)
   character(wl), allocatable :: labels(:)

   ! Initialization

   ntype = 0
   remaining = .true.
   labels = mol%get_labels()
!   weights = mol%get_weights()
   allocate(znums(mol%natom))
   allocate(weights(mol%natom))

   ! Create block list

   do i = 1, mol%natom

      if (remaining(i)) then

         call readlabel(labels(i), znums(i), ztype)
         weights = weight_func(znums(i))

         ntype = ntype + 1
         typeidcs(i) = ntype
         typeaggs(ntype) = 1
         blkznum(ntype) = znums(i)
         blkynum(ntype) = ztype
         remaining(i) = .false.

         do j = 1, mol%natom
            if (remaining(j)) then
               if (labels(i) == labels(j)) then
!                  if (weights(i) == weights(j)) then
                     typeidcs(j) = ntype
                     znums(j) = znums(i)
                     weights(j) = weights(i)
                     remaining(j) = .false.
                     typeaggs(ntype) = typeaggs(ntype) + 1
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

   ! Order blocks by type number

   order(:ntype) = sorted_order(blkynum, ntype)
   idorder(:ntype) = inverse_permut(order(:ntype))
   typeaggs(:ntype) = typeaggs(order(:ntype))
   typeidcs = idorder(typeidcs)

   ! Order blocks by atomic number

   order(:ntype) = sorted_order(blkznum, ntype)
   idorder(:ntype) = inverse_permut(order(:ntype))
   typeaggs(:ntype) = typeaggs(order(:ntype))
   typeidcs = idorder(typeidcs)

   call mol%set_znums(znums)
   call mol%set_weights(weights)
   call mol%set_typeidcs(typeidcs) 
   call mol%set_typeaggs(ntype, typeaggs)

end subroutine

subroutine set_equiv_atoms(mol)
! Group atoms by MNA type at infinite level
   type(Molecule), intent(inout) :: mol

   integer i, nin, nequiv
   integer, dimension(mol%natom) :: equividcs, equivaggs
   integer, dimension(mol%natom) :: intype, uptype, basetype
   integer, dimension(mol%natom) :: order, idorder

   ! Determine MNA types iteratively

   nin = mol%get_ntype()
   intype = mol%get_typeidcs()
   basetype = [(i, i=1, mol%natom)]

   do

      call getmnatypes(mol, nin, intype, nequiv, equividcs, equivaggs, uptype)
      basetype(:nequiv) = basetype(uptype(:nequiv))

      if (all(equividcs == intype)) exit

      nin = nequiv
      intype = equividcs

   end do

   order(:nequiv) = sorted_order(basetype, nequiv)
   idorder(:nequiv) = inverse_permut(order(:nequiv))
   equivaggs(:nequiv) = equivaggs(order(:nequiv))
   equividcs = idorder(equividcs)

!    do i = 1, mol%natom
!        print *, i, elsym(znum(i)), equividcs(i)
!    end do

   call mol%set_equividcs(equividcs)
   call mol%set_equivaggs(nequiv, equivaggs)

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
    integer, dimension(maxcoord) :: order, idorder

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

    order(:ngroup) = sorted_order(grouptype, ngroup)
    idorder(:ngroup) = inverse_permut(order(:ngroup))
    grouptype(:ngroup) = grouptype(order(:ngroup))
    groupsize(:ngroup) = groupsize(order(:ngroup))
    groupid(:nelem) = idorder(groupid(:nelem))

!    print *, groupid(:ngroup)
!    print *, typess(order)

end subroutine

subroutine assort_neighbors(mol)
! Purpose: Categorize atoms by eqtypes
   type(Molecule), intent(inout) :: mol

   integer :: i, h
   integer :: coonums(mol%natom), neighbors(maxcoord, mol%natom)
   integer :: adjeqvid(maxcoord), atomorder(maxcoord), equividcs(mol%natom)

   coonums = mol%get_coonums()
   neighbors = mol%get_neighbors()
   equividcs = mol%get_equividcs()

   do i = 1, mol%natom
      if (.not. allocated(mol%atoms(i)%neieqvlens)) then
         allocate(mol%atoms(i)%neieqvlens(maxcoord))
      end if
      call groupbytype(coonums(i), neighbors(:, i), equividcs, &
            mol%atoms(i)%nneieqv, mol%atoms(i)%neieqvlens, adjeqvid)
      atomorder(:coonums(i)) = sorted_order(adjeqvid, coonums(i))
      neighbors(:coonums(i), i) = neighbors(atomorder(:coonums(i)), i)
   end do

   call mol%set_neighbors(coonums, neighbors)

end subroutine

subroutine getmnatypes(mol, nin, intype, nout, outype, outsize, uptype)
   type(Molecule), intent(inout) :: mol
   integer, intent(in) :: nin
   integer, dimension(:), intent(in) :: intype
   integer, intent(out) :: nout
   integer, dimension(:), intent(out) :: outype, outsize, uptype

   integer :: i, j
   integer :: coonums(mol%natom), neighbors(maxcoord, mol%natom)
   logical :: untyped(mol%natom)

   coonums = mol%get_coonums()
   neighbors = mol%get_neighbors()

   nout = 0
   untyped(:) = .true.

   do i = 1, mol%natom
      if (untyped(i)) then
         nout = nout + 1
         outype(i) = nout
         outsize(nout) = 1
         uptype(nout) = intype(i)
         do j = i + 1, mol%natom
!               print '(a, x, i0, x, i0)', trim(elsym(intype0(i))), i, j
            if (untyped(j)) then
               if (intype(j) == intype(i)) then
                  if (same_adjacency(nin, intype, coonums(i), neighbors(:, i), intype, &
                        coonums(j), neighbors(:, j))) then
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

subroutine calcequivmat(mol0, mol1, nadjmna0, adjmnalen0, adjmnalist0, &
   nadjmna1, adjmnalen1, adjmnalist1, equivmat)
! Purpose: Calculate the maximum common MNA level for all atom cross assignments
   type(Molecule), intent(in) :: mol0, mol1

   integer, dimension(:, :), intent(out) :: nadjmna0, nadjmna1
   integer, dimension(:, :, :), intent(out) :: adjmnalen0, adjmnalen1
   integer, dimension(:, :, :), intent(out) :: adjmnalist0, adjmnalist1
   integer, dimension(:, :), intent(out) :: equivmat

   integer :: natom, ntype
   integer, allocatable :: typeaggs(:)
   integer, dimension(:), allocatable :: coonums0, coonums1
   integer, dimension(:, :), allocatable :: neighbors0, neighbors1

   integer :: h, i, j, offset, level, nin, nout
   integer, dimension(mol0%natom) :: intype0, intype1, outype0, outype1
   integer, dimension(maxcoord) :: indices, atomorder

   natom = mol0%get_natom()
   ntype = mol0%get_ntype()
   typeaggs = mol0%get_typeaggs()
   coonums0 = mol0%get_coonums()
   coonums1 = mol1%get_coonums()
   neighbors0 = mol0%get_neighbors()
   neighbors1 = mol1%get_neighbors()

   nin = ntype
   intype0 = mol0%get_typeidcs()
   intype1 = mol1%get_typeidcs()
   level = 1

   do

      offset = 0
      do h = 1, ntype
         do i = offset + 1, offset + typeaggs(h)
            do j = offset + 1, offset + typeaggs(h)
               if (intype0(i) == intype1(j)) then
                  equivmat(j, i) = level
               end if
            end do
         end do
         offset = offset + typeaggs(h)
      end do

      do i = 1, natom
         call groupbytype(coonums0(i), neighbors0(:, i), intype0, nadjmna0(i, level), adjmnalen0(:, i, level), indices)
         atomorder(:coonums0(i)) = sorted_order(indices, coonums0(i))
         adjmnalist0(:coonums0(i), i, level) = neighbors0(atomorder(:coonums0(i)), i)
      end do

      do i = 1, natom
         call groupbytype(coonums1(i), neighbors1(:, i), intype1, nadjmna1(i, level), adjmnalen1(:, i, level), indices)
         atomorder(:coonums1(i)) = sorted_order(indices, coonums1(i))
         adjmnalist1(:coonums1(i), i, level) = neighbors1(atomorder(:coonums1(i)), i)
      end do

      call getmnacrosstypes(mol0%natom, nin, intype0, coonums0, neighbors0, intype1, coonums1, neighbors1, &
         nout, outype0, outype1)

      if (all(outype0 == intype0) .and. all((outype1 == intype1))) exit

      nin = nout
      intype0 = outype0
      intype1 = outype1
      level = level + 1

   end do

!    print *, natom, nout

end subroutine

subroutine getmnacrosstypes(natom, nin, intype0, coonums0, neighbors0, intype1, coonums1, neighbors1, &
      nout, outype0, outype1)
   integer, intent(in) :: natom, nin
   integer, dimension(:), intent(in) :: intype0, intype1, coonums0, coonums1
   integer, dimension(:, :), intent(in) :: neighbors0, neighbors1
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
                  if (same_adjacency(nin, intype0, coonums0(i), neighbors0(:, i), intype0, coonums0(j), neighbors0(:, j))) then
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
               if (same_adjacency(nin, intype0, coonums0(archetype(i)), neighbors0(:, archetype(i)), intype1, coonums1(j), &
                            neighbors1(:, j))) then
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
                  if (same_adjacency(nin, intype1, coonums1(i), neighbors1(:, i), intype1, coonums1(j), neighbors1(:, j))) then
                     outype1(j) = nout
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

end subroutine

function same_adjacency(ntype, atomtype0, coonums0, neighbors0, atomtype1, coonums1, neighbors1) result(sameadj)
   integer, intent(in) :: ntype, coonums0, coonums1
   integer, dimension(:), intent(in) :: neighbors0, neighbors1
   integer, dimension(:) :: atomtype0, atomtype1
   logical :: sameadj

   integer :: i0, i1
   integer, dimension(ntype) :: n0, n1
!   real :: atoms0(3, maxcoord), atoms1(3, maxcoord)
!   integer :: typelist0(maxcoord, nin), typelist1(maxcoord, nin)

   sameadj = .true.

   if (coonums0 /= coonums1) then
      sameadj = .false.
      return
   end if

!   If coordination number is the same check if coordinated atoms are the same

   n0(:) = 0
   n1(:) = 0

   do i0 = 1, coonums0
      n0(atomtype0(neighbors0(i0))) = n0(atomtype0(neighbors0(i0))) + 1
!       typelist0(n0(atomtype0(neighbors0(i0))), atomtype0(neighbors0(i0))) = i0
   end do

   do i1 = 1, coonums1
      n1(atomtype1(neighbors1(i1))) = n1(atomtype1(neighbors1(i1))) + 1
!       typelist1(n1(atomtype1(neighbors1(i1))), atomtype1(neighbors1(i1))) = i1
   end do

   if (any(n0 /= n1)) then
      sameadj = .false.
      return
   end if

!   print *, typelist0(:coonums0), '/', typelist1(:coonums1)

end function

end module
