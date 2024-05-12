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

   integer :: i, j, nblock
   integer :: znumi, ztypei, znumj, ztypej
   integer :: blkid(mol%natom)
   logical :: remaining(mol%natom)
   integer, dimension(mol%natom) :: blklens, blkznum, blkynum
   integer, dimension(mol%natom) :: order, idorder

   ! Initialization

   nblock = 0
   remaining = .true.

   ! Create block list

   do i = 1, mol%natom

      if (remaining(i)) then

         call readlabel(mol%atoms(i)%label, znumi, ztypei)

         nblock = nblock + 1
         blkid(i) = nblock
         blklens(nblock) = 1
         blkznum(nblock) = znumi
         blkynum(nblock) = ztypei
         mol%atoms(i)%znum = znumi
         mol%atoms(i)%weight = weight_func(znumi)
         remaining(i) = .false.

         do j = 1, mol%natom
            if (remaining(j)) then
               call readlabel(mol%atoms(j)%label, znumj, ztypej)
               if (znumi == znumj .and. ztypei == ztypej) then
                  if (weight_func(znumi) == weight_func(znumj)) then
                     blkid(j) = nblock
                     mol%atoms(j)%znum = znumj
                     mol%atoms(j)%weight = weight_func(znumj)
                     remaining(j) = .false.
                     blklens(nblock) = blklens(nblock) + 1
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

   ! Order blocks by type number

   order(:nblock) = sorted_order(blkynum, nblock)
   idorder(:nblock) = inverse_permut(order(:nblock))
   blklens(:nblock) = blklens(order(:nblock))
   blkid = idorder(blkid)

   ! Order blocks by atomic number

   order(:nblock) = sorted_order(blkznum, nblock)
   idorder(:nblock) = inverse_permut(order(:nblock))
   blklens(:nblock) = blklens(order(:nblock))
   blkid = idorder(blkid)

   allocate (mol%blklens(nblock))

   mol%nblock = nblock
   mol%blklens = blklens(:nblock)
   mol%atoms(:)%blkid = blkid(:)

end subroutine

subroutine set_equiv_atoms(mol)
! Group atoms by MNA type at infinite level
   type(Molecule), intent(inout) :: mol

   integer i, nin, nequiv
   integer, dimension(mol%natom) :: eqvid, eqvlens
   integer, dimension(mol%natom) :: intype, uptype, basetype
   integer, dimension(mol%natom) :: order, idorder

   ! Determine MNA types iteratively

   nin = mol%nblock
   intype = mol%get_blkids()
   basetype = [(i, i=1, mol%natom)]

   do

      call getmnatypes(mol, nin, intype, nequiv, eqvid, eqvlens, uptype)
      basetype(:nequiv) = basetype(uptype(:nequiv))

      if (all(eqvid == intype)) exit

      nin = nequiv
      intype = eqvid

   end do

   order(:nequiv) = sorted_order(basetype, nequiv)
   idorder(:nequiv) = inverse_permut(order(:nequiv))
   eqvlens(:nequiv) = eqvlens(order(:nequiv))
   eqvid = idorder(eqvid)

!    do i = 1, mol%natom
!        print *, i, elsym(znum(i)), eqvid(i)
!    end do

   allocate (mol%eqvlens(nequiv))

   mol%nequiv = nequiv
   mol%eqvlens = eqvlens(:nequiv)
   mol%atoms(:)%eqvid = eqvid(:)

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
   integer, dimension(maxcoord) :: adjeqvid, atomorder

   do i = 1, mol%natom
      allocate(mol%atoms(i)%adjeqvlens(maxcoord))
      call groupbytype(mol%atoms(i)%nadj, mol%atoms(i)%adjlist, mol%atoms(:)%eqvid, &
            mol%atoms(i)%nadjeqv, mol%atoms(i)%adjeqvlens, adjeqvid)
      atomorder(:mol%atoms(i)%nadj) = sorted_order(adjeqvid, mol%atoms(i)%nadj)
      mol%atoms(i)%adjlist(:mol%atoms(i)%nadj) = mol%atoms(i)%adjlist(atomorder(:mol%atoms(i)%nadj))
   end do

end subroutine

subroutine getmnatypes(mol, nin, intype, nout, outype, outsize, uptype)
   type(Molecule), intent(inout) :: mol
   integer, intent(in) :: nin
   integer, dimension(:), intent(in) :: intype
   integer, intent(out) :: nout
   integer, dimension(:), intent(out) :: outype, outsize, uptype

   integer :: i, j
   logical :: untyped(mol%natom)

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
                  if (same_adjacency(nin, intype, mol%atoms(i)%nadj, mol%atoms(i)%adjlist, intype, &
                        mol%atoms(j)%nadj, mol%atoms(j)%adjlist)) then
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

   integer :: natom, nblock
   integer, allocatable :: blklens(:)
   integer, dimension(:), allocatable :: nadj0, nadj1
   integer, dimension(:, :), allocatable :: adjlist0, adjlist1

   integer :: h, i, j, offset, level, nin, nout
   integer, dimension(mol0%natom) :: intype0, intype1, outype0, outype1
   integer, dimension(maxcoord) :: indices, atomorder

   natom = mol0%get_natom()
   nblock = mol0%get_nblock()
   blklens = mol0%get_blklens()
   nadj0 = mol0%get_nadj()
   nadj1 = mol1%get_nadj()
   adjlist0 = mol0%get_adjlist()
   adjlist1 = mol1%get_adjlist()

   nin = nblock
   intype0 = mol0%get_blkids()
   intype1 = mol1%get_blkids()
   level = 1

   do

      offset = 0
      do h = 1, nblock
         do i = offset + 1, offset + blklens(h)
            do j = offset + 1, offset + blklens(h)
               if (intype0(i) == intype1(j)) then
                  equivmat(j, i) = level
               end if
            end do
         end do
         offset = offset + blklens(h)
      end do

      do i = 1, natom
         call groupbytype(nadj0(i), adjlist0(:, i), intype0, nadjmna0(i, level), adjmnalen0(:, i, level), indices)
         atomorder(:nadj0(i)) = sorted_order(indices, nadj0(i))
         adjmnalist0(:nadj0(i), i, level) = adjlist0(atomorder(:nadj0(i)), i)
      end do

      do i = 1, natom
         call groupbytype(nadj1(i), adjlist1(:, i), intype1, nadjmna1(i, level), adjmnalen1(:, i, level), indices)
         atomorder(:nadj1(i)) = sorted_order(indices, nadj1(i))
         adjmnalist1(:nadj1(i), i, level) = adjlist1(atomorder(:nadj1(i)), i)
      end do

      call getmnacrosstypes(mol0%natom, nin, intype0, nadj0, adjlist0, intype1, nadj1, adjlist1, &
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
                  if (same_adjacency(nin, intype0, nadj0(i), adjlist0(:, i), intype0, nadj0(j), adjlist0(:, j))) then
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
               if (same_adjacency(nin, intype0, nadj0(archetype(i)), adjlist0(:, archetype(i)), intype1, nadj1(j), &
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
                  if (same_adjacency(nin, intype1, nadj1(i), adjlist1(:, i), intype1, nadj1(j), adjlist1(:, j))) then
                     outype1(j) = nout
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

end subroutine

function same_adjacency(ntype, atomtype0, nadj0, adjlist0, atomtype1, nadj1, adjlist1) result(sameadj)
   integer, intent(in) :: ntype, nadj0, nadj1
   integer, dimension(:), intent(in) :: adjlist0, adjlist1
   integer, dimension(:) :: atomtype0, atomtype1
   logical :: sameadj

   integer :: i0, i1
   integer, dimension(ntype) :: n0, n1
!   real :: atoms0(3, maxcoord), atoms1(3, maxcoord)
!   integer :: typelist0(maxcoord, nin), typelist1(maxcoord, nin)

   sameadj = .true.

   if (nadj0 /= nadj1) then
      sameadj = .false.
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
      sameadj = .false.
      return
   end if

!   print *, typelist0(:nadj0), '/', typelist1(:nadj1)

end function

end module
