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

module remapping
use kinds
use types
use flags
use bounds
use random
use strutils
use adjacency
use rotation
use translation
use tracking
use backtracking
use alignment
use assorting
use printing
use biasing
use assignment

implicit none

private
public optimize_mapping
public find_reactive_sites

contains

subroutine optimize_mapping(mol0, mol1, maplist, countlist, nrec)
!subroutine optimize_mapping( &
!   mol0, &
!   mol1, &
!   ntype, &
!   typeaggs, &
!   nequiv0, &
!   equivaggs0, &
!   nequiv1, &
!   equivaggs1, &
!   maplist, &
!   countlist, &
!   nrec)

   type(Molecule), intent(inout) :: mol0, mol1
   integer :: ntype
   integer, dimension(mol0%natom) :: typeaggs
   integer :: nequiv0, nequiv1
   integer, dimension(mol0%natom) :: equivaggs0
   integer, dimension(mol1%natom) :: equivaggs1
!   integer, intent(in) :: ntype
!   integer, dimension(:), intent(in) :: typeaggs
!   integer, intent(in) :: nequiv0, nequiv1
!   integer, dimension(:), intent(in) :: equivaggs0, equivaggs1
   integer, intent(out) :: maplist(:, :)
   integer, intent(out) :: countlist(:)
   integer, intent(out) :: nrec

   integer :: natom
   real(wp), dimension(:), allocatable :: weights
   real(wp), dimension(:, :), allocatable :: coords0, coords1
   logical, dimension(:, :), allocatable :: adjmat0, adjmat1
   logical :: visited, overflow
   integer :: nfrag0, nfrag1
   integer irec, mrec, ntrial, nstep, steps
   integer :: adjd, recadjd(maxrec)
   integer :: coonums0(mol0%natom), neighbors0(maxcoord, mol0%natom)
   integer :: coonums1(mol1%natom), neighbors1(maxcoord, mol1%natom)
   integer, dimension(mol0%natom, maxlevel) :: nadjmna0
   integer, dimension(mol1%natom, maxlevel) :: nadjmna1
   integer, dimension(maxcoord, mol0%natom, maxlevel) :: adjmnalen0
   integer, dimension(maxcoord, mol1%natom, maxlevel) :: adjmnalen1
   integer, dimension(maxcoord, mol0%natom, maxlevel) :: adjmnalist0
   integer, dimension(maxcoord, mol1%natom, maxlevel) :: adjmnalist1
   integer, dimension(mol0%natom) :: neqvnei0
   integer, dimension(mol1%natom) :: neqvnei1
   integer, dimension(maxcoord, mol0%natom) :: neieqvlens0
   integer, dimension(maxcoord, mol1%natom) :: neieqvlens1
   integer, dimension(mol0%natom) :: fragroot0
   integer, dimension(mol1%natom) :: fragroot1
   integer :: equivmat(mol0%natom, mol1%natom)
   integer, dimension(mol0%natom) :: mapping, newmapping
   real(wp) :: rmsd, totalrot
   real(wp), dimension(4) :: rotquat, prodquat
   real(wp), dimension(maxrec) :: recrmsd, avgsteps, avgtotalrot, avgrealrot
   real(wp) :: biasmat(mol0%natom, mol1%natom)
   real(wp) :: workcoords1(3, mol1%natom)

   integer i
!   integer h, i, n, offset
!   real(wp), dimension(3, natom) :: randcoords0, randcoords1
!   integer :: votes(natom, natom)

!  Temporary variable copies

   natom = mol0%natom
   weights = mol0%get_weights()
   coords0 = mol0%get_coords()
   coords1 = mol1%get_coords()
   adjmat0 = mol0%get_adjmat()
   adjmat1 = mol1%get_adjmat()

   ntype = mol0%get_ntype()
   typeaggs = mol0%get_typeaggs()

   coonums0 = mol0%get_coonums()
   coonums1 = mol1%get_coonums()
   neighbors0 = mol0%get_neighbors()
   neighbors1 = mol1%get_neighbors()

   nequiv0 = mol0%get_nequiv()
   nequiv1 = mol1%get_nequiv()
   equivaggs0 = mol0%get_equivaggs()
   equivaggs1 = mol1%get_equivaggs()

   neqvnei0(:) = mol0%atoms(:)%nneieqv
   do i = 1, mol0%natom
      neieqvlens0(:, i) = mol0%atoms(i)%neieqvlens(:)
   end do
   neqvnei1(:) = mol1%atoms(:)%nneieqv
   do i = 1, mol1%natom
      neieqvlens1(:, i) = mol1%atoms(i)%neieqvlens(:)
   end do

   ! Detect fragments and root atoms

!CZGC: new call
!   call findmolfrags(mol0, coonums0, neighbors0, nfrag0, fragroot0)
!   call findmolfrags(mol1, coonums1, neighbors1, nfrag1, fragroot1)
!CZGC: optionally
!   call findmolfrags(mol0, nfrag0, fragroot0)
!   call findmolfrags(mol1, nfrag1, fragroot1)
!CZGC: optionally 2 (nota: sólo pasan todos los tests si se llama findmolfrags
!                          aquí y no al inicio en fileio.f90)
   call findmolfrags(mol0)
   call findmolfrags(mol1)
   nfrag0 = mol0%nfrag
   fragroot0 = mol0%get_fragroot()
   nfrag1 = mol1%nfrag
   fragroot1 = mol1%get_fragroot()
!CZGC: old call
!   call findmolfrags(natom, coonums0, neighbors0, ntype, typeaggs, nequiv0, equivaggs0, nfrag0, fragroot0)
!   call findmolfrags(natom, coonums1, neighbors1, ntype, typeaggs, nequiv1, equivaggs1, nfrag1, fragroot1)

   ! Calculate MNA equivalence matrix

!CZGC: new call
   call calcequivmat(mol0, mol1, nadjmna0, adjmnalen0, adjmnalist0, nadjmna1, adjmnalen1, adjmnalist1, equivmat)
!CZGC: old call
!   call calcequivmat(natom, ntype, typeaggs, coonums0, neighbors0, nadjmna0, adjmnalen0, adjmnalist0, &
!      coonums1, neighbors1, nadjmna1, adjmnalen1, adjmnalist1, equivmat)

   ! Print equivalence matrix

!   offset = 0
!   do h = 1, ntype
!      print *, h
!      do i = offset + 1, offset + typeaggs(h)
!         print '(100(1x, i3))', equivmat(offset+1:offset+typeaggs(h), i)
!      end do
!      offset = offset + typeaggs(h)
!      print *
!   end do

   ! Calculate bias matrix


!CZGC: new call
!   call setcrossbias(mol0, mol1, equivmat, biasmat)
!CZGC: old call
   call setcrossbias(natom, ntype, typeaggs, coords0, coords1, equivmat, biasmat)

   ! Print biases

!   offset = 0
!   do h = 1, blockcount0
!      do i = offset+1, offset+typeaggs(h)
!         write (stdout, '(i0,":")', advance='no') i
!         do j = offset+1, offset+typeaggs(h)
!            if (biasmat(i, j) == 0) then
!               write (stdout, '(1x,i0)', advance='no') j
!            end if
!         end do
!         print *
!      end do
!      offset = offset + typeaggs(h)
!   end do

   ! Calculate assignment frequencies

!   votes(:, :) = 0
!   do n = 1, 100
!      do i = 1, natom
!         randcoords0(:, i) = randarray(3)
!         randcoords1(:, i) = randarray(3)
!      end do
!      call mapatoms(natom, randcoords0, randcoords1, ntype, typeaggs, biasmat, mapping)
!      print *, adjacencydiff(natom, adjmat0, adjmat1, mapping)
!      do i = 1, natom
!         votes(mapping(i), i) = votes(mapping(i), i) + 1
!      end do
!   end do
!   print *
!   offset = 0
!   do h = 1, ntype
!      print *, h
!      do i = offset + 1, offset + typeaggs(h)
!         print '(100(1x, i3))', votes(offset+1:offset+typeaggs(h), i)
!      end do
!      offset = offset + typeaggs(h)
!      print *
!   end do

   ! Initialize loop variables

   nrec = 0
   nstep = 0
   ntrial = 0
   countlist(1) = 0
   overflow = .false.

   ! Loop for map searching

   do while (countlist(1) < maxcount .and. (.not. trial_flag .or. ntrial < maxtrials))

      ntrial = ntrial + 1

      ! Work with a copy of coords1

      workcoords1 = coords1

      ! Aply a random rotation to workcoords1

      call rotate(natom, workcoords1, genrotquat(randarray(3)))

      ! Minimize the euclidean distance

      call mapatoms(natom, ntype, typeaggs, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
         nadjmna1, adjmnalen1, adjmnalist1, workcoords1, weights, equivmat, biasmat, mapping)
      rotquat = leastrotquat(natom, weights, coords0, workcoords1, mapping)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate(natom, workcoords1, rotquat)

      steps = 1

      do while (iter_flag)
         call mapatoms(natom, ntype, typeaggs, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
            nadjmna1, adjmnalen1, adjmnalist1, workcoords1, weights, equivmat, biasmat, newmapping)
         if (all(newmapping == mapping)) exit
         rotquat = leastrotquat(natom, weights, coords0, workcoords1, newmapping)
         prodquat = quatmul(rotquat, prodquat)
         call rotate(natom, workcoords1, rotquat)
         totalrot = totalrot + rotangle(rotquat)
         mapping = newmapping
         steps = steps + 1
      end do

      nstep = nstep + steps

      if (back_flag) then

         call minadjdiff(natom, weights, ntype, typeaggs, coords0, coonums0, neighbors0, adjmat0, nequiv0, &
            equivaggs0, workcoords1, coonums1, neighbors1, adjmat1, nequiv1, equivaggs1, mapping, nfrag0, fragroot0)

         call eqvatomperm(natom, weights, coords0, adjmat0, neighbors0, nequiv0, equivaggs0, &
            neqvnei0, neieqvlens0, workcoords1, adjmat1, mapping, nfrag0, fragroot0)

      end if

      adjd = adjacencydiff(natom, adjmat0, adjmat1, mapping)
      rmsd = sqrt(leastsquaredist(natom, weights, coords0, coords1, mapping)/sum(weights))

      ! Check for new best maplist

      visited = .false.

      do irec = 1, nrec
         if (all(mapping == maplist(:, irec))) then
            countlist(irec) = countlist(irec) + 1
            avgsteps(irec) = avgsteps(irec) + (steps - avgsteps(irec))/countlist(irec)
            avgrealrot(irec) = avgrealrot(irec) + (rotangle(prodquat) - avgrealrot(irec))/countlist(irec)
            avgtotalrot(irec) = avgtotalrot(irec) + (totalrot - avgtotalrot(irec))/countlist(irec)
            visited = .true.
            exit
         end if
      end do

      if (.not. visited) then
         mrec = nrec + 1
         do irec = nrec, 1, -1
            if (adjd < recadjd(irec) .or. (adjd == recadjd(irec) .and. rmsd < recrmsd(irec))) then
               mrec = irec
            else
               exit
            end if
         end do
         if (nrec < maxrec) then
            nrec = nrec + 1
         else
            overflow = .true.
         end if
         if (mrec <= maxrec) then
            do irec = nrec, mrec + 1, -1
               countlist(irec) = countlist(irec - 1)
               recrmsd(irec) = recrmsd(irec - 1)
               recadjd(irec) = recadjd(irec - 1)
               avgsteps(irec) = avgsteps(irec - 1)
               avgrealrot(irec) = avgrealrot(irec - 1)
               avgtotalrot(irec) = avgtotalrot(irec - 1)
               maplist(:, irec) = maplist(:, irec - 1)
            end do
            countlist(mrec) = 1
            recrmsd(mrec) = rmsd
            recadjd(mrec) = adjd
            avgsteps(mrec) = steps
            avgrealrot(mrec) = rotangle(prodquat)
            avgtotalrot(mrec) = totalrot
            maplist(:, mrec) = mapping
         end if
      end if

   end do

   if (stats_flag) then
      call print_stats(nrec, countlist, avgsteps, avgtotalrot, avgrealrot, recadjd, recrmsd)
      call print_final_stats(overflow, maxrec, nrec, ntrial, nstep)
   end if

end subroutine

subroutine find_reactive_sites(mol0, mol1, mapping)

   type(Molecule), intent(inout) :: mol0, mol1
   integer, dimension(:), intent(in) :: mapping

   integer :: natom
   integer :: h, i, j, k, m, n
   integer, dimension(:), allocatable :: znums0, znums1
   integer, dimension(:), allocatable :: typeidcs0, typeidcs1
   real(wp), dimension(:, :), allocatable :: coords0, coords1
   logical, dimension(:, :), allocatable :: adjmat0, adjmat1
   type(Block), dimension(:), allocatable :: blocks0, blocks1
   real(wp) :: rotquat(4)

   ! Align coordinates

   rotquat = leastrotquat(mol0%natom, mol0%get_weights(), mol0%get_coords(), mol1%get_coords(), mapping)
   call mol1%rotate_coords(rotquat)

   ! Initialization

   coords0 = mol0%get_coords()
   coords1 = mol1%get_coords()
   adjmat0 = mol0%get_adjmat()
   adjmat1 = mol1%get_adjmat()
   blocks0 = mol0%get_blocks()
   blocks1 = mol1%get_blocks()
   typeidcs0 = mol0%get_typeidcs()
   typeidcs1 = mol1%get_typeidcs()

   ! Remove reactive bonds

   do i = 1, mol0%natom
      m = typeidcs0(i)
      do j = i + 1, mol0%natom
         n = typeidcs0(j)
         if (adjmat0(i, j) .neqv. adjmat1(mapping(i), mapping(j))) then
            call remove_reactive_bond(i, j, mol0, mol1, mapping)
            do k = 1, mol0%typeaggs(m)
               if (sum((coords0(:, i) - coords1(:, mapping(blocks1(m)%atomidx(k))))**2) < 2.0 &
                  .or. sum((coords0(:, blocks0(m)%atomidx(k)) - coords1(:, mapping(i)))**2) < 2.0 &
               ) then
                  call remove_reactive_bond(blocks0(m)%atomidx(k), j, mol0, mol1, mapping)
               end if
            end do
            do k = 1, mol0%typeaggs(n)
               if (sum((coords0(:, j) - coords1(:, mapping(blocks1(n)%atomidx(k))))**2) < 2.0 &
                  .or. sum((coords0(:, blocks0(n)%atomidx(k)) - coords1(:, mapping(j)))**2) < 2.0 &
               ) then
                  call remove_reactive_bond(i, blocks1(n)%atomidx(k), mol0, mol1, mapping)
               end if
            end do
         end if
      end do
   end do

end subroutine

subroutine remove_reactive_bond(i, j, mol0, mol1, mapping)
! Purpose: Remove reactive bonds
   integer, intent(in) :: i, j
   integer, dimension(:), intent(in) :: mapping
   type(Molecule), intent(inout) :: mol0, mol1

   integer :: k
   integer, allocatable, dimension(:) :: znums0, znums1
   integer :: coonums0(mol0%natom), neighbors0(maxcoord, mol0%natom)
   integer :: coonums1(mol1%natom), neighbors1(maxcoord, mol1%natom)

   coonums0 = mol0%get_coonums()
   coonums1 = mol1%get_coonums()
   neighbors0 = mol0%get_neighbors()
   neighbors1 = mol1%get_neighbors()
   znums0 = mol0%get_znums()
   znums1 = mol1%get_znums()

   if (mol0%bonded(i, j)) then
      call mol0%remove_bond(i, j)
   end if

   if (mol1%bonded(mapping(i), mapping(j))) then
      call mol1%remove_bond(mapping(i), mapping(j))
   end if

   if (znums0(i) == 1) then
      do k = 1, coonums0(i)
         if (any([7, 8] == znums0(neighbors0(k, i)))) then
            call mol0%remove_bond(i, neighbors0(k, i))
         end if
      end do
   end if

   if (znums1(i) == 1) then
      do k = 1, coonums1(i)
         if (any([7, 8] == znums1(neighbors1(k, i)))) then
            if (mol1%bonded(mapping(i), mapping(neighbors1(k, i)))) then
               call mol1%remove_bond(mapping(i), mapping(neighbors1(k, i)))
            end if
         end if
      end do
   end if

   if (znums0(j) == 1) then
      do k = 1, coonums0(j)
         if (any([7, 8] == znums0(neighbors0(k, j)))) then
            call mol0%remove_bond(j, neighbors0(k, j))
         end if
      end do
   end if

   if (znums1(j) == 1) then
      do k = 1, coonums1(j)
         if (any([7, 8] == znums1(neighbors1(k, j)))) then
            if (mol1%bonded(mapping(j), mapping(neighbors1(k, j)))) then
               call mol1%remove_bond(mapping(j), mapping(neighbors1(k, j)))
            end if
         end if
      end do
   end if

end subroutine

end module
