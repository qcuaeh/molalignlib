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
use parameters
use options
use derived_types
use sorting
use utils
use permutation
use molecule
use lcrs_tree
use atom_types
use atom_mnas
use spatial_transforms
implicit none
contains

subroutine compute_mna_biases(atoms1, atoms2, atomtypes, biases)
! Iteratively compute MNA types
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_t), intent(in) :: atomtypes
   type(int_matrix), dimension(:), allocatable, intent(out) :: biases
   ! Local variables
   type(assigntree_node_t), pointer :: mnachain
   integer :: h, i, j, iatom, jatom
   integer :: num_splits
!   integer :: link_idx

   allocate(biases(atomtypes%num_parts))

   do h = 1, atomtypes%num_parts
      allocate(biases(h)%ee(atomtypes%parts(h)%num_items1, atomtypes%parts(h)%num_items2))
      biases(h)%ee = 0
   end do

   ! Initialize MNA chain with element types
   mnachain => chain_from_partition(atomtypes)

!   link_idx = 0
   do

!      write(stderr, *)
!      write(stderr, '(a)') repeat('-- link_idx '//str(link_idx)//' --', 6)

      ! Call compute_mna_partition and get the number of splits
      call compute_mna_partition(atoms1, atoms2, mnachain, num_splits)

      ! Exit loop if no splits occurred in the last iteration
      if (num_splits == 0) exit

      ! Update biases with MNAs at current level
      do h = 1, atomtypes%num_parts
         do j = 1, atomtypes%parts(h)%num_items2
            jatom = atomtypes%parts(h)%items2(j)
            do i = 1, atomtypes%parts(h)%num_items1
               iatom = atomtypes%parts(h)%items1(i)
               if (associated( &
                  mnachain%last_link%itemdir1(iatom)%ptr, &
                  mnachain%last_link%itemdir2(jatom)%ptr) &
               ) then
                  biases(h)%ee(i, j) = biases(h)%ee(i, j) + 1
               end if
            end do
         end do
      end do

!      link_idx = link_idx + 1
   end do

!   do h = 1, atomtypes%num_parts
!      write(stderr, *)
!      do j = 1, atomtypes%parts(h)%num_items2
!         write(stderr, '(*(i2))') biases(h)%ee(:atomtypes%parts(h)%num_items1, j)
!      end do
!   end do

   call delete_chain(mnachain)  ! Cleanup
end subroutine


! tests all permutations of the neighbors of a pair of atoms
subroutine build_minbiases(atomtypes, atoms1, atoms2, biases, minbiases)
   type(partition_t), intent(in) :: atomtypes
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(int_matrix), dimension(:), intent(in) :: biases
   type(real_matrix), dimension(:), allocatable, intent(out) :: minbiases

   integer :: h, ha, i, j, iatom, jatom, maxbias
   integer, dimension(:), allocatable :: adjlistat1, adjlistat2, atomperm
   integer, dimension(:), allocatable :: itemsat1, itemsat2
   real(rk), dimension(:,:), allocatable :: coords1, coords2

   type(partition_t) :: adjeltypes

   logical :: more
   integer :: rank, n
   real(rk) :: rmsd, min_rmsd, center1(3), center2(3)

   allocate(minbiases(atomtypes%num_parts))
   coords1 = get_coords( atoms1)
   coords2 = get_coords( atoms2)


   do h = 1, atomtypes%num_parts   ! run over all eltype partitions
!write (stderr,*) "eltype: ", h

      ! initialize minbiases matrix elements
      allocate(minbiases(h)%ee(atomtypes%parts(h)%num_items1, atomtypes%parts(h)%num_items2))
      minbiases(h)%ee = 0.5
      maxbias = maxval(biases(h)%ee)
!write (stderr,*) "maxbias: ", maxbias
      do i = 1, atomtypes%parts(h)%num_items1   ! run over atoms in atoms1
         iatom = atomtypes%parts(h)%items1(i)
         center1(:) = coords1(:,iatom)
!        center1(:) = 0

         ! neighbors for iatom
         allocate(adjlistat1(size(atoms1(iatom)%adjlist)))
         adjlistat1 = atoms1(iatom)%adjlist
!write (stderr,*) "adjlist1: ", adjlistat1

         do j = 1, atomtypes%parts(h)%num_items2   ! run over atoms in atoms2
            jatom = atomtypes%parts(h)%items2(j)
            center2(:) = coords2(:,jatom)
!            center2(:) = 0

            if (biases(h)%ee(i,j) == maxbias .or. biases(h)%ee(i,j) >= 1) then   ! compatible neighbors
!*** ¿el orden de los átomos en atomtypes es el mismo que en biases?
               ! neighbors for jatom
               allocate(adjlistat2(size(atoms2(jatom)%adjlist)))
               adjlistat2 = atoms2(jatom)%adjlist
!write (stderr,*) "adjlist2: ", adjlistat2

               call collect_atomtypes(atoms1(adjlistat1), atoms2(adjlistat2), adjeltypes)   ! eltype for adjlists

               ! run over adjlist's atomtypes
               do ha = 1, adjeltypes%num_parts
!write (stderr, *) "adj part: ", ha
                  n = adjeltypes%parts(ha)%num_items1   ! same as num_items2?
                  if (n >= 2) then

                     ! process neighbors of the same eltype
                     allocate(atomperm(n), itemsat1(n), itemsat2(n))
                     itemsat1 = adjlistat1(adjeltypes%parts(ha)%items1)
                     itemsat2 = adjlistat2(adjeltypes%parts(ha)%items2)
!write (stderr, *) "itemsat1: ", itemsat1
!write (stderr, *) "itemsat2: ", itemsat2

                     ! initialize permutations
                     more = .false.
                     call perm1_next3(n, atomperm, more, rank)

                     ! initial min_rmsd value; constant initial value?
                     min_rmsd = least_total_sqdist(itemsat1, itemsat2(atomperm), coords1, coords2, center1, center2)

                     ! run over permutations of neighbos with the same eltype
                     do while (more)
                        rmsd = least_total_sqdist(itemsat1, itemsat2(atomperm), coords1, coords2, center1, center2)   ! which center is needed?º
!write (stderr, '(A,f0.4)') "rmsd: ", rmsd
                        if (rmsd < min_rmsd) min_rmsd = rmsd
                        call perm1_next3(n, atomperm, more, rank)
                     end do
!write (stderr, '(A,f0.4)') "min_rmsd: ", min_rmsd
                     deallocate (atomperm, itemsat1, itemsat2)
                  end if
               end do
               deallocate(adjlistat2)
            end if

            ! record the minimum RMSD for the i-j assignment
            minbiases(h)%ee(i,j) = min_rmsd

         end do
         deallocate(adjlistat1)
      end do
   end do

!   do h = 1, atomtypes%num_parts
!      write(stderr, *)
!      do j = 1, atomtypes%parts(h)%num_items2
!         write(stderr, '(*(f4.1,1X))') minbiases(h)%ee(:atomtypes%parts(h)%num_items1, j)
!      end do
!   end do

end subroutine

end module
