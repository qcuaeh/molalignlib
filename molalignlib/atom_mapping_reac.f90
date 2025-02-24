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

module atom_mapping_reac
use parameters
use globals
use random
use molecule
use chemdata
use rotation
use rigid_body
use adjacency
use alignment
use lap_driver
use bipartition
use biasing
use printing
use registry
use tracking
use backtracking
use reactivity

implicit none

contains

subroutine remap_reactive_bonds( mol1, mol2, eltypes, mnatypes, results)
   type(mol_type), intent(inout) :: mol1, mol2
   type(bipartition_type), intent(in) :: eltypes, mnatypes
   type(registry_type), target, intent(out) :: results

   ! Local variables
   integer :: num_atoms1
   integer :: num_trials, num_steps
   integer, pointer :: lead_count
   real(rk) :: eigquat(4), totquat(4)
   integer, dimension(:), allocatable :: atomperm, auxperm
   real(rk), dimension(:,:), allocatable :: coords1, coords2
   type(intlist_type), dimension(:), allocatable :: molfrags1, molfrags2
   type(intmatrix_type), allocatable :: mnadiffs(:)
   integer :: adjd

   num_atoms1 = size(mol1%atoms)
   coords1 = mol1%get_weightcoords()
   coords2 = mol2%get_weightcoords()
   call results%initialize(max_records)

   allocate (atomperm(num_atoms1))
   allocate (auxperm(num_atoms1))

   ! Find molecular fragments
   call find_molfrags( mol1, eltypes%partition1(), molfrags1)
   call find_molfrags( mol2, eltypes%partition2(), molfrags2)

   ! Reflect atoms
   if (mirror_flag) then
      call reflect_coords( coords2)
   end if

   ! Translate atoms to their centroids
   call translate_coords( coords1, -centroid(coords1))
   call translate_coords( coords2, -centroid(coords2))

   ! Find unfeasible assignments
   call bias_procedure( mol1, mol2, eltypes, mnadiffs)

   ! Initialize random number generator
   call random_initialize()

   ! Optimize atom permutation

   num_trials = 0
   lead_count => results%records(1)%count

   do while (lead_count < max_count .and. num_trials < max_trials)

      num_trials = num_trials + 1

      ! Aply a random rotation to coords2
      call rotate_coords(coords2, randrotquat())

      ! Assign atoms with current orientation
      call assign_atoms_biased(eltypes, coords1, coords2, mnadiffs, atomperm)
      totquat = leasteigquat(atomperm, coords1, coords2)
      call rotate_coords(coords2, totquat)
      num_steps = 1

      do while (iter_flag)
         call assign_atoms_biased(eltypes, coords1, coords2, mnadiffs, auxperm)
         if (all(auxperm == atomperm)) exit
         atomperm = auxperm
         eigquat = leasteigquat(atomperm, coords1, coords2)
         call rotate_coords(coords2, eigquat)
         totquat = quatmul(eigquat, totquat)
         num_steps = num_steps + 1
      end do

      call minadjdiff( eltypes, mnatypes, molfrags1, mol1, mol2, coords1, coords2, atomperm)

      ! Update results
      adjd = adjacencydiff(atomperm, mol1%adjmat, mol2%adjmat)
      call results%push_adjd(atomperm, num_steps, angle(totquat), adjd)

   end do

   results%num_trials = num_trials

   ! Remove bonds from reactive sites and reoptimize assignment
   atomperm = results%records(1)%atomperm
!   write (stderr, '(i0)') adjacencydiff( atomperm, mol1%adjmat, mol2%adjmat)
   call remove_reactive_bonds( mol1, mol2, molfrags1, molfrags2, mnatypes, atomperm)
!   write (stderr, '(i0)') adjacencydiff( atomperm, mol1%adjmat, mol2%adjmat)

end subroutine

end module
