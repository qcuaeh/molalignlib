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

module recording
use parameters
use permutation
use adjacency
use euclidean
use sorting
implicit none
private

public record_t
public registry_t
public insert_record_homo
public insert_record_hetero
public print_records
public init_registry

type :: record_t
   integer :: count
   real(rk) :: steps
   real(rk) :: permdist
   real(rk) :: rotation(4)
   integer, dimension(:), allocatable :: atomperm1
   integer, dimension(:,:), allocatable :: moldiffs    ! List of differing bond pairs [2, nbonds]
end type

type :: registry_t
   logical :: overflow
   integer :: total_steps
   integer :: num_trials
   integer :: occ_records
   type(record_t), dimension(:), allocatable :: records
end type

contains

subroutine init_registry(registry, num_records)
   type(registry_t), intent(inout) :: registry
   integer, intent(in) :: num_records

   if (num_records < 1) then
      error stop 'num_records < 1'
   end if

   allocate (registry%records(num_records))
   registry%occ_records = 0
   registry%num_trials = 0
   registry%total_steps = 0
   registry%overflow = .false.
   registry%records%count = 0
   registry%records%permdist = huge(rk)
end subroutine

subroutine insert_record_homo(registry, atomperm1, permdist, steps, rotation)
   ! Group by atomperm1, sort by permdist
   type(registry_t), target, intent(inout) :: registry
   integer, dimension(:), intent(in) :: atomperm1
   integer, intent(in) :: steps
   real(rk), intent(in) :: permdist
   real(rk), intent(in) :: rotation(4)
   ! Local variables
   type(record_t), pointer :: record
   integer :: i, j, insert_pos

   registry%num_trials = registry%num_trials + 1
   registry%total_steps = registry%total_steps + steps

   ! Check for existing record with same permutation
   do i = 1, registry%occ_records
      record => registry%records(i)
      if (all(atomperm1 == record%atomperm1)) then
         record%count = record%count + 1
         record%steps = record%steps + (steps - record%steps) / record%count
         return
      end if
   end do

   ! Find insertion point: sort by permdist
   insert_pos = registry%occ_records + 1
   do i = 1, registry%occ_records
      if (permdist < registry%records(i)%permdist) then
         insert_pos = i
         exit
      end if
   end do

   ! Only insert if position is within bounds
   if (insert_pos <= size(registry%records)) then
      ! Shift records to make room (if full, last one gets dropped)
      do j = min(registry%occ_records, size(registry%records) - 1), insert_pos, -1
         registry%records(j + 1) = registry%records(j)
      end do

      ! Initialize new record (no moldiffs for permutation grouping)
      record => registry%records(insert_pos)
      record%atomperm1 = atomperm1
      record%count = 1
      record%permdist = permdist
      record%rotation = rotation
      record%steps = steps

      ! Update record count and overflow status
      if (registry%occ_records < size(registry%records)) then
         registry%occ_records = registry%occ_records + 1
      else
         registry%overflow = .true.
      end if
   end if
end subroutine

subroutine insert_record_hetero(registry, atomperm1, moldiffs, permdist, steps, rotation)
   ! Group by differing bonds, sort by permdiff then permdist
   ! moldiffs must be pre-sorted
   type(registry_t), target, intent(inout) :: registry
   integer, dimension(:), intent(in) :: atomperm1
   integer, dimension(:,:), intent(in) :: moldiffs  ! Must be sorted [2, nbonds]
   integer, intent(in) :: steps
   real(rk), intent(in) :: permdist
   real(rk), intent(in) :: rotation(4)
   ! Local variables
   type(record_t), pointer :: record
   integer :: i, j, insert_pos, num_diff
   logical :: bonds_match

   registry%num_trials = registry%num_trials + 1
   registry%total_steps = registry%total_steps + steps

   num_diff = size(moldiffs, 2)

   ! Check for existing record with same set of differing bonds
   ! Both arrays are sorted, so simple equality check works
   do i = 1, registry%occ_records
      record => registry%records(i)

      ! Check if the differing bond sets match
      if (allocated(record%moldiffs)) then
         if (num_diff == size(record%moldiffs, 2)) then
            ! Compare all bond pairs
            bonds_match = all(moldiffs(1, :) == record%moldiffs(1, :)) .and. &
                          all(moldiffs(2, :) == record%moldiffs(2, :))

            if (bonds_match) then
               ! Found matching bond difference set - increment count
               record%count = record%count + 1

               ! Keep the permutation with minimum permdist
               if (permdist < record%permdist) then
                  ! Found better permutation for this bond difference set
                  record%atomperm1 = atomperm1
                  record%permdist = permdist
                  record%rotation = rotation
                  record%steps = (record%steps * (record%count - 1) + steps) / record%count
               else
                  ! Keep existing permutation, but update average steps
                  record%steps = record%steps + (steps - record%steps) / record%count
               end if
               return
            end if
         end if
      end if
   end do

   ! Find insertion point: sort by number of differing bonds first, then by permdist
   insert_pos = registry%occ_records + 1
   do i = 1, registry%occ_records
      if (allocated(registry%records(i)%moldiffs)) then
         if (num_diff < size(registry%records(i)%moldiffs, 2) .or. &
             (num_diff == size(registry%records(i)%moldiffs, 2) .and. &
              permdist < registry%records(i)%permdist)) then
            insert_pos = i
            exit
         end if
      else
         ! Unallocated moldiffs means permutation record (shouldn't happen in practice)
         ! Treat as having 0 differing bonds
         if (num_diff < 0) then
            insert_pos = i
            exit
         end if
      end if
   end do

   ! Only insert if position is within bounds
   if (insert_pos <= size(registry%records)) then
      ! Shift records to make room (if full, last one gets dropped)
      do j = min(registry%occ_records, size(registry%records) - 1), insert_pos, -1
         registry%records(j + 1) = registry%records(j)
      end do

      ! Initialize new record (bonds already sorted)
      record => registry%records(insert_pos)
      record%atomperm1 = atomperm1
      record%count = 1
      record%moldiffs = moldiffs  ! Automatic allocation
      record%permdist = permdist
      record%rotation = rotation
      record%steps = steps

      ! Update record count and overflow status
      if (registry%occ_records < size(registry%records)) then
         registry%occ_records = registry%occ_records + 1
      else
         registry%overflow = .true.
      end if
   end if
end subroutine

subroutine print_records(registry)
   type(registry_t), intent(in) :: registry
   type(record_t) :: record
   character(49) :: line
   integer :: i, num_diff

   line = repeat('-', 49)
   write (stderr, '(2x,a,4x,a,5x,a,4x,a,5x,a,6x,a)') '#', 'Count', 'Steps', 'Rotθ', 'Δadj', 'Δxyz'
   write (stderr, '(a)') line
   do i = 1, registry%occ_records
      record = registry%records(i)

      ! Get number of differing bonds
      if (allocated(record%moldiffs)) then
         num_diff = size(record%moldiffs, 2)
      else
         num_diff = 0
      end if

      write (stderr, '(i3,4x,i4,4x,f5.1,5x,f5.1,3x,i4,4x,f8.4)') &
         i, record%count, record%steps, angle(record%rotation), num_diff, record%permdist
   end do
   write (stderr, '(a)') line

   write (stderr, *)
   write (stderr, '(a,1x,i0)') 'Random trials =', registry%num_trials
   write (stderr, '(a,1x,i0)') 'Minimization steps =', registry%total_steps

   if (registry%overflow) then
      write (stderr, '(a,1x,i0)') 'Visited local minima >', registry%occ_records
   else
      write (stderr, '(a,1x,i0)') 'Visited local minima =', registry%occ_records
   end if

   write (stderr, *)
end subroutine

end module
