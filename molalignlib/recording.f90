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
implicit none
private

public record_t
public registry_t
public insert_record
public print_records
public init_adjacency_grouped_registry
public init_permutation_grouped_registry

type :: record_t
   integer :: count
   integer :: permdiff                                      ! Used when track_adjacency=TRUE
   real(rk) :: permdist                                 ! Used when track_position=TRUE
   real(rk) :: rotation(4)
   real(rk) :: steps                               ! Used when track_position=TRUE
   integer, dimension(:), allocatable :: atomperm
end type

type :: registry_t
   logical :: track_position                            ! Track position-based metrics (steps, permdist, rotation)
   logical :: track_adjacency                           ! Track adjacency-based metrics (adjacency differences)
   logical :: group_by_adjacency                        ! Group records by adjacency difference instead of permutation
   logical :: overflow
   integer :: total_steps                               ! Used when track_position=TRUE
   integer :: num_trials
   integer :: occ_records
   type(record_t), dimension(:), allocatable :: records
end type

contains

subroutine init_adjacency_grouped_registry(registry, num_records)
   ! Groups by adjacency difference and keeps optimal permutation
   ! that minimizes squared distance sum for each adjacency value
   type(registry_t), intent(inout) :: registry
   integer, intent(in) :: num_records

   if (num_records < 1) then
      error stop 'num_records < 1'
   end if

   allocate (registry%records(num_records))
   registry%track_position = .true.                     ! Track position to find optimal permutation
   registry%track_adjacency = .true.                    ! Track adjacency differences
   registry%group_by_adjacency = .true.                 ! Group by adjacency, not by permutation
   registry%occ_records = 0
   registry%num_trials = 0
   registry%total_steps = 0
   registry%overflow = .false.
   registry%records%count = 0
   registry%records%permdiff = huge(ik)
   registry%records%permdist = huge(rk)
end subroutine

subroutine init_permutation_grouped_registry(registry, num_records)
   ! Groups by permutation, tracks both adjacency and position metrics
   type(registry_t), intent(inout) :: registry
   integer, intent(in) :: num_records

   if (num_records < 1) then
      error stop 'num_records < 1'
   end if

   allocate (registry%records(num_records))
   registry%track_position = .true.
   registry%track_adjacency = .true.
   registry%group_by_adjacency = .false.                ! Group by permutation
   registry%occ_records = 0
   registry%num_trials = 0
   registry%total_steps = 0
   registry%overflow = .false.
   registry%records%count = 0
   registry%records%permdiff = huge(ik)
   registry%records%permdist = huge(rk)
end subroutine

subroutine insert_record(registry, atomperm, permdiff, permdist, steps, rotation)
   ! Handles both adjacency-grouped and permutation-grouped modes
   type(registry_t), target, intent(inout) :: registry
   integer, dimension(:), intent(in) :: atomperm
   integer, intent(in) :: steps
   integer, intent(in) :: permdiff
   real(rk), intent(in) :: permdist
   real(rk), intent(in) :: rotation(4)
   ! Local variables
   type(record_t), pointer :: record
   integer :: i, j

   registry%num_trials = registry%num_trials + 1
   registry%total_steps = registry%total_steps + steps

   if (registry%group_by_adjacency) then
      ! Adjacency-grouped mode: Group by adjacency difference, keep optimal permutation

      ! Check for existing record with same adjacency difference
      do i = 1, registry%occ_records
         record => registry%records(i)

         if (permdiff == record%permdiff) then
            ! Found matching adjacency - increment count
            record%count = record%count + 1

            ! Keep the permutation with minimum permdist
            if (permdist < record%permdist) then
               ! Found better permutation for this adjacency difference
               record%atomperm = atomperm
               record%permdist = permdist
               record%rotation = rotation
               record%steps = (record%steps * (record%count - 1) + steps) / record%count
            else
               ! Keep existing permutation, but update average steps
               record%steps = record%steps + (steps - record%steps) / record%count
            end if
            return
         end if
      end do
   else
      ! Permutation-grouped mode: Group by permutation

      ! Check for existing record with same permutation
      do i = 1, registry%occ_records
         record => registry%records(i)
         if (all(atomperm == record%atomperm)) then
            record%count = record%count + 1
            record%steps = record%steps + (steps - record%steps) / record%count
            return
         end if
      end do
   end if

   ! Find insertion point: sort by adjacency first, then by permdist
   do i = 1, size(registry%records)
      record => registry%records(i)
      if (permdiff < record%permdiff .or. (permdiff == record%permdiff .and. permdist < record%permdist)) then
         ! Shift records to make room
         do j = size(registry%records), i + 1, -1
            registry%records(j) = registry%records(j - 1)
         end do

         ! Initialize new record
         record%atomperm = atomperm
         record%count = 1
         record%permdiff = permdiff
         record%permdist = permdist
         record%rotation = rotation
         record%steps = steps

         exit
      end if
   end do

   ! Update record count and overflow status
   if (.not. registry%overflow) then
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
   integer :: i

   if (registry%track_adjacency) then
      ! Both adjacency-grouped and permutation-grouped modes
      line = repeat('-', 49)
      write (stderr, '(2x,a,4x,a,5x,a,4x,a,5x,a,6x,a)') '#', 'Count', 'Steps', 'RotΘ', 'Δadj', 'Δxyz'
      write (stderr, '(a)') line
      do i = 1, registry%occ_records
         record = registry%records(i)
         write (stderr, '(i3,4x,i4,4x,f5.1,5x,f5.1,3x,i4,4x,f8.4)') &
            i, record%count, record%steps, angle(record%rotation), record%permdiff, record%permdist
      end do
      write (stderr, '(a)') line
   else
      ! Position only mode (should not occur with current init procedures)
      line = repeat('-', 42)
      write (stderr, '(2x,a,4x,a,5x,a,5x,a,7x,a)') '#', 'Count', 'Steps', 'RotΘ', 'Δxyz'
      write (stderr, '(a)') line(1:42)
      do i = 1, registry%occ_records
         record = registry%records(i)
         write (stderr, '(i3,4x,i4,4x,f5.1,5x,f5.1,4x,f8.4)') &
            i, record%count, record%steps, angle(record%rotation), record%permdist
      end do
      write (stderr, '(a)') line(1:42)
   end if

   write (stderr, *)
   write (stderr, '(a,1x,i0)') 'Random trials =', registry%num_trials
   write (stderr, '(a,1x,i0)') 'Minimization steps =', registry%total_steps

   if (registry%overflow) then
      write (stderr, '(a,1x,i0)') 'Visited local minima >', registry%occ_records
   else
      write (stderr, '(a,1x,i0)') 'Visited local minima =', registry%occ_records
   end if

   if (registry%group_by_adjacency) then
      write (stderr, '(a)') 'Note: Count tracks repetitions of same Δadj'
   end if

   write (stderr, *)
end subroutine

end module
