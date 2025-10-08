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

module registration
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
public init_rmsd_registry
public init_adjd_registry
public init_dual_registry

type :: record_t
   integer :: count
   integer :: adjd                                      ! Used when use_adjacency=TRUE
   real(rk) :: permdist                                 ! Used when use_position=TRUE
   real(rk) :: rotation(4)
   real(rk) :: aver_steps                               ! Used when use_position=TRUE
   integer, dimension(:), allocatable :: atomperm
end type

type :: registry_t
   logical :: use_position                              ! Enable position-based processing (steps, permdist, rotation)
   logical :: use_adjacency                             ! Enable adjacency-based processing (adjacency differences)
   logical :: overflow
   integer :: total_steps                               ! Used when use_position=TRUE
   integer :: num_trials
   integer :: occ_records
   type(record_t), dimension(:), allocatable :: records
end type

contains

subroutine init_rmsd_registry(registry, num_records)
   type(registry_t), intent(inout) :: registry
   integer, intent(in) :: num_records

   if (num_records < 1) then
      error stop 'num_records < 1'
   end if

   allocate (registry%records(num_records))
   registry%use_position = .true.
   registry%use_adjacency = .false.
   registry%occ_records = 0
   registry%num_trials = 0
   registry%total_steps = 0
   registry%overflow = .false.
   registry%records%count = 0
   registry%records%permdist = huge(rk)
end subroutine

subroutine init_adjd_registry(registry, num_records)
   type(registry_t), intent(inout) :: registry
   integer, intent(in) :: num_records

   if (num_records < 1) then
      error stop 'num_records < 1'
   end if

   allocate (registry%records(num_records))
   registry%use_position = .false.
   registry%use_adjacency = .true.
   registry%occ_records = 0
   registry%num_trials = 0
   registry%overflow = .false.
   registry%records%count = 0
   registry%records%adjd = huge(ik)
end subroutine

subroutine init_dual_registry(registry, num_records)
   type(registry_t), intent(inout) :: registry
   integer, intent(in) :: num_records

   if (num_records < 1) then
      error stop 'num_records < 1'
   end if

   allocate (registry%records(num_records))
   registry%use_position = .true.
   registry%use_adjacency = .true.
   registry%occ_records = 0
   registry%num_trials = 0
   registry%total_steps = 0
   registry%overflow = .false.
   registry%records%count = 0
   registry%records%adjd = huge(ik)
   registry%records%permdist = huge(rk)
end subroutine

subroutine insert_record(registry, atomperm, num_steps, adjd, permdist, rotation)
   type(registry_t), target, intent(inout) :: registry
   integer, dimension(:), intent(in) :: atomperm
   integer, intent(in) :: num_steps
   integer, intent(in), optional :: adjd
   real(rk), intent(in), optional :: permdist
   real(rk), intent(in), optional :: rotation(4)
   ! Local variables
   type(record_t), pointer :: record
   logical :: should_insert
   integer :: i, j

   ! Validate required parameters based on flags
   if (.not. registry%use_position .and. .not. registry%use_adjacency) then
      error stop 'Registry not properly initialized - no processing mode enabled'
   end if

   ! Handle position-based processing
   if (registry%use_position) then
      if (.not. present(permdist)) error stop 'permdist required when use_position=TRUE'
      if (.not. present(rotation)) error stop 'rotation required when use_position=TRUE'
   end if

   ! Handle adjacency-based processing
   if (registry%use_adjacency) then
      if (.not. present(adjd)) error stop 'adjd required when use_adjacency=TRUE'
   end if

   registry%num_trials = registry%num_trials + 1
   registry%total_steps = registry%total_steps + num_steps

   ! Check for existing records to update
   do i = 1, registry%occ_records
      record => registry%records(i)
      if (all(atomperm == record%atomperm)) then
         record%count = record%count + 1
         record%aver_steps = record%aver_steps + (num_steps - record%aver_steps) / record%count
         return
      end if
   end do

   ! Find insertion point and insert new record
   do i = 1, size(registry%records)
      record => registry%records(i)

      ! Determine insertion criteria based on mode
      if (registry%use_adjacency .and. registry%use_position) then
         should_insert = (adjd < record%adjd .or. (adjd == record%adjd .and. permdist < record%permdist))
      else if (registry%use_position) then
         should_insert = (permdist < record%permdist)
      else
         should_insert = (adjd < record%adjd)
      end if

      if (should_insert) then
         ! Shift records to make room
         do j = size(registry%records), i + 1, -1
            registry%records(j) = registry%records(j - 1)
         end do

         ! Initialize new record
         record%atomperm = atomperm
         record%count = 1

         ! Set fields based on enabled processing modes
         if (registry%use_position) then
            record%permdist = permdist
            record%rotation = rotation
            record%aver_steps = num_steps
         end if

         if (registry%use_adjacency) then
            record%adjd = adjd
         end if

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
   ! Local variables
   type(record_t) :: record
   character(49) :: line
   integer :: i

   ! Validate registry state
   if (.not. registry%use_position .and. .not. registry%use_adjacency) then
      error stop 'Registry not properly initialized - no processing mode enabled'
   end if

   if (registry%use_position .and. .not. registry%use_adjacency) then
      ! Position only
      line = repeat('-', 42)
      write (stderr, '(2x,a,4x,a,5x,a,5x,a,7x,a)') '#', 'Count', 'Steps', 'Rotθ', 'Δxyz'
      write (stderr, '(a)') line(1:42)
      do i = 1, registry%occ_records
         record = registry%records(i)
         write (stderr, '(i3,4x,i4,4x,f5.1,5x,f5.1,4x,f8.4)') &
            i, record%count, record%aver_steps, angle(record%rotation), record%permdist
      end do
      write (stderr, '(a)') line(1:42)
   else if (registry%use_position .and. registry%use_adjacency) then
      ! Both position and adjacency
      line = repeat('-', 49)
      write (stderr, '(2x,a,4x,a,5x,a,4x,a,5x,a,6x,a)') '#', 'Count', 'Steps', 'Rotθ', 'Δadj', 'Δxyz'
      write (stderr, '(a)') line
      do i = 1, registry%occ_records
         record = registry%records(i)
         write (stderr, '(i3,4x,i4,4x,f5.1,5x,f5.1,3x,i4,4x,f8.4)') &
            i, record%count, record%aver_steps, angle(record%rotation), record%adjd, record%permdist
      end do
      write (stderr, '(a)') line
   else if (registry%use_adjacency .and. .not. registry%use_position) then
      ! Adjacency only
      line = repeat('-', 25)
      write (stderr, '(2x,a,4x,a,4x,a)') '#', 'Count', 'Δadj'
      write (stderr, '(a)') line(1:25)
      do i = 1, registry%occ_records
         record = registry%records(i)
         write (stderr, '(i3,4x,i4,4x,i4)') i, record%count, record%adjd
      end do
      write (stderr, '(a)') line(1:25)
   end if

   write (stderr, *)
   write (stderr, '(a,1x,i0)') 'Random trials =', registry%num_trials
   if (registry%use_position) then
      write (stderr, '(a,1x,i0)') 'Minimization steps =', registry%total_steps
   end if
   if (registry%overflow) then
      write (stderr, '(a,1x,i0)') 'Visited local minima >', registry%occ_records
   else
      write (stderr, '(a,1x,i0)') 'Visited local minima =', registry%occ_records
   end if
   write (stderr, *)

   flush(stderr)
end subroutine

end module
