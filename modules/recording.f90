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
public insert_record_atomperm
public insert_record_moldiff
public print_records
public allocate_registry
public reset_registry

type :: record_t
   integer :: freq
   integer :: permdiff
   real(rk) :: permdist
   real(rk) :: steps
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

subroutine allocate_registry(registry, num_records)
   type(registry_t), intent(inout) :: registry
   integer, intent(in) :: num_records

   if (num_records < 1) then
      error stop 'num_records is less than 1'
   end if

   allocate (registry%records(num_records))
end subroutine

subroutine reset_registry(registry)
   type(registry_t), intent(inout) :: registry

   registry%occ_records = 0
   registry%num_trials = 0
   registry%total_steps = 0
   registry%overflow = .false.
   registry%records%freq = 0
   registry%records%steps = 0
   registry%records%permdiff = huge(ik)
   registry%records%permdist = huge(rk)
end subroutine

subroutine insert_record_atomperm(registry, atomperm1, steps, rotation, permdiff, permdist)
   ! Group by atomperm1, sort by permdist
   type(registry_t), target, intent(inout) :: registry
   integer, dimension(:), intent(in) :: atomperm1
   real(rk), intent(in) :: steps, rotation(4)
   integer, intent(in) :: permdiff
   real(rk), intent(in) :: permdist
   ! Local variables
   type(record_t), pointer :: record
   integer :: i, j, insert_pos

   registry%num_trials = registry%num_trials + 1
   registry%total_steps = registry%total_steps + steps

   ! Check for existing record with same permutation
   do i = 1, registry%occ_records
      record => registry%records(i)
      if (all(atomperm1 == record%atomperm1)) then
         record%freq = record%freq + 1
         record%steps = record%steps + (steps - record%steps)/record%freq
         return
      end if
   end do

   ! Find insertion point: sort by permdist
   insert_pos = registry%occ_records + 1
   do i = 1, registry%occ_records
      record => registry%records(i)
      if (permdiff < record%permdiff) then
         insert_pos = i
         exit
      else if (permdiff == record%permdiff) then
         if (permdist < record%permdist) then
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

      ! Initialize new record (no moldiffs for permutation grouping)
      record => registry%records(insert_pos)
      record%atomperm1 = atomperm1
      record%freq = 1
      record%permdiff = permdiff
      record%permdist = permdist
      record%rotation = rotation
      record%steps = steps

      ! Update record freq and overflow status
      if (registry%occ_records < size(registry%records)) then
         registry%occ_records = registry%occ_records + 1
      else
         registry%overflow = .true.
      end if
   end if
end subroutine

subroutine insert_record_moldiff(registry, moldiffs, atomperm1, steps, rotation, permdist)
   type(registry_t), target, intent(inout) :: registry
   integer, dimension(:,:), intent(in) :: moldiffs
   integer, dimension(:), intent(in) :: atomperm1
   integer, intent(in) :: steps
   real(rk), intent(in) :: rotation(4)
   real(rk), intent(in) :: permdist
   type(record_t), pointer :: record
   type(record_t) :: temp_record
   integer :: i, j, insert_pos, permdiff, match_pos
   logical :: bonds_match, found_match

   registry%num_trials = registry%num_trials + 1
   registry%total_steps = registry%total_steps + steps

   permdiff = size(moldiffs, 2)
   found_match = .false.
   match_pos = 0

   do i = 1, registry%occ_records
      record => registry%records(i)

      if (permdiff == record%permdiff) then
         bonds_match = all(moldiffs(1, :) == record%moldiffs(1, :)) .and. &
                       all(moldiffs(2, :) == record%moldiffs(2, :))

         if (bonds_match) then
            found_match = .true.
            match_pos = i
            record%freq = record%freq + 1
            record%steps = record%steps + (steps - record%steps)/record%freq

            if (permdist < record%permdist) then
               temp_record = record
               temp_record%atomperm1 = atomperm1
               temp_record%permdist = permdist
               temp_record%rotation = rotation

               do j = i, registry%occ_records - 1
                  registry%records(j) = registry%records(j + 1)
               end do
               registry%occ_records = registry%occ_records - 1

               insert_pos = registry%occ_records + 1
               do j = 1, registry%occ_records
                  if (permdiff < registry%records(j)%permdiff .or. &
                      (permdiff == registry%records(j)%permdiff .and. &
                       permdist < registry%records(j)%permdist)) then
                     insert_pos = j
                     exit
                  end if
               end do

               if (insert_pos <= size(registry%records)) then
                  do j = min(registry%occ_records, size(registry%records) - 1), insert_pos, -1
                     registry%records(j + 1) = registry%records(j)
                  end do

                  registry%records(insert_pos) = temp_record

                  if (registry%occ_records < size(registry%records)) then
                     registry%occ_records = registry%occ_records + 1
                  else
                     registry%overflow = .true.
                  end if
               end if
            end if
            return
         end if
      end if
   end do

   insert_pos = registry%occ_records + 1
   do i = 1, registry%occ_records
      if (permdiff < registry%records(i)%permdiff .or. &
          (permdiff == registry%records(i)%permdiff .and. &
           permdist < registry%records(i)%permdist)) then
         insert_pos = i
         exit
      end if
   end do

   if (insert_pos <= size(registry%records)) then
      do j = min(registry%occ_records, size(registry%records) - 1), insert_pos, -1
         registry%records(j + 1) = registry%records(j)
      end do

      record => registry%records(insert_pos)
      record%atomperm1 = atomperm1
      record%freq = 1
      record%moldiffs = moldiffs
      record%permdiff = permdiff
      record%permdist = permdist
      record%rotation = rotation
      record%steps = steps

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
   character(:), allocatable :: line
   integer :: i

   line = repeat('-', 39)
   write (stdout, '(2x,a,4x,a,5x,a,4x,a,5x,a,6x,a)') '#', 'Freq', 'Steps', 'Δadj', 'Δxyz'
   write (stdout, '(a)') line
   do i = 1, registry%occ_records
      record = registry%records(i)
      write (stdout, '(i3,4x,i4,4x,f5.1,3x,i4,4x,f8.4)') &
         i, record%freq, record%steps, record%permdiff, record%permdist
   end do
   write (stdout, '(a)') line

   write (stdout, *)
   write (stdout, '(a,1x,i0)') 'Random trials:', registry%num_trials
   write (stdout, '(a,1x,i0)') 'Minimization steps:', registry%total_steps

   if (registry%overflow) then
      write (stdout, '(a,1x,i0)') 'Visited local minima: >', registry%occ_records
   else
      write (stdout, '(a,1x,i0)') 'Visited local minima:', registry%occ_records
   end if

   write (stdout, *)
end subroutine

end module
