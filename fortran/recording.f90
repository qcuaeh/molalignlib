! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

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
! Registry of the distinct local minima found by the stochastic searches,
! ranked by adjacency difference and then by squared distance
use parameters
use permutation
use adjacency
use euclidean
use sorting
implicit none
private

public record_t
public registry_t
public insert_record_mapping
public insert_record_moldiff
public print_records
public allocate_registry
public reset_registry

type :: record_t
   integer(ik) :: freq         ! times the minimum was found
   integer(ik) :: mapdiff      ! adjacency difference
   real(rk) :: mapdist         ! sum of squared distances
   real(rk) :: steps           ! mean number of minimization steps
   real(rk) :: rotation(4)
   integer(ik), dimension(:), allocatable :: mapping1
   integer(ik), dimension(:,:), allocatable :: moldiffs  ! differing bonds [2, n]
end type

type :: registry_t
   logical(lk) :: overflow     ! some minimum did not fit in records
   integer(ik) :: total_steps
   integer(ik) :: n_trials
   integer(ik) :: n_records  ! occupied records
   type(record_t), dimension(:), allocatable :: records
end type

contains

subroutine allocate_registry(registry, max_records)
   type(registry_t), intent(inout) :: registry
   integer(ik), intent(in) :: max_records

   if (max_records < 1) then
      error stop 'max_records is less than 1'
   end if

   allocate (registry%records(max_records))
end subroutine

subroutine reset_registry(registry)
   type(registry_t), intent(inout) :: registry

   registry%n_records = 0
   registry%n_trials = 0
   registry%total_steps = 0
   registry%overflow = .FALSE.
   registry%records%freq = 0
   registry%records%steps = 0
   registry%records%mapdiff = huge(registry%records(1)%mapdiff)
   registry%records%mapdist = huge(registry%records(1)%mapdist)
end subroutine

subroutine insert_record_mapping(registry, mapping1, steps, rotation, mapdiff, mapdist)
! Count a trial that ended at mapping1. A known permutation only updates its
! frequency and mean steps; a new one is inserted in rank order (mapdiff,
! then mapdist), dropping the last record if the registry is full.
   type(registry_t), target, intent(inout) :: registry
   integer(ik), dimension(:), intent(in) :: mapping1
   real(rk), intent(in) :: steps, rotation(4)
   integer(ik), intent(in) :: mapdiff
   real(rk), intent(in) :: mapdist
   ! Local variables
   type(record_t), pointer :: record
   integer(ik) :: i, j, insert_pos

   registry%n_trials = registry%n_trials + 1
   registry%total_steps = registry%total_steps + steps

   ! Known permutation
   do i = 1, registry%n_records
      record => registry%records(i)
      if (all(mapping1 == record%mapping1)) then
         record%freq = record%freq + 1
         record%steps = record%steps + (steps - record%steps)/record%freq
         return
      end if
   end do

   ! Rank of the new permutation
   insert_pos = registry%n_records + 1
   do i = 1, registry%n_records
      record => registry%records(i)
      if (mapdiff < record%mapdiff) then
         insert_pos = i
         exit
      else if (mapdiff == record%mapdiff) then
         if (mapdist < record%mapdist) then
            insert_pos = i
            exit
         end if
      end if
   end do

   if (insert_pos <= size(registry%records)) then
      ! Make room (if full, the last record is dropped)
      do j = min(registry%n_records, size(registry%records) - 1), insert_pos, -1
         registry%records(j + 1) = registry%records(j)
      end do

      record => registry%records(insert_pos)
      record%mapping1 = mapping1
      record%freq = 1
      record%mapdiff = mapdiff
      record%mapdist = mapdist
      record%rotation = rotation
      record%steps = steps

      if (registry%n_records < size(registry%records)) then
         registry%n_records = registry%n_records + 1
      else
         registry%overflow = .TRUE.
      end if
   else
      registry%overflow = .TRUE.
   end if
end subroutine

subroutine insert_record_moldiff(registry, moldiffs, mapping1, steps, rotation, mapdist)
! As insert_record_mapping, but minima are grouped by the set of differing
! bonds moldiffs instead of by permutation. A group keeps the permutation
! with the lowest mapdist, and is moved to its new rank when it improves.
   type(registry_t), target, intent(inout) :: registry
   integer(ik), dimension(:,:), intent(in) :: moldiffs
   integer(ik), dimension(:), intent(in) :: mapping1
   integer(ik), intent(in) :: steps
   real(rk), intent(in) :: rotation(4)
   real(rk), intent(in) :: mapdist
   type(record_t), pointer :: record
   type(record_t) :: temp_record
   integer(ik) :: i, j, insert_pos, mapdiff, match_pos
   logical(lk) :: bonds_match, found_match

   registry%n_trials = registry%n_trials + 1
   registry%total_steps = registry%total_steps + steps

   mapdiff = size(moldiffs, 2)
   found_match = .FALSE.
   match_pos = 0

   do i = 1, registry%n_records
      record => registry%records(i)

      if (mapdiff == record%mapdiff) then
         bonds_match = all(moldiffs(1, :) == record%moldiffs(1, :)) .and. &
                       all(moldiffs(2, :) == record%moldiffs(2, :))

         if (bonds_match) then
            found_match = .TRUE.
            match_pos = i
            record%freq = record%freq + 1
            record%steps = record%steps + (steps - record%steps)/record%freq

            if (mapdist < record%mapdist) then
               temp_record = record
               temp_record%mapping1 = mapping1
               temp_record%mapdist = mapdist
               temp_record%rotation = rotation

               do j = i, registry%n_records - 1
                  registry%records(j) = registry%records(j + 1)
               end do
               registry%n_records = registry%n_records - 1

               insert_pos = registry%n_records + 1
               do j = 1, registry%n_records
                  if (mapdiff < registry%records(j)%mapdiff .or. &
                      (mapdiff == registry%records(j)%mapdiff .and. &
                       mapdist < registry%records(j)%mapdist)) then
                     insert_pos = j
                     exit
                  end if
               end do

               if (insert_pos <= size(registry%records)) then
                  do j = min(registry%n_records, size(registry%records) - 1), insert_pos, -1
                     registry%records(j + 1) = registry%records(j)
                  end do

                  registry%records(insert_pos) = temp_record

                  if (registry%n_records < size(registry%records)) then
                     registry%n_records = registry%n_records + 1
                  else
                     registry%overflow = .TRUE.
                  end if
               else
                  registry%overflow = .TRUE.
               end if
            end if
            return
         end if
      end if
   end do

   insert_pos = registry%n_records + 1
   do i = 1, registry%n_records
      if (mapdiff < registry%records(i)%mapdiff .or. &
          (mapdiff == registry%records(i)%mapdiff .and. &
           mapdist < registry%records(i)%mapdist)) then
         insert_pos = i
         exit
      end if
   end do

   if (insert_pos <= size(registry%records)) then
      do j = min(registry%n_records, size(registry%records) - 1), insert_pos, -1
         registry%records(j + 1) = registry%records(j)
      end do

      record => registry%records(insert_pos)
      record%mapping1 = mapping1
      record%freq = 1
      record%moldiffs = moldiffs
      record%mapdiff = mapdiff
      record%mapdist = mapdist
      record%rotation = rotation
      record%steps = steps

      if (registry%n_records < size(registry%records)) then
         registry%n_records = registry%n_records + 1
      else
         registry%overflow = .TRUE.
      end if
   else
      registry%overflow = .TRUE.
   end if
end subroutine

subroutine print_records(registry)
! Print the ranked minima and the search statistics
   type(registry_t), intent(in) :: registry
   type(record_t) :: record
   character(:), allocatable :: line
   integer(ik) :: i

   line = repeat('-', 39)
   write (stdout, '(2x,a,4x,a,5x,a,4x,a,5x,a,6x,a)') '#', 'Freq', 'Steps', 'Δadj', 'Δxyz'
   write (stdout, '(a)') line
   do i = 1, registry%n_records
      record = registry%records(i)
      write (stdout, '(i3,4x,i4,4x,f5.1,3x,i4,4x,f8.4)') &
         i, record%freq, record%steps, record%mapdiff, record%mapdist
   end do
   write (stdout, '(a)') line

   write (stdout, *)
   write (stdout, '(a,1x,i0)') 'Random trials:', registry%n_trials
   write (stdout, '(a,1x,i0)') 'Minimization steps:', registry%total_steps

   if (registry%overflow) then
      write (stdout, '(a,1x,i0)') 'Visited local minima: >', registry%n_records
   else
      write (stdout, '(a,1x,i0)') 'Visited local minima:', registry%n_records
   end if

   write (stdout, *)
end subroutine

end module
