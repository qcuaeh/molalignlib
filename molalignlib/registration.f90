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
use adjacency
use spatial_transforms
implicit none
private

public record_t
public registry_t
public push_record
public print_records
public init_rmsd_registry
public init_adjd_registry
public init_dual_registry

type :: record_t
   integer :: count
   integer :: adjd                                      ! Used when use_adjacency=TRUE
   real(rk) :: rmsd                                     ! Used when use_position=TRUE
   real(rk) :: aver_steps                               ! Used when use_position=TRUE
   real(rk) :: aver_rotangle                            ! Used when use_position=TRUE
   integer, dimension(:), allocatable :: atomperm
   real(rk), dimension(:,:), allocatable :: coords2     ! Used when use_position=TRUE
   logical, dimension(:,:), allocatable :: adjmat2      ! Used when use_adjacency=TRUE
end type

type :: registry_t
   logical :: use_position                              ! Enable position-based processing (RMSD, steps, rotation)
   logical :: use_adjacency                             ! Enable adjacency-based processing (adjacency differences)
   logical :: overflow
   integer :: total_steps                               ! Used when use_position=TRUE
   integer :: num_trials
   integer :: num_records
   real(rk), dimension(:,:), allocatable :: coords1     ! Used when use_position=TRUE
   logical, dimension(:,:), allocatable :: adjmat1      ! Used when use_adjacency=TRUE
   type(record_t), dimension(:), allocatable :: records
end type

contains

subroutine init_rmsd_registry(self, max_records, coords1)
   class(registry_t), intent(inout) :: self
   integer, intent(in) :: max_records
   real(rk), dimension(:,:), intent(in) :: coords1

   if (max_records < 1) then
      error stop 'max_records < 1'
   end if

   self%use_position = .true.
   self%use_adjacency = .false.
   self%num_records = 0
   self%num_trials = 0
   self%total_steps = 0
   self%overflow = .false.
   self%coords1 = coords1

   allocate (self%records(max_records))

   self%records%count = 0
   self%records%rmsd = huge(self%records(1)%rmsd)
end subroutine

subroutine init_dual_registry(self, max_records, coords1, adjmat1)
   class(registry_t), intent(inout) :: self
   integer, intent(in) :: max_records
   real(rk), dimension(:,:), intent(in) :: coords1
   logical, dimension(:,:), intent(in) :: adjmat1

   if (max_records < 1) then
      error stop 'max_records < 1'
   end if

   self%use_position = .true.
   self%use_adjacency = .true.
   self%num_records = 0
   self%num_trials = 0
   self%total_steps = 0
   self%overflow = .false.
   self%coords1 = coords1
   self%adjmat1 = adjmat1

   allocate (self%records(max_records))

   self%records%count = 0
   self%records%adjd = huge(self%records(1)%adjd)
end subroutine

subroutine init_adjd_registry(self, max_records, adjmat1)
   class(registry_t), intent(inout) :: self
   integer, intent(in) :: max_records
   logical, dimension(:,:), intent(in) :: adjmat1

   if (max_records < 1) then
      error stop 'max_records < 1'
   end if

   self%use_position = .false.
   self%use_adjacency = .true.
   self%num_records = 0
   self%num_trials = 0
   self%overflow = .false.
   self%adjmat1 = adjmat1

   allocate (self%records(max_records))

   self%records%count = 0
   self%records%adjd = huge(self%records(1)%adjd)
end subroutine

subroutine push_record(self, atomperm, coords2, adjmat2, num_steps, rotation)
   class(registry_t), target, intent(inout) :: self
   integer, dimension(:), intent(in) :: atomperm
   real(rk), dimension(:,:), intent(in), optional :: coords2
   logical, dimension(:,:), intent(in), optional :: adjmat2
   integer, intent(in), optional :: num_steps
   real(rk), intent(in), optional :: rotation(4)
   ! Local variables
   type(record_t), pointer :: record
   real(rk) :: rmsd
   integer :: adjd, steps
   real(rk) :: rot_angle
   integer :: i, j
   logical :: should_insert

   ! Validate required parameters based on flags
   if (.not. self%use_position .and. .not. self%use_adjacency) then
      error stop 'Registry not properly initialized - no processing mode enabled'
   end if

   self%num_trials = self%num_trials + 1
   
   ! Handle position-based processing
   if (self%use_position) then
      if (.not. present(num_steps)) error stop 'num_steps required when use_position=TRUE'
      if (.not. present(rotation)) error stop 'rotation required when use_position=TRUE'
      if (.not. present(coords2)) error stop 'coords2 required when use_position=TRUE'
      steps = num_steps
      rot_angle = angle(rotation)
      self%total_steps = self%total_steps + steps
   end if
   
   ! Handle adjacency-based processing
   if (self%use_adjacency) then
      if (.not. present(adjmat2)) error stop 'adjmat2 required when use_adjacency=TRUE'
      adjd = adjacencydiff(atomperm, self%adjmat1, adjmat2)
   end if

   ! Check for existing records to update (CRITICAL: topoatomperm has different logic!)
   do i = 1, self%num_records
      record => self%records(i)
      if (allocated(record%atomperm)) then
         if ((self%use_adjacency .and. .not. self%use_position .and. adjd == record%adjd) .or. &
             (.not. (self%use_adjacency .and. .not. self%use_position) .and. all(atomperm == record%atomperm))) then
            record%count = record%count + 1
            if (self%use_position) then
               record%aver_steps = record%aver_steps + (steps - record%aver_steps) / record%count
               record%aver_rotangle = record%aver_rotangle + (rot_angle - record%aver_rotangle) / record%count
            end if
            return
         end if
      end if
   end do

   ! Calculate RMSD if needed (after the early return check)
   if (self%use_position) then
      rmsd = sqrt(total_sqdist(atomperm, self%coords1, coords2))
   end if

   ! Find insertion point and insert new record
   do i = 1, size(self%records)
      record => self%records(i)
      
      ! Determine insertion criteria based on mode
      if (self%use_adjacency .and. self%use_position) then
         should_insert = (adjd < record%adjd .or. (adjd == record%adjd .and. rmsd < record%rmsd))
      else if (self%use_position) then
         should_insert = (rmsd < record%rmsd)
      else
         should_insert = (adjd < record%adjd)
      end if
      
      if (should_insert) then
         ! Shift records to make room
         do j = size(self%records), i + 1, -1
            self%records(j) = self%records(j - 1)
         end do
         
         ! Initialize new record
         record%atomperm = atomperm
         record%count = 1
         
         ! Set fields based on enabled processing modes
         if (self%use_position) then
            record%rmsd = rmsd
            record%coords2 = coords2
            record%aver_steps = steps
            record%aver_rotangle = rot_angle
         end if
         
         if (self%use_adjacency) then
            record%adjd = adjd
            record%adjmat2 = adjmat2
         end if
         
         exit
      end if
   end do

   ! Update record count and overflow status
   if (.not. self%overflow) then
      if (self%num_records < size(self%records)) then
         self%num_records = self%num_records + 1
      else
         self%overflow = .true.
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
      write (stdout, '(2x,a,4x,a,4x,a,4x,a,7x,a)') '#', 'Count', 'Steps', 'Rot-θ', 'RMSD'
      write (stdout, '(a)') line(1:42)
      do i = 1, registry%num_records
         record = registry%records(i)
         write (stdout, '(i3,4x,i4,4x,f5.1,5x,f5.1,4x,f8.4)') &
            i, record%count, record%aver_steps, record%aver_rotangle, record%rmsd
      end do
      write (stdout, '(a)') line(1:42)
   
   else if (registry%use_position .and. registry%use_adjacency) then
      ! Both position and adjacency
      line = repeat('-', 49)
      write (stdout, '(2x,a,4x,a,4x,a,4x,a,4x,a,6x,a)') '#', 'Count', 'Steps', 'Rot-θ', 'Δadj', 'RMSD'
      write (stdout, '(a)') line
      do i = 1, registry%num_records
         record = registry%records(i)
         write (stdout, '(i3,4x,i4,4x,f5.1,5x,f5.1,3x,i4,4x,f8.4)') &
            i, record%count, record%aver_steps, record%aver_rotangle, record%adjd, record%rmsd
      end do
      write (stdout, '(a)') line
   
   else if (registry%use_adjacency .and. .not. registry%use_position) then
      ! Adjacency only
      line = repeat('-', 25)
      write (stdout, '(2x,a,4x,a,4x,a)') '#', 'Count', 'Δadj'
      write (stdout, '(a)') line(1:25)
      do i = 1, registry%num_records
         record = registry%records(i)
         write (stdout, '(i3,4x,i4,4x,i4)') i, record%count, record%adjd
      end do
      write (stdout, '(a)') line(1:25)
   end if

   write (stdout, '(a,1x,i0)') 'Random trials =', registry%num_trials
   
   if (registry%use_position) then
      write (stdout, '(a,1x,i0)') 'Minimization steps =', registry%total_steps
   end if
   
   if (registry%overflow) then
      write (stdout, '(a,1x,i0)') 'Visited local minima >', registry%num_records
   else
      write (stdout, '(a,1x,i0)') 'Visited local minima =', registry%num_records
   end if
   flush(stdout)
end subroutine

end module
