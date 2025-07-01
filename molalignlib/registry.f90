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

module registry
use parameters
use spatial_transforms
use adjacency

implicit none

type :: atomperm_record
   integer :: count
   real(rk) :: rmsd
   real(rk) :: aver_steps
   real(rk) :: aver_rotangle
   integer, dimension(:), allocatable :: atomperm
   real(rk), dimension(:,:), allocatable :: coords2
end type

type :: atomperm_registry
   logical :: overflow
   integer :: total_steps
   integer :: num_trials
   integer :: num_records
   real(rk), dimension(:,:), allocatable :: coords1
   type(atomperm_record), dimension(:), allocatable :: records
end type

type :: bondatomperm_record
   integer :: count
   integer :: adjd
   real(rk) :: rmsd
   real(rk) :: aver_steps
   real(rk) :: aver_rotangle
   integer, dimension(:), allocatable :: atomperm
   real(rk), dimension(:,:), allocatable :: coords2
   logical, dimension(:,:), allocatable :: adjmat2
end type

type :: bondatomperm_registry
   logical :: overflow
   integer :: total_steps
   integer :: num_trials
   integer :: num_records
   real(rk), dimension(:,:), allocatable :: coords1
   logical, dimension(:,:), allocatable :: adjmat1
   type(bondatomperm_record), dimension(:), allocatable :: records
end type

type :: topoatomperm_record
   integer :: count
   integer :: adjd
   integer, dimension(:), allocatable :: atomperm
   logical, dimension(:,:), allocatable :: adjmat2
end type

type :: topoatomperm_registry
   logical :: overflow
   integer :: num_trials
   integer :: num_records
   logical, dimension(:,:), allocatable :: adjmat1
   type(topoatomperm_record), dimension(:), allocatable :: records
end type

interface registry_init
   module procedure atomperm_registry_init
   module procedure bondatomperm_registry_init
   module procedure topoatomperm_registry_init
end interface

interface registry_push
   module procedure atomperm_registry_push
   module procedure bondatomperm_registry_push
   module procedure topoatomperm_registry_push
end interface

interface registry_print
   module procedure atomperm_registry_print
   module procedure bondatomperm_registry_print
end interface

contains

subroutine atomperm_registry_init(self, max_records, coords1)
   class(atomperm_registry), intent(inout) :: self
   integer, intent(in) :: max_records
   real(rk), dimension(:,:), intent(in) :: coords1

   if (max_records < 1) then
      error stop 'max_records < 1'
   end if

   self%num_records = 0
   self%num_trials = 0
   self%total_steps = 0
   self%overflow = .false.
   self%coords1 = coords1

   allocate (self%records(max_records))

   self%records%count = 0
   self%records%rmsd = huge(self%records(1)%rmsd)
end subroutine

subroutine atomperm_registry_push(self, coords2, atomperm, num_steps, rotation)
   class(atomperm_registry), target, intent(inout) :: self
   real(rk), dimension(:,:), intent(in) :: coords2
   integer, dimension(:), intent(in) :: atomperm
   real(rk), intent(in) :: rotation(4)
   integer, intent(in) :: num_steps
   ! Local variables
   type(atomperm_record), pointer :: record
   real(rk) :: rmsd
   integer :: i, j

   self%num_trials = self%num_trials + 1
   self%total_steps = self%total_steps + num_steps

   do i = 1, self%num_records
      record => self%records(i)
      if (allocated(record%atomperm)) then
         if (all(atomperm == record%atomperm)) then
            record%count = record%count + 1
            record%aver_steps = record%aver_steps + (num_steps - record%aver_steps) / record%count
            record%aver_rotangle = record%aver_rotangle + (angle(rotation) - record%aver_rotangle) / record%count
            return
         end if
      end if
   end do

   rmsd = sqrt(total_sqdist(atomperm, self%coords1, coords2))

   do i = 1, size(self%records)
      record => self%records(i)
      if (rmsd < record%rmsd) then
         do j = size(self%records), i + 1, -1
            self%records(j) = self%records(j - 1)
         end do
         record%atomperm = atomperm
         record%count = 1
         record%rmsd = rmsd
         record%coords2 = coords2
         record%aver_steps = num_steps
         record%aver_rotangle = angle(rotation)
         exit
      end if
   end do

   if (.not. self%overflow) then
      if (self%num_records < size(self%records)) then
         self%num_records = self%num_records + 1
      else
         self%overflow = .true.
      end if
   end if
end subroutine

subroutine bondatomperm_registry_init(self, max_records, coords1, adjmat1)
   class(bondatomperm_registry), intent(inout) :: self
   integer, intent(in) :: max_records
   real(rk), dimension(:,:), intent(in) :: coords1
   logical, dimension(:,:), intent(in) :: adjmat1

   if (max_records < 1) then
      error stop 'max_records < 1'
   end if

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

subroutine bondatomperm_registry_push(self, coords2, adjmat2, atomperm, num_steps, rotation)
   class(bondatomperm_registry), target, intent(inout) :: self
   real(rk), dimension(:,:), intent(in) :: coords2
   logical, dimension(:,:), intent(in) :: adjmat2
   integer, dimension(:), intent(in) :: atomperm
   real(rk), intent(in) :: rotation(4)
   integer, intent(in) :: num_steps
   ! Local variables
   type(bondatomperm_record), pointer :: record
   real(rk) :: rmsd
   integer :: adjd
   integer :: i, j

   self%num_trials = self%num_trials + 1
   self%total_steps = self%total_steps + num_steps

   do i = 1, self%num_records
      record => self%records(i)
      if (allocated(record%atomperm)) then
         if (all(atomperm == record%atomperm)) then
            record%count = record%count + 1
            record%aver_steps = record%aver_steps + (num_steps - record%aver_steps) / record%count
            record%aver_rotangle = record%aver_rotangle + (angle(rotation) - record%aver_rotangle) / record%count
            return
         end if
      end if
   end do

   adjd = adjacencydiff(atomperm, self%adjmat1, adjmat2)
   rmsd = sqrt(total_sqdist(atomperm, self%coords1, coords2))

   do i = 1, size(self%records)
      record => self%records(i)
      if (adjd < record%adjd .or. (adjd == record%adjd .and. rmsd < record%rmsd)) then
         do j = size(self%records), i + 1, -1
            self%records(j) = self%records(j - 1)
         end do
         record%atomperm = atomperm
         record%count = 1
         record%rmsd = rmsd
         record%adjd = adjd
         record%coords2 = coords2
         record%adjmat2 = adjmat2
         record%aver_steps = num_steps
         record%aver_rotangle = angle(rotation)
         exit
      end if
   end do

   if (.not. self%overflow) then
      if (self%num_records < size(self%records)) then
         self%num_records = self%num_records + 1
      else
         self%overflow = .true.
      end if
   end if
end subroutine

subroutine topoatomperm_registry_init(self, max_records, adjmat1)
   class(topoatomperm_registry), intent(inout) :: self
   integer, intent(in) :: max_records
   logical, dimension(:,:), intent(in) :: adjmat1

   if (max_records < 1) then
      error stop 'max_records < 1'
   end if

   self%num_records = 0
   self%num_trials = 0
   self%overflow = .false.
   self%adjmat1 = adjmat1

   allocate (self%records(max_records))

   self%records%count = 0
   self%records%adjd = huge(self%records(1)%adjd)
end subroutine

subroutine topoatomperm_registry_push(self, adjmat2, atomperm)
   class(topoatomperm_registry), target, intent(inout) :: self
   logical, dimension(:,:), intent(in) :: adjmat2
   integer, dimension(:), intent(in) :: atomperm
   ! Local variables
   type(topoatomperm_record), pointer :: record
   integer :: adjd
   integer :: i, j

   self%num_trials = self%num_trials + 1

   adjd = adjacencydiff(atomperm, self%adjmat1, adjmat2)

   do i = 1, self%num_records
      record => self%records(i)
      if (allocated(record%atomperm)) then
         if (adjd == record%adjd) then
            record%count = record%count + 1
            return
         end if
      end if
   end do

   do i = 1, size(self%records)
      record => self%records(i)
      if (adjd < record%adjd) then
         do j = size(self%records), i + 1, -1
            self%records(j) = self%records(j - 1)
         end do
         record%atomperm = atomperm
         record%count = 1
         record%adjd = adjd
         record%adjmat2 = adjmat2
         exit
      end if
   end do

   if (.not. self%overflow) then
      if (self%num_records < size(self%records)) then
         self%num_records = self%num_records + 1
      else
         self%overflow = .true.
      end if
   end if
end subroutine

subroutine atomperm_registry_print(results)
   type(atomperm_registry), intent(in) :: results
   ! Parameters
   character(*), parameter :: line = repeat('-', 42)
   ! Local variables
   type(atomperm_record) :: record
   integer :: i

   write (stdout, '(2x,a,4x,a,4x,a,4x,a,7x,a)') '#', 'Count', 'Steps', 'Rot-θ', 'RMSD'
   write (stdout, '(a)') line
   do i = 1, results%num_records
      record = results%records(i)
      write (stdout, '(i3,4x,i4,4x,f5.1,5x,f5.1,4x,f8.4)') &
         i, record%count, record%aver_steps, record%aver_rotangle, record%rmsd
   end do
   write (stdout, '(a)') line
   write (stdout, '(a,1x,i0)') 'Random trials =', results%num_trials
   write (stdout, '(a,1x,i0)') 'Minimization steps =', results%total_steps
   if (results%overflow) then
      write (stdout, '(a,1x,i0)') 'Visited local minima >', results%num_records
   else
      write (stdout, '(a,1x,i0)') 'Visited local minima =', results%num_records
   end if
   flush(stdout)
end subroutine

subroutine bondatomperm_registry_print(results)
   type(bondatomperm_registry), intent(in) :: results
   ! Parameters
   character(*), parameter :: line = repeat('-', 49)
   ! Local variables
   type(bondatomperm_record) :: record
   integer :: i

   write (stdout, '(2x,a,4x,a,4x,a,4x,a,4x,a,6x,a)') '#', 'Count', 'Steps', 'Rot-θ', 'Δadj', 'RMSD'
   write (stdout, '(a)') line
   do i = 1, results%num_records
      record = results%records(i)
      write (stdout, '(i3,4x,i4,4x,f5.1,5x,f5.1,3x,i4,4x,f8.4)') &
         i, record%count, record%aver_steps, record%aver_rotangle, record%adjd
   end do
   write (stdout, '(a)') line
   write (stdout, '(a,1x,i0)') 'Random trials =', results%num_trials
   write (stdout, '(a,1x,i0)') 'Minimization steps =', results%total_steps
   if (results%overflow) then
      write (stdout, '(a,1x,i0)') 'Visited local minima >', results%num_records
   else
      write (stdout, '(a,1x,i0)') 'Visited local minima =', results%num_records
   end if
   flush(stdout)
end subroutine

end module
