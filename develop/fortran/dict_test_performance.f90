program test_dict_part2b
   use dict_mod
   use, intrinsic :: iso_fortran_env, only: int64
   implicit none

   integer, parameter :: MAX_KEY_VALUE = 100
   integer, parameter :: MAX_KEY_LENGTH = 8
   integer, parameter :: NUM_KEYS = 1000
   integer, parameter :: PROGRESS_INTERVAL = 100000
   real, parameter :: DICT_SIZE_FACTOR = 1.0

   type(dict) :: d
   integer, allocatable :: keys(:,:)
   integer :: i, j, value, num_repetitions
   integer(int64) :: start_time, end_time, count_rate, count_max
   real :: insertion_time, retrieval_time
   integer :: total_operations, current_operation
   character(len=32) :: arg

   ! Check for command-line argument for number of repetitions
   if (command_argument_count() > 0) then
      call get_command_argument(1, arg)
      read(arg, *) num_repetitions
   else
      num_repetitions = 1
   end if

   print '(A,I0)', "Number of keys to be generated: ", NUM_KEYS

   ! Create dictionary
   d = create_dict(int(NUM_KEYS * DICT_SIZE_FACTOR))
   print '(A,I0)', "Dictionary size: ", d%size

   ! Allocate keys array
   allocate(keys(MAX_KEY_LENGTH, NUM_KEYS))

   total_operations = 2 * num_repetitions * NUM_KEYS
   current_operation = 0

   insertion_time = 0.0
   retrieval_time = 0.0

   ! Print initial progress
   write(*, '(A)', advance='no') "Progress: 0%"
   flush(6)

   do j = 1, num_repetitions
      call d%reset()  ! Reset the dictionary for each repetition

      ! Generate all keys for this repetition
      do i = 1, NUM_KEYS
         call generate_random_key(keys(:, i), MAX_KEY_LENGTH, MAX_KEY_VALUE)
      end do

      ! Addition loop
      call system_clock(start_time, count_rate, count_max)
      do i = 1, NUM_KEYS
         value = i
         call d%add(keys(:keys(1,i), i), value)

         current_operation = current_operation + 1
         if (mod(current_operation, PROGRESS_INTERVAL) == 0) then
            write(*, '(A,I3,A)', advance='no') char(13)//"Progress: ", &
               int(real(current_operation) / real(total_operations) * 100), "%"
            flush(6)
         end if
      end do
      call system_clock(end_time)
      insertion_time = insertion_time + real(end_time - start_time) / real(count_rate)

      ! Retrieval loop
      call system_clock(start_time, count_rate, count_max)
      do i = 1, NUM_KEYS
         value = d%get(keys(:keys(1,i), i))

         current_operation = current_operation + 1
         if (mod(current_operation, PROGRESS_INTERVAL) == 0) then
            write(*, '(A,I3,A)', advance='no') char(13)//"Progress: ", &
               int(real(current_operation) / real(total_operations) * 100), "%"
            flush(6)
         end if
      end do
      call system_clock(end_time)
      retrieval_time = retrieval_time + real(end_time - start_time) / real(count_rate)
   end do

   ! Print final progress
   write(*, '(A)') char(13)//"Progress: 100%"

   deallocate(keys)

   print '(A,I0)', "Repetitions: ", num_repetitions
   print '(A,F10.3,A)', "Total insertion time: ", insertion_time, " seconds"
   print '(A,F10.3,A)', "Total retrieval time: ", retrieval_time, " seconds"

contains

   subroutine generate_random_key(key, max_length, max_value)
      integer, intent(out) :: key(:)
      integer, intent(in) :: max_length, max_value
      integer :: length, i
      real :: r

      call random_number(r)
      length = int(r * max_length) + 1
      key(1) = length  ! Store the length in the first element

      do i = 2, length + 1
         call random_number(r)
         key(i) = int(r * max_value) + 1
      end do

      ! Fill the rest of the array with zeros (if any)
      key(length+2:) = 0
   end subroutine generate_random_key

end program test_dict_part2b