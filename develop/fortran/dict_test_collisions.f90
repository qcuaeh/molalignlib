program test_dict_part2a
   use dict_mod
   implicit none

   integer, parameter :: MAX_KEY_VALUE = 30
   integer, parameter :: TEST_KEY_LENGTH = 8
   integer, parameter :: PROGRESS_INTERVAL = 100000
   real, parameter :: DICT_SIZE_FACTOR = 1.0

   type(dict) :: d
   integer :: i, test_count, success_count
   integer :: key(TEST_KEY_LENGTH)
   integer :: value
   integer :: expected_dict_size
   logical :: has_next

   ! Calculate the number of combinations
   test_count = binomial(MAX_KEY_VALUE + TEST_KEY_LENGTH - 1, TEST_KEY_LENGTH)
   print '(A,I0)', "Total combinations to be generated: ", test_count

   ! Create dictionary
   d = create_dict(int(test_count * DICT_SIZE_FACTOR))
   print '(A,I0)', "Dictionary size: ", d%size

   success_count = 0
   expected_dict_size = 0

   ! Print initial progress
   write(*, '(A)', advance='no') "Progress: 0%"
   flush(6)

   ! Initialize the first combination
   key = 1
   i = 0

   ! Generate and add all combinations
   do while (.true.)
      i = i + 1
      value = i  ! Use a simple value based on the combination index

      ! Uncomment the following line to print generated keys
      ! print *, "Generated key:", key

      call d%add(key, value)
      expected_dict_size = expected_dict_size + 1

      if (d%has(key) .and. d%get(key) == value) then
         success_count = success_count + 1
      end if

      ! Report progress
      if (mod(i, PROGRESS_INTERVAL) == 0) then
         write(*, '(A,I3,A)', advance='no') char(13)//"Progress: ", &
            int(real(i) / real(test_count) * 100), "%"
         flush(6)
      end if

      ! Generate next combination
      has_next = next_combination(key, MAX_KEY_VALUE)
      if (.not. has_next) exit
   end do

   ! Print final progress
   write(*, '(A)') char(13)//"Progress: 100%"

   ! Print results
   print '(A)', "Insertion Test Summary:"
   print '(A,I0)', "Total combinations added: ", test_count
   print '(A,I0)', "Successful insertions: ", success_count
   print '(A,I0)', "Final dictionary size: ", sum(d%occupations)

   ! Print collision statistics
   print '(A)', "Collision Statistics:"
   do i = 1, size(d%collisions)
      print '(A,I2,A,I10,A,F6.2,A)', "  Level ", i, " collisions: ", d%collisions(i), &
         " (", real(d%collisions(i)) / real(test_count) * 100, "%)"
   end do

   ! Print occupation statistics
   print '(A)', "Occupation Statistics:"
   do i = 0, d%num_collision_levels
      print '(A,I2,A,I10,A,F6.2,A)', "  Level ", i, " occupations: ", d%occupations(i), &
         " (", real(d%occupations(i)) / real(d%size/2**(i+1)) * 100, "%)"
   end do

   ! Check if all insertions were successful
   if (success_count == test_count) then
      print '(A)', "All insertions successful!"
   else
      print '(A)', "Warning: Some insertions failed or were not correctly retrieved"
      print '(A,F6.2,A)', "Success rate: ", real(success_count) / real(test_count) * 100, "%"
   end if

contains

   function binomial(n, k) result(res)
      integer, intent(in) :: n, k
      integer :: res, i
      real :: temp

      temp = 1.0
      do i = 1, k
         temp = temp * (n - i + 1) / i
      end do
      res = nint(temp)
   end function binomial

   function next_combination(a, n) result(has_next)
      integer, intent(inout) :: a(:)
      integer, intent(in) :: n
      logical :: has_next
      integer :: i, j

      j = size(a)
      
      ! Find the rightmost element that can be incremented
      do i = j, 1, -1
         if (a(i) < n) then
            a(i) = a(i) + 1
            ! Reset all elements to the right
            a(i+1:j) = a(i)
            has_next = .true.
            return
         end if
      end do
      
      has_next = .false.
   end function next_combination

end program test_dict_part2a