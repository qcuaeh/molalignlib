module dict_mod
   use, intrinsic :: iso_fortran_env, only: int32, int64
   implicit none
   private

   public :: dict, create_dict

   integer, parameter :: MAX_KEY_SIZE = 16
   integer, parameter :: HASH_CONSTANT = 1640531527

   type, public :: dict
      integer, allocatable, public :: keys(:,:)
      integer, allocatable :: key_lengths(:)
      integer, allocatable :: values(:)
      logical, allocatable :: occupied(:)
      integer, public :: size
      integer, public :: num_collision_levels
      integer, allocatable, public :: collisions(:)
      integer, allocatable, public :: occupations(:)
   contains
      procedure :: add => dict_add
      procedure :: get => dict_get
      procedure :: has => dict_has
      procedure :: reset => dict_reset
      final :: dict_finalize
   end type dict

contains

   function create_dict(least_size) result(d)
      integer, intent(in) :: least_size
      type(dict) :: d
      integer :: i, collision_level_size, base_size

      base_size = 1
      do while (base_size < least_size)
         base_size = base_size * 2
      end do
      
      d%size = 2 * base_size

      d%num_collision_levels = 0
      collision_level_size = base_size / 2
      do while (collision_level_size > 0)
         d%num_collision_levels = d%num_collision_levels + 1
         collision_level_size = collision_level_size / 2
      end do

      allocate(d%keys(MAX_KEY_SIZE, d%size))
      allocate(d%key_lengths(d%size))
      allocate(d%values(d%size))
      allocate(d%occupied(d%size))
      allocate(d%collisions(d%num_collision_levels))
      allocate(d%occupations(0:d%num_collision_levels))

      d%occupied = .false.
      d%collisions = 0
      d%occupations = 0
   end function create_dict

subroutine dict_add(this, key, value)
      class(dict), intent(inout) :: this
      integer, intent(in) :: key(:)
      integer, intent(in) :: value
      integer :: hash
      integer :: collision_level_index, collision_level, collision_level_offset
      integer :: collision_level_size, index, i, j, k
      integer :: first_pos, last_pos, level_size, prev_pos
      integer :: level_offset, num_gaps
      integer, allocatable :: positions(:), gaps(:)
      real :: mean_gap, std_gap, min_gap, max_gap
      real :: uniformity

      hash = compute_hash(key)
      collision_level_index = modulo(hash, this%size/2)

      collision_level = 0
      collision_level_offset = 0
      collision_level_size = this%size/2

      do
         index = collision_level_index + collision_level_offset + 1

         if (this%occupied(index)) then
            if (are_keys_equivalent(this%keys(:this%key_lengths(index), index), key)) then
               this%values(index) = value
               return
            end if

            collision_level = collision_level + 1
            if (collision_level > this%num_collision_levels) then
               print *
               print '(16(i0,1x))', key
               print '(A)', "Dictionary Full - Occupation Statistics:"
               do i = 0, this%num_collision_levels
                  first_pos = -1
                  last_pos = -1
                  level_size = this%size/2**(i+1)
                  
                  ! Calculate level offset
                  if (i == 0) then
                     level_offset = 0
                  else
                     level_offset = this%size/2
                     do j = 1, i-1
                        level_offset = level_offset + this%size/2**(j+1)
                     end do
                  end if

                  ! Collect occupied positions
                  allocate(positions(this%occupations(i)))
                  k = 0
                  do j = 1, level_size
                     index = level_offset + j
                     if (this%occupied(index)) then
                        k = k + 1
                        positions(k) = j
                        if (first_pos == -1) first_pos = j
                        last_pos = j
                     end if
                  end do

                  print '(A,I2,A,I10,A,F6.2,A)', "  Level ", i, " occupations: ", &
                     this%occupations(i), " (", &
                     real(this%occupations(i)) / real(level_size) * 100, "%)"
                  if (this%occupations(i) > 0) then
                     if (this%occupations(i) == 1) then
                        print '(A,I0)', "          Single position: ", first_pos
                     else
                        print '(A,I0,A,I0)', "          Range: ", first_pos, " to ", last_pos
                        
                        ! Calculate gaps between occupied positions
                        allocate(gaps(this%occupations(i)-1))
                        do j = 1, this%occupations(i)-1
                           gaps(j) = positions(j+1) - positions(j)
                        end do

                        ! Calculate gap statistics
                        min_gap = real(minval(gaps))
                        max_gap = real(maxval(gaps))
                        mean_gap = sum(real(gaps)) / size(gaps)
                        std_gap = sqrt(sum((real(gaps) - mean_gap)**2) / size(gaps))
                        uniformity = mean_gap / std_gap ! Higher means more uniform

                        print '(A,F5.2)', "          Uniformity: ", uniformity
                        print '(A,F5.2)', "          Min/Max gap ratio: ", min_gap/max_gap
                        deallocate(gaps)
                     end if
                  end if
                  deallocate(positions)
               end do
               error stop "Dictionary is full"
            end if

            this%collisions(collision_level) = this%collisions(collision_level) + 1
            collision_level_offset = collision_level_offset + collision_level_size
            collision_level_size = collision_level_size / 2
            collision_level_index = collision_level_index / 2
         else
            this%keys(:size(key), index) = key
            this%key_lengths(index) = size(key)
            this%values(index) = value
            this%occupied(index) = .true.
            this%occupations(collision_level) = this%occupations(collision_level) + 1
            return
         end if
      end do
   end subroutine dict_add

   function dict_get(this, key) result(value)
      class(dict), intent(in) :: this
      integer, intent(in) :: key(:)
      integer :: value
      integer :: hash
      integer :: collision_level_index, collision_level, collision_level_offset
      integer :: collision_level_size, index

      hash = compute_hash(key)
      collision_level_index = modulo(hash, this%size/2)

      collision_level = 0
      collision_level_offset = 0
      collision_level_size = this%size/2

      do
         index = collision_level_index + collision_level_offset + 1

         if (this%occupied(index)) then
            if (are_keys_equivalent(this%keys(:this%key_lengths(index), index), key)) then
               value = this%values(index)
               return
            end if
         else
            error stop "Key not found"
         end if

         collision_level = collision_level + 1
         if (collision_level > this%num_collision_levels) then
            error stop "Key not found"
         end if

         collision_level_offset = collision_level_offset + collision_level_size
         collision_level_size = collision_level_size / 2
         collision_level_index = collision_level_index / 2
      end do
   end function dict_get

   function dict_has(this, key) result(exists)
      class(dict), intent(in) :: this
      integer, intent(in) :: key(:)
      logical :: exists
      integer :: hash
      integer :: collision_level_index, collision_level, collision_level_offset
      integer :: collision_level_size, index

      hash = compute_hash(key)
      collision_level_index = modulo(hash, this%size/2)

      collision_level = 0
      collision_level_offset = 0
      collision_level_size = this%size/2

      do
         index = collision_level_index + collision_level_offset + 1

         if (this%occupied(index)) then
            if (are_keys_equivalent(this%keys(:this%key_lengths(index), index), key)) then
               exists = .true.
               return
            end if
         else
            exists = .false.
            return
         end if

         collision_level = collision_level + 1
         if (collision_level > this%num_collision_levels) then
            exists = .false.
            return
         end if

         collision_level_offset = collision_level_offset + collision_level_size
         collision_level_size = collision_level_size / 2
         collision_level_index = collision_level_index / 2
      end do
   end function dict_has

   subroutine dict_reset(this)
      class(dict), intent(inout) :: this
      this%occupied = .false.
      this%collisions = 0
      this%occupations = 0
   end subroutine dict_reset

   subroutine dict_finalize(this)
      type(dict), intent(inout) :: this
      if (allocated(this%keys)) deallocate(this%keys)
      if (allocated(this%key_lengths)) deallocate(this%key_lengths)
      if (allocated(this%values)) deallocate(this%values)
      if (allocated(this%occupied)) deallocate(this%occupied)
      if (allocated(this%collisions)) deallocate(this%collisions)
      if (allocated(this%occupations)) deallocate(this%occupations)
   end subroutine dict_finalize

   function are_keys_equivalent(key1, key2) result(equivalent)
      integer, intent(in) :: key1(:), key2(:)
      logical :: equivalent
      integer :: i

      if (size(key1) /= size(key2)) then
         equivalent = .false.
         return
      end if

      equivalent = .true.
      do i = 1, size(key1)
         if (count(key1 == key1(i)) /= count(key2 == key1(i))) then
            equivalent = .false.
            return
         end if
      end do
   end function are_keys_equivalent

   function compute_hash(key) result(hash)
      integer, intent(in) :: key(:)
      integer :: hash
      integer :: i

      hash = 1
      do i = 1, size(key)
         hash = hash * (HASH_CONSTANT + 2 * key(i))
      end do
      hash = hash / 2
   end function compute_hash

end module dict_mod