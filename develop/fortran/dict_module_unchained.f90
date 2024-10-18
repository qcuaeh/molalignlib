module dict_mod
   use, intrinsic :: iso_fortran_env, only: int32, int64
   implicit none
   private

   public :: dict, create_dict

   integer, parameter :: MAX_KEY_SIZE = 16
   integer, parameter :: HASH_CONSTANT = 1779033703

   type, public :: dict
      integer, allocatable, public :: keys(:,:)
      integer, allocatable :: key_lengths(:)
      integer, allocatable :: values(:)
      logical, allocatable :: occupied(:)
      integer, public :: actual_capacity
      integer, public :: num_collision_levels
      integer, allocatable, public :: collisions(:)
      integer, allocatable, public :: occupations(:)
   contains
      procedure :: put => dict_put
      procedure :: get => dict_get
      procedure :: exists => dict_exists
      procedure :: reset => dict_reset
      final :: dict_finalize
   end type dict

contains

   function create_dict(requested_capacity) result(d)
      integer, intent(in) :: requested_capacity
      type(dict) :: d
      integer :: i, level_capacity

      d%actual_capacity = 1
      do while (d%actual_capacity < requested_capacity)
         d%actual_capacity = d%actual_capacity * 2
      end do

      d%num_collision_levels = 0
      level_capacity = d%actual_capacity / 2
      do while (level_capacity > 0)
         d%num_collision_levels = d%num_collision_levels + 1
         level_capacity = level_capacity / 2
      end do

      allocate(d%keys(MAX_KEY_SIZE, 2 * d%actual_capacity))
      allocate(d%key_lengths(2 * d%actual_capacity))
      allocate(d%values(2 * d%actual_capacity))
      allocate(d%occupied(2 * d%actual_capacity))
      allocate(d%collisions(d%num_collision_levels))
      allocate(d%occupations(0:d%num_collision_levels))

      d%occupied = .false.
      d%collisions = 0
      d%occupations = 0
   end function create_dict

   subroutine dict_put(this, key, value)
      class(dict), intent(inout) :: this
      integer, intent(in) :: key(:)
      integer, intent(in) :: value
      integer :: hash
      integer :: level_index, level, level_offset, level_capacity, bucket_index
      logical :: inserted

      hash = compute_hash(key)
      level_index = modulo(hash, this%actual_capacity)

      level = 0
      level_offset = 0
      level_capacity = this%actual_capacity
      inserted = .false.

      do while (.not. inserted)
         bucket_index = level_index + level_offset + 1

         if (this%occupied(bucket_index)) then
            if (are_keys_equivalent(this%keys(:this%key_lengths(bucket_index), bucket_index), key)) then
               this%values(bucket_index) = value
               inserted = .true.
            else
               level = level + 1
               if (level > this%num_collision_levels) then
                  error stop "Dictionary is full"
               end if
               this%collisions(level) = this%collisions(level) + 1
               level_offset = level_offset + level_capacity
               level_capacity = level_capacity / 2
               level_index = level_index / 2
            end if
         else
            this%keys(:size(key), bucket_index) = key
            this%key_lengths(bucket_index) = size(key)
            this%values(bucket_index) = value
            this%occupied(bucket_index) = .true.
            this%occupations(level) = this%occupations(level) + 1
            inserted = .true.
         end if
      end do
   end subroutine dict_put

   function dict_get(this, key) result(value)
      class(dict), intent(in) :: this
      integer, intent(in) :: key(:)
      integer :: value
      integer :: hash
      integer :: level_index, level, level_offset, level_capacity, bucket_index

      hash = compute_hash(key)
      level_index = modulo(hash, this%actual_capacity)

      level = 0
      level_offset = 0
      level_capacity = this%actual_capacity

      do while (level <= this%num_collision_levels)
         bucket_index = level_index + level_offset + 1
         if (this%occupied(bucket_index)) then
            if (are_keys_equivalent(this%keys(:this%key_lengths(bucket_index), bucket_index), key)) then
               value = this%values(bucket_index)
               return
            end if
         else
            exit
         end if
         level = level + 1
         level_offset = level_offset + level_capacity
         level_capacity = level_capacity / 2
         level_index = level_index / 2
      end do

      error stop "Key not found"
   end function dict_get

   function dict_exists(this, key) result(exists)
      class(dict), intent(in) :: this
      integer, intent(in) :: key(:)
      logical :: exists
      integer :: hash
      integer :: level_index, level, level_offset, level_capacity, bucket_index

      hash = compute_hash(key)
      level_index = modulo(hash, this%actual_capacity)

      level = 0
      level_offset = 0
      level_capacity = this%actual_capacity

      do while (level <= this%num_collision_levels)
         bucket_index = level_index + level_offset + 1
         if (this%occupied(bucket_index)) then
            if (are_keys_equivalent(this%keys(:this%key_lengths(bucket_index), bucket_index), key)) then
               exists = .true.
               return
            end if
         else
            exit
         end if
         level = level + 1
         level_offset = level_offset + level_capacity
         level_capacity = level_capacity / 2
         level_index = level_index / 2
      end do

      exists = .false.
   end function dict_exists

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