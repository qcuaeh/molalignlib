module dict_mod
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private

    public :: dict, create_dict

    type dict_item
        integer, allocatable :: key(:)
        integer :: value
        type(dict_item), pointer :: next => null()
    end type dict_item

    type :: dict
        private
        type(dict_item), pointer :: buckets(:) => null()
        integer :: count = 0
        integer :: capacity
        integer :: capacity_mask
    contains
        procedure, public :: insert => dict_insert
        procedure, public :: retrieve => dict_retrieve
        procedure, public :: exists => dict_exists
        procedure, public :: get_count => dict_get_count
        procedure, public :: reset => dict_reset
        final :: dict_finalize
    end type dict

contains

    function create_dict(capacity) result(d)
        integer, intent(in) :: capacity
        type(dict) :: d
        integer :: actual_capacity
        
        actual_capacity = 1
        do while (actual_capacity < capacity)
            actual_capacity = ishft(actual_capacity, 1)
        end do
        
        allocate(d%buckets(actual_capacity))
        d%capacity = actual_capacity
        d%capacity_mask = actual_capacity - 1
    end function create_dict

    subroutine dict_insert(this, key, value)
        class(dict), intent(inout) :: this
        integer, intent(in) :: key(:)
        integer, intent(in) :: value
        integer(int64) :: hash
        integer :: bucket_index
        type(dict_item), pointer :: new_item, current
        logical :: collision_occurred

        hash = compute_hash(key)
        bucket_index = iand(int(hash, kind=int64), int(this%capacity_mask, kind=int64)) + 1

        collision_occurred = .false.
        current => this%buckets(bucket_index)%next
        do while (associated(current))
            if (are_keys_equivalent(current%key, key)) then
                current%value = value
                return  ! Not a collision, just updating an existing key
            end if
            collision_occurred = .true.  ! We've found a different key in this bucket
            current => current%next
        end do

        allocate(new_item)
        allocate(new_item%key(size(key)))
        new_item%key = key
        new_item%value = value
        new_item%next => this%buckets(bucket_index)%next
        this%buckets(bucket_index)%next => new_item

        this%count = this%count + 1

        if (collision_occurred) then
            print *, "Warning: Hash collision occurred at bucket", bucket_index
        end if
    end subroutine dict_insert

    function dict_retrieve(this, key) result(value)
        class(dict), intent(in) :: this
        integer, intent(in) :: key(:)
        integer :: value
        integer(int64) :: hash
        integer :: bucket_index
        type(dict_item), pointer :: current

        hash = compute_hash(key)
        bucket_index = iand(int(hash, kind=int64), int(this%capacity_mask, kind=int64)) + 1

        current => this%buckets(bucket_index)%next
        do while (associated(current))
            if (are_keys_equivalent(current%key, key)) then
                value = current%value
                return
            end if
            current => current%next
        end do

        error stop "Key not found"
    end function dict_retrieve

    function dict_exists(this, key) result(exists)
        class(dict), intent(in) :: this
        integer, intent(in) :: key(:)
        logical :: exists
        integer(int64) :: hash
        integer :: bucket_index
        type(dict_item), pointer :: current

        hash = compute_hash(key)
        bucket_index = iand(int(hash, kind=int64), int(this%capacity_mask, kind=int64)) + 1

        current => this%buckets(bucket_index)%next
        do while (associated(current))
            if (are_keys_equivalent(current%key, key)) then
                exists = .true.
                return
            end if
            current => current%next
        end do

        exists = .false.
    end function dict_exists

    function dict_get_count(this) result(count)
        class(dict), intent(in) :: this
        integer :: count
        count = this%count
    end function dict_get_count

    subroutine dict_reset(this)
        class(dict), intent(inout) :: this
        integer :: i
        type(dict_item), pointer :: current, next

        do i = 1, this%capacity
            current => this%buckets(i)%next
            do while (associated(current))
                next => current%next
                deallocate(current%key)
                deallocate(current)
                current => next
            end do
            nullify(this%buckets(i)%next)
        end do

        this%count = 0
    end subroutine dict_reset

    subroutine dict_finalize(this)
        type(dict), intent(inout) :: this
        call this%reset()
        deallocate(this%buckets)
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
        integer(int64) :: hash
        integer :: i

        hash = 1_int64
        do i = 1, size(key)
            hash = hash * (2654435769 + 2*key(i))
        end do
        hash = modulo(hash/2, 2_int64**63 - 1_int64)  ! Keep within positive 63-bit range
    end function compute_hash

end module dict_mod

program test_dict
    use dict_mod
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none

    integer, parameter :: MAX_DICT_SIZE = 2000000
    integer, parameter :: NUM_TEST_PAIRS = 50
    integer, parameter :: NUM_TEST_CASES = 25
    integer, parameter :: MAX_KEY_LENGTH = 16
    integer, parameter :: MAX_KEY_VALUE = 100
    integer, parameter :: LONG_TEST_KEY_LENGTH = 3  ! New constant for long test key length

    type test_pair
        integer, allocatable :: key(:)
        integer :: value
    end type test_pair

    type(dict) :: d
    type(test_pair) :: first_set(NUM_TEST_PAIRS), second_set(NUM_TEST_PAIRS)
    type(test_pair) :: shuffled_first_set(NUM_TEST_PAIRS)
    integer :: i, j, test_count, success_count, total_success, total_tests
    real :: random_value
    integer(int64) :: start_time, end_time, count_rate, count_max
        integer :: key(LONG_TEST_KEY_LENGTH)
        integer :: temp, value


    ! Initialize random number generator
    call random_seed()

    ! Generate first set of key-value pairs
    call generate_unique_set(first_set, NUM_TEST_PAIRS)

    ! Generate second set of key-value pairs
    call generate_unique_set(second_set, NUM_TEST_PAIRS)

    ! Ensure second set is distinct from first set
    do i = 1, NUM_TEST_PAIRS
        do while (is_equivalent_to_any(first_set, NUM_TEST_PAIRS, second_set(i)%key))
            deallocate(second_set(i)%key)
            call generate_key_value_pair(second_set(i))
        end do
    end do

    ! Create shuffled version of first set
    do i = 1, NUM_TEST_PAIRS
        allocate(shuffled_first_set(i)%key(size(first_set(i)%key)))
        shuffled_first_set(i)%key = first_set(i)%key
        call shuffle_array(shuffled_first_set(i)%key)
        shuffled_first_set(i)%value = first_set(i)%value
    end do

    ! Create dictionary
    d = create_dict(MAX_DICT_SIZE)

    total_success = 0
    total_tests = 0

    ! Test 1: Insertion and existence of keys (first set)
    test_count = 0
    success_count = 0
    do i = 1, min(NUM_TEST_CASES, NUM_TEST_PAIRS)
        test_count = test_count + 1
        call d%insert(shuffled_first_set(i)%key, shuffled_first_set(i)%value)
        if (d%exists(first_set(i)%key)) then
            success_count = success_count + 1
        end if
    end do
    print *, "Test 1 Summary: ", success_count, "/", test_count, " tests passed"
    total_success = total_success + success_count
    total_tests = total_tests + test_count

    ! Test 2: Non-existence of keys (second set)
    test_count = 0
    success_count = 0
    do i = 1, min(NUM_TEST_CASES, NUM_TEST_PAIRS)
        test_count = test_count + 1
        if (.not. d%exists(second_set(i)%key)) then
            success_count = success_count + 1
        end if
    end do
    print *, "Test 2 Summary: ", success_count, "/", test_count, " tests passed"
    total_success = total_success + success_count
    total_tests = total_tests + test_count

    ! Test 3: Retrieval of values
    test_count = 0
    success_count = 0
    do i = 1, min(NUM_TEST_CASES, NUM_TEST_PAIRS)
        test_count = test_count + 1
        if (d%retrieve(first_set(i)%key) == first_set(i)%value) then
            success_count = success_count + 1
        end if
    end do
    print *, "Test 3 Summary: ", success_count, "/", test_count, " tests passed"
    total_success = total_success + success_count
    total_tests = total_tests + test_count

    ! Test 4: Dictionary size
    test_count = 1
    if (d%get_count() == NUM_TEST_CASES) then
        success_count = 1
    else
        success_count = 0
    end if
    print *, "Test 4 Summary: ", success_count, "/", test_count, " tests passed"
    total_success = total_success + success_count
    total_tests = total_tests + test_count

    ! Test 5: Overwriting values
    test_count = 0
    success_count = 0
    do i = 1, min(NUM_TEST_CASES, NUM_TEST_PAIRS)
        test_count = test_count + 1
        call d%insert(first_set(i)%key, first_set(i)%value * 2)
        if (d%retrieve(first_set(i)%key) == first_set(i)%value * 2) then
            success_count = success_count + 1
        end if
    end do
    print *, "Test 5 Summary: ", success_count, "/", test_count, " tests passed"
    total_success = total_success + success_count
    total_tests = total_tests + test_count

    ! Test 6: Reset function
    call d%reset()
    test_count = 2
    success_count = 0
    if (d%get_count() == 0) then
        success_count = success_count + 1
    end if
    if (.not. d%exists(first_set(1)%key)) then
        success_count = success_count + 1
    end if
    print *, "Test 6 Summary: ", success_count, "/", test_count, " tests passed"
    total_success = total_success + success_count
    total_tests = total_tests + test_count

! Test 7: Longer insertion test with all possible unique keys
    print *, "Test 7: Longer insertion test with all possible unique keys"
    call d%reset()
    test_count = MAX_KEY_VALUE**LONG_TEST_KEY_LENGTH
    success_count = 0

    call system_clock(start_time, count_rate, count_max)

    do i = 0, test_count - 1
        temp = i
        do j = LONG_TEST_KEY_LENGTH, 1, -1
            key(j) = mod(temp, MAX_KEY_VALUE) + 1
            temp = temp / MAX_KEY_VALUE
        end do

        value = i + 1  ! Use a simple value based on the key's index

        call d%insert(key, value)
        if (d%exists(key) .and. d%retrieve(key) == value) then
            success_count = success_count + 1
        end if
    end do

    call system_clock(end_time)

    print *, "Test 7 Summary: ", success_count, "/", test_count, " tests passed"
    print *, "Time taken for long insertion test: ", real(end_time - start_time) / real(count_rate), " seconds"
    print *, "Final dictionary size:", d%get_count()

    total_success = total_success + success_count
    total_tests = total_tests + test_count

    ! Final score
    print *, "Final Score: ", total_success, "/", total_tests, " tests passed"

contains

    subroutine generate_unique_set(set, n)
        type(test_pair), intent(out) :: set(:)
        integer, intent(in) :: n
        integer :: i, j
        logical :: is_unique

        do i = 1, n
            is_unique = .false.
            do while (.not. is_unique)
                call generate_key_value_pair(set(i))
                is_unique = .true.
                do j = 1, i-1
                    if (are_keys_equivalent(set(j)%key, set(i)%key)) then
                        is_unique = .false.
                        exit
                    end if
                end do
            end do
        end do
    end subroutine generate_unique_set

    subroutine generate_key_value_pair(pair)
        type(test_pair), intent(out) :: pair
        integer :: j
        real :: random_value

        call random_number(random_value)
        allocate(pair%key(int(random_value * (MAX_KEY_LENGTH - 1)) + 1))
        do j = 1, size(pair%key)
            call random_number(random_value)
            pair%key(j) = int(random_value * MAX_KEY_VALUE) + 1
        end do
        call random_number(random_value)
        pair%value = int(random_value * 1000) + 1
    end subroutine generate_key_value_pair

    function is_equivalent_to_any(pairs, n, key) result(equivalent)
        type(test_pair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        integer, intent(in) :: key(:)
        logical :: equivalent
        integer :: i

        equivalent = .false.
        do i = 1, n
            if (are_keys_equivalent(pairs(i)%key, key)) then
                equivalent = .true.
                return
            end if
        end do
    end function is_equivalent_to_any

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

    subroutine shuffle_array(arr)
        integer, intent(inout) :: arr(:)
        integer :: i, j, temp
        real :: r

        do i = size(arr), 2, -1
            call random_number(r)
            j = int(r * i) + 1
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
        end do
    end subroutine shuffle_array

end program test_dict