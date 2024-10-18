program test_dict_part1
    use dict_mod
    implicit none

    integer, parameter :: NUM_TEST_PAIRS = 200
    integer, parameter :: NUM_TEST_CASES = 100
    integer, parameter :: MAX_KEY_LENGTH = 16
    integer, parameter :: MAX_KEY_VALUE = 100
    real, parameter :: DICT_SIZE_FACTOR = 1.0

    type test_pair
        integer, allocatable :: key(:)
        integer :: value
    end type test_pair

    type(dict) :: d
    type(test_pair) :: first_set(NUM_TEST_PAIRS), second_set(NUM_TEST_PAIRS)
    type(test_pair) :: shuffled_first_set(NUM_TEST_PAIRS)
    integer :: i, test_count, success_count, total_success, total_tests
    integer :: expected_dict_size
    real :: random_value

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
    d = create_dict(int(NUM_TEST_PAIRS * DICT_SIZE_FACTOR))
    print '(A,I0)', "Actual dictionary capacity: ", d%actual_capacity

    total_success = 0
    total_tests = 0
    expected_dict_size = 0

    ! Test 1: Put and existence of keys (first set)
    print *, "Test 1: Put and existence of keys (first set)"
    test_count = 0
    success_count = 0
    do i = 1, NUM_TEST_CASES
        test_count = test_count + 1
        ! Uncomment the following line to print generated keys
        ! print *, "Putting key:", shuffled_first_set(i)%key
        call d%put(shuffled_first_set(i)%key, shuffled_first_set(i)%value)
        expected_dict_size = expected_dict_size + 1
        if (d%exists(first_set(i)%key)) then
            success_count = success_count + 1
        end if
    end do
    print '(A,I0,A,I0,A)', "Test 1 Summary: ", success_count, "/", test_count, " tests passed"
    total_success = total_success + success_count
    total_tests = total_tests + test_count

    ! Test 2: Non-existence of keys (second set)
    print *, "Test 2: Non-existence of keys (second set)"
    test_count = 0
    success_count = 0
    do i = 1, NUM_TEST_CASES
        test_count = test_count + 1
        ! Uncomment the following line to print generated keys
        ! print *, "Checking non-existence of key:", second_set(i)%key
        if (.not. d%exists(second_set(i)%key)) then
            success_count = success_count + 1
        end if
    end do
    print '(A,I0,A,I0,A)', "Test 2 Summary: ", success_count, "/", test_count, " tests passed"
    total_success = total_success + success_count
    total_tests = total_tests + test_count

    ! Test 3: Get of values
    print *, "Test 3: Get of values"
    test_count = 0
    success_count = 0
    do i = 1, NUM_TEST_CASES
        test_count = test_count + 1
        ! Uncomment the following line to print generated keys
        ! print *, "Getting value for key:", first_set(i)%key
        if (d%get(first_set(i)%key) == first_set(i)%value) then
            success_count = success_count + 1
        end if
    end do
    print '(A,I0,A,I0,A)', "Test 3 Summary: ", success_count, "/", test_count, " tests passed"
    total_success = total_success + success_count
    total_tests = total_tests + test_count

    ! Test 4: Overwriting values
    print *, "Test 4: Overwriting values"
    test_count = 0
    success_count = 0
    do i = 1, NUM_TEST_CASES
        test_count = test_count + 1
        ! Uncomment the following line to print generated keys
        ! print *, "Overwriting value for key:", first_set(i)%key
        call d%put(first_set(i)%key, first_set(i)%value * 2)
        if (d%get(first_set(i)%key) == first_set(i)%value * 2) then
            success_count = success_count + 1
        end if
    end do
    print '(A,I0,A,I0,A)', "Test 4 Summary: ", success_count, "/", test_count, " tests passed"
    total_success = total_success + success_count
    total_tests = total_tests + test_count

    ! Test 5: Reset function
    print *, "Test 5: Reset function"
    call d%reset()
    test_count = 2
    success_count = 0
    if (sum(d%occupations) == 0) then
        success_count = success_count + 1
        print *, "Test 5.1: PASS (Dictionary count reset to 0)"
    else
        print *, "Test 5.1: FAIL (Dictionary count not reset to 0)"
    end if
    if (.not. d%exists(first_set(1)%key)) then
        success_count = success_count + 1
        print *, "Test 5.2: PASS (Key no longer exists after reset)"
    else
        print *, "Test 5.2: FAIL (Key still exists after reset)"
    end if
    print '(A,I0,A,I0,A)', "Test 5 Summary: ", success_count, "/", test_count, " tests passed"
    total_success = total_success + success_count
    total_tests = total_tests + test_count

    ! Final score
    print '(A,I0,A,I0,A)', "Final Score: ", total_success, "/", total_tests, " tests passed"

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

end program test_dict_part1