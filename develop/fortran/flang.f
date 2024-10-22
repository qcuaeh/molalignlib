integer, pointer, rank=1 : integer_array_pointer
(integer, pointer), rank=1 : integer_pointer_array
type=Permutation, pointer, rank=1 : permutation_array_pointer
(type=Permutation, pointer), rank=1 : permutation_pointer_array
(integer, allocatable), allocatable : integer_nested_list ! List of lists of integers
((integer, allocatable), allocatable), allocatable : integer_more_nested_list ! List of lists of lists of integers
integer, kind=int64, allocatable, rank=2 : int64_array2d
integer, kind=int64, rank=1, dimensions=[2] : int64_array_3x4

allocate( int64_array2d, [3, 4])