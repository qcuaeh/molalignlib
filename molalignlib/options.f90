module options
implicit none

! Flags
logical :: random_flag
logical :: mirror_flag
logical :: label_flag
logical :: iterate_flag
logical :: heavy_flag
logical :: mass_flag
logical :: align_flag
logical :: remap_flag
logical :: stoch_flag
logical :: adaptive_flag
logical :: coords_flag
logical :: stdin_flag
logical :: stats_flag
logical :: tree_flag
logical :: bond_flag
logical :: mapping_flag
logical :: fileout_flag

! Bounds
integer :: count_thres
integer :: max_trials
integer :: num_records

end module
