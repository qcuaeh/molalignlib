module options
implicit none

! Flags
logical :: random_flag
logical :: mirror_flag
logical :: label_flag
logical :: heavy_flag
logical :: mass_flag
logical :: align_flag
logical :: remap_flag
logical :: stoch_flag
logical :: adaptive_flag
logical :: coords_flag
logical :: stats_flag
logical :: tree_flag
logical :: rebond_flag
logical :: atomorder_flag
logical :: full_flag

! Bounds
integer :: max_trials
integer :: num_records
integer :: ato_thres
integer :: iso_thres
integer :: confo_thres

end module
