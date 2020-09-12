program ecgrid_swm_driver

  use log_mod
  use namelist_mod
  use static_mod
  use state_mod
  use ecgrid_mod
  use mountain_zonal_flow_test_mod
  use rossby_haurwitz_wave_test_mod
  use cross_pole_flow_test_mod
  use jet_zonal_flow_test_mod

  implicit none

  character(256) :: namelist_file_path
  integer i 
  
  interface
    subroutine set_initial_condition_interface(static, state)
      import static_type, state_type
      type(static_type), intent(inout) :: static
      type(state_type ), intent(inout) :: state
    end subroutine
  end interface
  procedure(set_initial_condition_interface), pointer :: set_initial_condition

  call get_command_argument(1, namelist_file_path)
  if (namelist_file_path == '') then
    call log_error('You should give a namelist file path!')
  end if 

  call parse_namelist(namelist_file_path)
  call ecgrid_init()

  select case (test_case)
  case ('mountain_zonal_flow')
    set_initial_condition => mountain_zonal_flow_test_set_initial_condition
  case ('rossby_haurwitz_wave')
    set_initial_condition => rossby_haurwitz_wave_test_set_initial_condition
  case ('cross_pole_flow')
    set_initial_condition => cross_pole_flow_test_set_initial_condition
  case ('jet_zonal_flow')
    set_initial_condition => jet_zonal_flow_test_set_initial_condition
  case default
    call log_error('Unknown test case ' // trim(test_case) // '!')
  end select

  call set_initial_condition(static, states(1))

  call ecgrid_run()

  call ecgrid_final()

end program ecgrid_swm_driver