module namelist_mod
  use const_mod
  use time_mod, start_time => start_time_array, end_time => end_time_array

  implicit none

  character(256) :: case_desc = 'N/A'
  character(256) :: case_name = 'N/A'
  character(30)  :: test_case = 'N/A'
  character(30)  :: history_interval(1) = 'N/A'

  integer :: num_lon
  integer :: num_lat


  character(30) :: time_scheme  = 'pc2'
  integer       :: coriolis_scheme = 2

  namelist /ecgrid_swm_control/ &
    case_name                 , &
    case_desc                 , &
    test_case                 , &
    num_lon                   , &
    num_lat                   , &
    dt_in_seconds             , &
    run_days                  , &
    history_interval          , &
    time_scheme               , &
    coriolis_scheme              

contains
  
  subroutine parse_namelist(file_path)
    
    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=ecgrid_swm_control)
    close(10) 

  end subroutine parse_namelist

end module namelist_mod