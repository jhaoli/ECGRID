module history_mod

  use flogger
  use fiona
  use string
  use const_mod
  use namelist_mod
  use allocator_mod
  use time_mod, dt => dt_in_seconds
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod

  implicit none

  private
  public history_init
  public history_final
  public history_write_state
  public history_write_debug

  ! A-grid variables
  real(r8), allocatable, dimension(:,:) :: u
  real(r8), allocatable, dimension(:,:) :: v
  real(r8), allocatable, dimension(:,:) :: h
  real(r8), allocatable, dimension(:,:) :: hs

contains

  subroutine history_init()
    
    character(10) time_value, time_units
    real(r8) seconds

    if (history_interval(1) == 'N/A') call log_error('Parameter history_interval is not set!')
    if (case_name == 'N/A') call log_error('Parameter case_name is not set!')

    time_value = split_string(history_interval(1), ' ', 1)
    time_units = split_string(history_interval(1), ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      call log_error('Invalid history interval ' // trim(history_interval(1)) // '!')
    end select
    
    call fiona_init(time_units, start_time_str)

    call fiona_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name))
    call fiona_add_att('h0', 'time_step_size', dt)
    call fiona_add_dim('h0', 'time' , add_var=.true.)
    call fiona_add_dim('h0', 'lon'  , size=mesh%num_full_lon, add_var=.true.)
    call fiona_add_dim('h0', 'lat'  , size=mesh%num_full_lat, add_var=.true.)
    call fiona_add_dim('h0', 'ilon' , size=mesh%num_half_lon, add_var=.true.)
    call fiona_add_dim('h0', 'ilat' , size=mesh%num_half_lat, add_var=.true.)
    call fiona_add_var('h0', 'u'    , long_name='u wind component'           , units='m s-1' , dim_names=['lon ', 'lat ', 'time'])
    call fiona_add_var('h0', 'v'    , long_name='v wind component'           , units='m s-1' , dim_names=['lon ', 'lat ', 'time'])
    call fiona_add_var('h0', 'h'    , long_name='height'                     , units='m    ' , dim_names=['lon ', 'lat ', 'time'])
    call fiona_add_var('h0', 'hs'   , long_name='surface height'             , units='m    ' , dim_names=['lon ', 'lat ', 'time'])
    call fiona_add_var('h0', 'pv'   , long_name='potential vorticity'        , units='s-1  ' , dim_names=['ilon', 'ilat', 'time'])
    call fiona_add_var('h0', 'tm'   , long_name='total mass         '        , units='m    ' , dim_names=['time'])
    call fiona_add_var('h0', 'te'   , long_name='total energy       '        , units='m4 s-4', dim_names=['time'], data_type='real(8)')
    call fiona_add_var('h0', 'tpe'  , long_name='total potential enstrophy'  , units='m2 s-5', dim_names=['time'], data_type='real(8)')
    call fiona_add_var('h0', 'tpv'  , long_name='total potential vorticity'  , units='s-1   ', dim_names=['time'], data_type='real(8)')

    call fiona_create_dataset('h1', desc=case_desc, file_prefix=trim(case_name))
    call fiona_add_att('h1', 'time_step_size', dt)
    call fiona_add_dim('h1', 'time'     , add_var=.true.)
    call fiona_add_dim('h1', 'lon'      , size=mesh%num_full_lon, add_var=.true.)
    call fiona_add_dim('h1', 'lat'      , size=mesh%num_full_lat, add_var=.true.)
    call fiona_add_dim('h1', 'ilon'     , size=mesh%num_half_lon, add_var=.true.)
    call fiona_add_dim('h1', 'ilat'     , size=mesh%num_half_lat, add_var=.true.)
    call fiona_add_var('h1', 'qhv'      , long_name='nonliear zonal Coriolis force'                , units='m s-2'  , dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'qhu'      , long_name='nonliear meridional Coriolis force'           , units='m s-2'  , dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'dpedlon'  , long_name='zonal geopotential energy gradient force'     , units=''       , dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'dkedlon'  , long_name='zonal kinetic energy gradient force'          , units=''       , dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'dpedlat'  , long_name='meridional geopotential energy gradient force', units=''       , dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'dkedlat'  , long_name='meridional kinetic energy gradient force'     , units=''       , dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'dmfdlon'  , long_name='zonal mass flux divergence'                   , units=''       , dim_names=['lon ', 'lat ', 'time'])
    call fiona_add_var('h1', 'dmfdlat'  , long_name='meridional mass flux divergence'              , units=''       , dim_names=['lon ', 'lat ', 'time'])
    call fiona_add_var('h1', 'us'       , long_name='staggered u wind component'                   , units='m s-1'  , dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'vs'       , long_name='staggered v wind component'                   , units='m s-1'  , dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'hv'       , long_name='mass on vertices'                             , units='m'      , dim_names=['ilon', 'ilat', 'time'])
    call fiona_add_var('h1', 'mf_lon_n' , long_name='normal mass flux on U grid'                   , units='m2 s-1' , dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'mf_lon_t' , long_name='tangent mass flux on U grid'                  , units='m2 s-1' , dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'mf_lat_n' , long_name='normal mass flux on V grid'                   , units='m2 s-1' , dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'mf_lat_t' , long_name='tangent mass flux on V grid'                  , units='m2 s-1' , dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'pv_lon'   , long_name='pv on U grid'                                 , units='m-1 s-1', dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'pv_lat'   , long_name='pv on V grid'                                 , units='m-1 s-1', dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'ke'       , long_name='kinetic energy on cell grid'                  , units='',        dim_names=['lon ', 'lat ', 'time'])

    if (.not. allocated(u )) call allocate_array(mesh, u )
    if (.not. allocated(v )) call allocate_array(mesh, v )
    if (.not. allocated(h )) call allocate_array(mesh, h )
    if (.not. allocated(hs)) call allocate_array(mesh, hs)
    
    call time_add_alert('history_write', seconds=seconds)
    
    call log_notice('History module is initialized.')

  end subroutine history_init

  subroutine history_final()

    if (allocated(u )) deallocate(u )
    if (allocated(v )) deallocate(v )
    if (allocated(h )) deallocate(h )
    if (allocated(hs)) deallocate(hs)

  end subroutine history_final
  
  subroutine history_write_state(static, state)
    
    type(static_type), intent(in) :: static
    type(state_type ), intent(in) :: state

    integer i, j

    ! Convert wind from C grid to A grid 
    do j = state%mesh%full_lat_start_idx, state%mesh%full_lat_end_idx
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef V_POLE
        v(i,j) = 0.5_r8 * (state%v(i,j) + state%v(i,j+1))
#else
        v(i,j) = 0.5_r8 * (state%v(i,j) + state%v(i,j-1))
#endif
        u(i,j) = 0.5_r8 * (state%u(i,j) + state%u(i-1,j))
        h(i,j) = (state%gd(i,j) + static%ghs(i,j)) / g
        hs(i,j) = static%ghs(i,j) / g
      end do
    end do
    
    call fiona_start_output('h0', elapsed_seconds, new_file=time_step == 0)
    call fiona_output('h0', 'lon' , mesh%full_lon_deg(1:mesh%num_full_lon))
    call fiona_output('h0', 'lat' , mesh%full_lat_deg(1:mesh%num_full_lat))
    call fiona_output('h0', 'ilon', mesh%half_lon_deg(1:mesh%num_half_lon))
    call fiona_output('h0', 'ilat', mesh%half_lat_deg(1:mesh%num_half_lat))
    call fiona_output('h0', 'u'   , u                (1:mesh%num_full_lon,1:mesh%num_full_lat))
    call fiona_output('h0', 'v'   , v                (1:mesh%num_full_lon,1:mesh%num_full_lat))
    call fiona_output('h0', 'h'   , h                (1:mesh%num_full_lon,1:mesh%num_full_lat))
    call fiona_output('h0', 'hs'  , hs               (1:mesh%num_full_lon,1:mesh%num_full_lat))
    call fiona_output('h0', 'pv'  , state%pv         (1:mesh%num_half_lon,1:mesh%num_half_lat))
    call fiona_output('h0', 'tm'  , state%total_m)
    call fiona_output('h0', 'te'  , state%total_e)
    call fiona_output('h0', 'tpe' , state%total_pe)
    call fiona_output('h0', 'tpv' , state%total_pv)
    call fiona_end_output('h0')

  end subroutine history_write_state

  subroutine history_write_debug(static, state, tend)

    type(static_type), intent(in) :: static
    type(state_type ), intent(in) :: state
    type(tend_type  ), intent(in) :: tend

    call fiona_start_output('h1', elapsed_seconds, new_file=time_step == 0)
    call fiona_output('h1', 'lon'      , mesh%full_lon_deg(1:mesh%num_full_lon))
    call fiona_output('h1', 'lat'      , mesh%full_lat_deg(1:mesh%num_full_lat))
    call fiona_output('h1', 'ilon'     , mesh%half_lon_deg(1:mesh%num_half_lon))
    call fiona_output('h1', 'ilat'     , mesh%half_lat_deg(1:mesh%num_half_lat))
    call fiona_output('h1', 'qhv'      , tend%qhv         (1:mesh%num_half_lon,1:mesh%num_full_lat))
    call fiona_output('h1', 'qhu'      , tend%qhu         (1:mesh%num_full_lon,1:mesh%num_half_lat))
    call fiona_output('h1', 'dpedlon'  , tend%dpedlon     (1:mesh%num_half_lon,1:mesh%num_full_lat))
    call fiona_output('h1', 'dkedlon'  , tend%dkedlon     (1:mesh%num_half_lon,1:mesh%num_full_lat))
    call fiona_output('h1', 'dpedlat'  , tend%dpedlat     (1:mesh%num_half_lon,1:mesh%num_full_lat))
    call fiona_output('h1', 'dkedlat'  , tend%dkedlat     (1:mesh%num_half_lon,1:mesh%num_full_lat))
    call fiona_output('h1', 'dmfdlon'  , tend%dmfdlon     (1:mesh%num_full_lon,1:mesh%num_full_lat))
    call fiona_output('h1', 'dmfdlat'  , tend%dmfdlat     (1:mesh%num_full_lon,1:mesh%num_full_lat))
    call fiona_output('h1', 'us'       , state%u          (1:mesh%num_half_lon,1:mesh%num_full_lat))
    call fiona_output('h1', 'vs'       , state%v          (1:mesh%num_full_lon,1:mesh%num_half_lat))
    call fiona_output('h1', 'hv'       , state%m_vtx      (1:mesh%num_half_lon,1:mesh%num_half_lat))
    call fiona_output('h1', 'mf_lon_n' , state%mf_lon_n   (1:mesh%num_half_lon,1:mesh%num_full_lat))
    call fiona_output('h1', 'mf_lat_n' , state%mf_lat_n   (1:mesh%num_full_lon,1:mesh%num_half_lat))
    call fiona_output('h1', 'pv_lon'   , state%pv_lon     (1:mesh%num_half_lon,1:mesh%num_full_lat))
    call fiona_output('h1', 'pv_lat'   , state%pv_lat     (1:mesh%num_full_lon,1:mesh%num_half_lat))
    call fiona_output('h1', 'ke'       , state%ke         (1:mesh%num_full_lon,1:mesh%num_full_lat))
    call fiona_end_output('h1')

  end subroutine history_write_debug

end module history_mod