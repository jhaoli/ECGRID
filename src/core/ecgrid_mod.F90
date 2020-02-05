module ecgrid_mod
  
  use flogger
  use const_mod
  use namelist_mod
  use parallel_mod
  use time_mod, dt => dt_in_seconds, old => old_time_idx, new => new_time_idx
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod
  use history_mod
  use operators_mod
  use debug_mod

  implicit none

  private

  public ecgrid_init
  public ecgrid_run
  public ecgrid_final

  integer, parameter :: all_pass = 0

  interface
    subroutine integrator_interface(dt, static, tends, states, old, new, pass)     
      import r8, static_type, tend_type, state_type
      real(r8)         , intent(in   ) :: dt
      type(static_type), intent(in   ) :: static
      type(tend_type  ), intent(inout) :: tends (:)
      type(state_type ), intent(inout) :: states(:)
      integer          , intent(in   ) :: old
      integer          , intent(in   ) :: new
      integer          , intent(in   ) :: pass
    end subroutine integrator_interface

    subroutine splitter_interface(dt, static, tends, states)
      import r8, static_type, tend_type, state_type
      real(r8)         , intent(in   ) :: dt
      type(static_type), intent(in   ) :: static
      type(tend_type  ), intent(inout) :: tends (:)
      type(state_type ), intent(inout) :: states(:)
    end subroutine splitter_interface
  end interface

  procedure(integrator_interface), pointer :: integrator
  procedure(splitter_interface), pointer :: splitter

contains
  
  subroutine ecgrid_init()
    
    call log_init()
    call time_init()
    call mesh%init(num_lon, num_lat)
    call state_init_root()
    call static_init_root()
    call tend_init_root()
    call history_init()

    select case (time_scheme)
    case ('pc2')
      integrator => predict_correct
    case default
      integrator => predict_correct
      call log_notice('Use pc2 integrator.')
    end select
   
    select case (split_scheme)
    ! case ('pc2')
    !   splitter => csp2_splitting
    case default
      splitter => no_splitting
      call log_notice('No fast-slow split.')
    end select
    
    call time_add_alert('print', days=1.0_r8)

  end subroutine ecgrid_init
  
  subroutine ecgrid_run()

    call operators_prepare(states(old))
    call diagnose(states(old))
    call output(states(old), tends(old))

    do while (.not. time_is_finished())
      call time_integrate(dt, static, tends, states)
      if (time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
      call time_advance()
      call operators_prepare(states(old))
      call diagnose(states(old))
      call output(states(old), tends(old))
    end do

  end subroutine ecgrid_run

  subroutine ecgrid_final()
  
    call parallel_final()

  end subroutine ecgrid_final

  subroutine output(state, tend)
    
    type(state_type), intent(in) :: state
    type(tend_type ), intent(in) :: tend

    if (time_is_alerted('history_write')) then
      call history_write_state(static, state)
      ! call history_write_debug(static, state, tend)
    end if
  end subroutine output

  subroutine diagnose(state)
   
    type(state_type), intent(inout) :: state
    
    type(mesh_type), pointer :: mesh
    integer i, j
    real(r8) vor

    mesh => state%mesh

    state%total_m = 0.0_r8
    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%total_m = state%total_m + state%gd(i,j) * mesh%full_cos_lat(j) * radius**2 * mesh%dlon * mesh%dlat
      end do
    end do

    call log_add_diag('total_m' , state%total_m)
 
    state%total_e = 0.0_r8
    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%total_e = state%total_e + (state%gd(i,j) * state%ke(i,j) + state%gd(i,j) * (0.5_r8 * state%gd(i,j) + static%ghs(i,j))) * mesh%full_cos_lat(j) * radius**2 * mesh%dlon * mesh%dlat 
      end do
    end do
    call log_add_diag('total_e' , state%total_e)

    state%total_pv = 0.0_r8
    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
#ifdef V_POLE
        vor = (state%v(i+1,j) - state%v(i,j)) / mesh%dlon - &
              (state%u(i,j  ) * mesh%full_cos_lat(j  ) -&
               state%u(i,j-1) * mesh%full_cos_lat(j-1)) / mesh%dlat                       
        ! state%total_pv = state%total_pv + state%m_vtx(i,j) * state%pv(i,j) * mesh%half_cos_lat(j) * radius**2 * 
#else
        vor = (state%v(i+1,j) - state%v(i,j)) / mesh%dlon - &
              (state%u(i,j+1) * mesh%full_cos_lat(j+1) -&
               state%u(i,j  ) * mesh%full_cos_lat(j  )) / mesh%dlat
#endif
        state%total_pv = state%total_pv + vor * radius * mesh%dlon * mesh%dlat
      end do
    end do
    ! call log_add_diag('total_pv' , state%total_pv)

    state%total_pe = 0.0_r8
    do j = mesh%half_lat_start_idx, mesh%half_lat_end_idx
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%total_pe = state%total_pe + 0.5_r8 * state%m_vtx(i,j) * state%pv(i,j)**2 * mesh%half_cos_lat(j) * radius**2 * mesh%dlon * mesh%dlat
      end do
    end do
    call log_add_diag('total_pe' , state%total_pe)

  end subroutine diagnose

  subroutine space_operators(static, state, tend, dt, pass)
    
    type(static_type), intent(in   ) :: static
    type(state_type) , intent(inout) :: state
    type(tend_type)  , intent(inout) :: tend
    real(r8)         , intent(in   ) :: dt 
    integer          , intent(in   ) :: pass

    type(mesh_type), pointer :: mesh 
    integer i, j

    call operators_prepare(state)

    mesh => state%mesh
    select case(pass)
    case (all_pass)
      call calc_qhu_qhv_2(state, tend, dt)
      call calc_dkedlon_dkedlat(state, tend, dt)
      call calc_dpedlon_dpedlat(static, state, tend, dt)
      call calc_dmfdlon_dmfdlat(state, tend, dt)
  
      do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
        do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
          tend%du(i,j) =   tend%qhv(i,j) - tend%dpedlon(i,j) - tend%dkedlon(i,j)
        end do
      end do

      do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
        do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
          tend%dv(i,j) = - tend%qhu(i,j) - tend%dpedlat(i,j) - tend%dkedlat(i,j)
        end do
      end do

      do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
        do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
          tend%dgd(i,j) = - tend%dmfdlon(i,j) - tend%dmfdlat(i,j)
        end do
      end do
    end select
  
    ! call debug_check_space_operators(static, state, tend)

  end subroutine space_operators

  subroutine time_integrate(dt, static, tends, states)
    
    real(r8)         , intent(in   ) :: dt
    type(static_type), intent(in   ) :: static
    type(tend_type  ), intent(inout) :: tends(:)
    type(state_type ), intent(inout) :: states(:)

    call splitter(dt, static, tends, states)
   
  end subroutine time_integrate

  subroutine no_splitting(dt, static, tends, states)
  
    real(r8)         , intent(in)    :: dt
    type(static_type), intent(in)    :: static
    type(tend_type  ), intent(inout) :: tends (:)
    type(state_type ), intent(inout) :: states(:) 

    call integrator(dt, static, tends, states, old, new, all_pass)

  end subroutine no_splitting

  subroutine predict_correct(dt, static, tends, states, old, new, pass)
    
    real(r8)         , intent(in   ) :: dt
    type(static_type), intent(in   ) :: static
    type(tend_type  ), intent(inout) :: tends (:)
    type(state_type ), intent(inout) :: states(:)
    integer          , intent(in   ) :: old
    integer          , intent(in   ) :: new
    integer          , intent(in   ) :: pass

    ! Do first predict step.
    call space_operators(static, states(old), tends(old), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, tends(old), states(old), states(new))

    ! Do second predict step.
    call space_operators(static, states(new), tends(old), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, tends(old), states(old), states(new))

    ! Do correct stepe
    call space_operators(static, states(new), tends(new),          dt, pass)
    call update_state(         dt, tends(new), states(old), states(new))

  end subroutine predict_correct
  
  subroutine update_state(dt, tend, old_state, new_state)
    
    real(r8)        , intent(in   ) :: dt
    type(tend_type ), intent(in   ) :: tend
    type(state_type), intent(in   ) :: old_state
    type(state_type), intent(inout) :: new_state

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => old_state%mesh

    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
    	do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        new_state%gd(i,j) = old_state%gd(i,j) + dt * tend%dgd(i,j)
      end do
    end do
    call parallel_fill_halo(mesh, new_state%gd)

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        new_state%u(i,j) = old_state%u(i,j) + dt * tend%du(i,j)
      end do
    end do
    call parallel_fill_halo(mesh, new_state%u)

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        new_state%v(i,j) = old_state%v(i,j) + dt * tend%dv(i,j)
      end do
    end do
    call parallel_fill_halo(mesh, new_state%v)

  end subroutine update_state

end module ecgrid_mod