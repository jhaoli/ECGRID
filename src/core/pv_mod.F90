module pv_mod

  use const_mod
  use mesh_mod
  use namelist_mod
  use parallel_mod
  use state_mod

  implicit none

  private
  
  public calc_pv_on_vertex
  public calc_pv_on_edge_midpoint

contains

  subroutine calc_pv_on_vertex(state)
  
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j
    real(r8) half_cos_lat, full_cos_lat 

    mesh => state%mesh

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        ! cos_lat = (mesh%half_f(j+1) * mesh%full_cos_lat(j+1) + mesh%half_f(j) * mesh%full_cos_lat(j)) * 0.5_r8
        ! cos_lat = (mesh%full_cos_lat(j+1) + mesh%full_cos_lat(j)) * 0.5_r8
        half_cos_lat = mesh%half_cos_lat(j)
#ifdef V_POLE
        full_cos_lat = (mesh%full_cos_lat(j) + mesh%full_cos_lat(j-1)) * 0.5_r8
        state%pv(i,j) = ((state%v(i+1,j) - state%v(i,j)) / mesh%dlon - &
                         (state%u(i,j  ) * mesh%full_cos_lat(j  ) -&
                          state%u(i,j-1) * mesh%full_cos_lat(j-1)) / mesh%dlat +&
                         radius * mesh%half_f(j) * half_cos_lat ) /&
                         (radius * state%m_vtx(i,j) * full_cos_lat)
#else    
        full_cos_lat = (mesh%full_cos_lat(j) + mesh%full_cos_lat(j+1)) * 0.5_r8
        state%pv(i,j) = ((state%v(i+1,j) - state%v(i,j)) / mesh%dlon - &
                         (state%u(i,j+1) * mesh%full_cos_lat(j+1) -&
                          state%u(i,j  ) * mesh%full_cos_lat(j  )) / mesh%dlat + &
                         radius * mesh%half_f(j) * half_cos_lat ) / &
                         (radius * state%m_vtx(i,j) * full_cos_lat)
#endif
      end do
    end do
#ifdef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%half_lat_start_idx
      state%vor_sp = 0.0_r8
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%vor_sp = state%vor_sp - state%u(i,j)
      end do
      state%vor_sp = state%vor_sp / mesh%num_half_lon / (radius * mesh%dlat * 0.5_r8)
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%pv(i,j) = (state%vor_sp + mesh%half_f(j)) / state%m_vtx(i,j)
      end do
    endif

    if (mesh%has_north_pole()) then
      j = mesh%half_lat_end_idx
      state%vor_np = 0.0_r8
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%vor_np = state%vor_np + state%u(i,j-1)
      end do
      state%vor_np = state%vor_np / mesh%num_half_lon / (radius * mesh%dlat * 0.5_r8)
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%pv(i,j) = (state%vor_sp + mesh%half_f(j)) / state%m_vtx(i,j)
      end do
    endif
#endif
    call parallel_fill_halo(mesh, state%pv)

  end subroutine calc_pv_on_vertex

  subroutine calc_pv_on_edge_midpoint(state)
    
    type(state_type), intent(inout) :: state
    
    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

    do j = mesh%half_lat_start_idx, mesh%half_lat_end_idx
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%pv_lat(i,j) = 0.5_r8 * (state%pv(i-1,j) + state%pv(i,j))
      end do
    end do
    call parallel_fill_halo(mesh, state%pv_lat)

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
#ifdef V_POLE
        state%pv_lon(i,j) = 0.5_r8 * (state%pv(i,j) + state%pv(i,j+1))
#else
        state%pv_lon(i,j) = 0.5_r8 * (state%pv(i,j) + state%pv(i,j-1))
#endif
      end do
    end do
    call parallel_fill_halo(mesh, state%pv_lon)

  end subroutine calc_pv_on_edge_midpoint

end module pv_mod