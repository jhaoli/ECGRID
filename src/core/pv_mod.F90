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
    real(r8) cos_lat 

    mesh => state%mesh

    do j = mesh%half_lat_start_idx, mesh%half_lat_end_idx
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        ! cos_lat = (mesh%full_cos_lat(j+1) + mesh%full_cos_lat(j)) * 0.5_r8
        cos_lat = mesh%half_cos_lat(j)
        state%pv(i,j) =((state%v(i+1,j) - state%v(i,j)) / mesh%dlon - &
                        (state%u(i,j+1) * mesh%full_cos_lat(j+1) - state%u(i,j) * mesh%full_cos_lat(j)) / mesh%dlat + &
                        radius * mesh%half_f(j) * cos_lat ) / &
                        (radius * state%m_vtx(i,j))
      end do
    end do

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
        state%pv_lon(i,j) = 0.5_r8 * (state%pv(i,j) + state%pv(i,j-1))
      end do
    end do
    call parallel_fill_halo(mesh, state%pv_lon)

  end subroutine calc_pv_on_edge_midpoint

end module pv_mod