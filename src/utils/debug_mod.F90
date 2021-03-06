module debug_mod

  use flogger
  use const_mod
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod

  implicit none
  
  private

  public debug_check_space_operators

contains
  
  subroutine debug_check_space_operators(static, state, tend)
    
    type(static_type), intent(in) :: static
    type(state_type ), intent(in) :: state
    type(tend_type  ), intent(in) :: tend

    integer i, j
    type(mesh_type), pointer :: mesh

    real(r8) ip_cf
    real(r8) ip_dke
    real(r8) ip_dpe
    real(r8) ip_dmf

    mesh => state%mesh
    ip_cf  = 0.0_r8
    ip_dke = 0.0_r8
    ip_dpe = 0.0_r8
    ip_dmf = 0.0_r8

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        ip_cf = ip_cf + tend%qhv(i,j) * state%mf_lon_n(i,j) * radius**2 * mesh%full_cos_lat(j) * mesh%dlon * mesh%dlat
        ip_dke = ip_dke + tend%dkedlon(i,j) * state%mf_lon_n(i,j) * radius**2 * mesh%full_cos_lat(j) * mesh%dlon * mesh%dlat
        ip_dpe = ip_dpe + tend%dpedlon(i,j) * state%mf_lon_n(i,j) * radius**2 * mesh%full_cos_lat(j) * mesh%dlon * mesh%dlat
      end do
    end do
    ! print*, '1', ip_dke, ip_dpe
    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        ip_cf = ip_cf - tend%qhu(i,j) * state%mf_lat_n(i,j) * radius**2 * mesh%half_cos_lat(j) * mesh%dlon * mesh%dlat
        ip_dke = ip_dke + tend%dkedlat(i,j) * state%mf_lat_n(i,j) * radius**2 * mesh%half_cos_lat(j) * mesh%dlon * mesh%dlat
        ip_dpe = ip_dpe + tend%dpedlat(i,j) * state%mf_lat_n(i,j) * radius**2 * mesh%half_cos_lat(j) * mesh%dlon * mesh%dlat
      end do
    end do
    ! ip_dke = 0.0_r8
    ! ip_dpe = 0.0_r8
    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        ! ip_dke = ip_dke - tend%dmfdlon(i,j) * state%ke(i,j) * radius**2 * mesh%full_cos_lat(j) * mesh%dlon * mesh%dlat
        ip_dke = ip_dke - (tend%dmfdlon(i,j) + tend%dmfdlat(i,j)) * state%ke(i,j) * radius**2 * mesh%full_cos_lat(j) * mesh%dlon * mesh%dlat
        ! ip_dpe = ip_dpe - tend%dmfdlon(i,j) * (state%gd(i,j) + static%ghs(i,j)) * radius**2 * mesh%full_cos_lat(j) * mesh%dlon * mesh%dlat
        ip_dpe = ip_dpe - (tend%dmfdlon(i,j) + tend%dmfdlat(i,j)) * (state%gd(i,j) + static%ghs(i,j)) * radius**2 * mesh%full_cos_lat(j) * mesh%dlon * mesh%dlat
        ip_dmf = ip_dmf + (tend%dmfdlon(i,j) + tend%dmfdlat(i,j)) * radius**2 * mesh%full_cos_lat(j) * mesh%dlon * mesh%dlat
      end do
    end do
    ! print*, '2', ip_dke, ip_dpe

    ! print*, ip_cf  / radius**2, ip_dmf / radius**2, ip_dke / radius**2 , ip_dpe / radius**2
    call log_add_diag('cf' , ip_cf  / radius**2)
    call log_add_diag('dke', ip_dke / radius**2)
    call log_add_diag('dpe', ip_dpe / radius**2)
    call log_add_diag('dmf', ip_dmf / radius**2)

  end subroutine debug_check_space_operators
end module debug_mod