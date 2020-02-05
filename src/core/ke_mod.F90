module ke_mod
  
  use const_mod
  use mesh_mod
  use namelist_mod
  use parallel_mod
  use state_mod

  implicit none

  private

  public calc_ke_on_cell

contains

  subroutine calc_ke_on_cell(state)
    
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh 
    integer i, j
    real pole

    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%ke(i,j) = ((state%u(i,j  )**2 + state%u(i-1,j)**2) + &
#ifdef V_POLE
                         (state%v(i,j+1)**2 * mesh%half_cos_lat(j+1) +&
                          state%v(i,j  )**2 * mesh%half_cos_lat(j  )) / mesh%full_cos_lat(j)) * 0.25_r8                      
#else          
                         (state%v(i,j-1)**2 * mesh%half_cos_lat(j-1) +&
                          state%v(i,j  )**2 * mesh%half_cos_lat(j  )) / mesh%full_cos_lat(j)) * 0.25_r8                  
#endif      
      end do
    end do

#ifndef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_start_idx
      pole = 0.0_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole + state%v(i,j)**2
      end do
      pole = pole / mesh%num_full_lon
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%ke(i,j) = pole
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_end_idx
      pole = 0.0_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole + state%v(i,j-1)**2
      end do
      pole = pole / mesh%num_full_lon
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%ke(i,j) = pole
      end do
    end if  
#endif
    call parallel_fill_halo(mesh, state%ke)

  end subroutine calc_ke_on_cell
end module ke_mod