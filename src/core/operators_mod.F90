module operators_mod

  use flogger
  use const_mod
  use namelist_mod
  use mesh_mod
  use static_mod
  use state_mod
  use tend_mod
  use parallel_mod
  use pv_mod
  use ke_mod

  implicit none

  private

  public operators_prepare_all
  public calc_qhu_qhv
  public calc_qhu_qhv_2
  public calc_dkedlon_dkedlat
  public calc_dpedlon_dpedlat
  public calc_dmfdlon_dmfdlat

contains

  subroutine operators_prepare_all(state)

    type(state_type), intent(inout) :: state

    call calc_m_lon_m_lat(state)
    call calc_m_vtx(state)
    call calc_mf_lon_n_mf_lat_n(state)
    call calc_pv_on_vertex(state)
    call calc_ke_on_cell(state)

  end subroutine operators_prepare_all

  subroutine calc_m_lon_m_lat(state)
    
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%m_lon(i,j) = (state%gd(i,j) + state%gd(i+1,j)) * 0.5_r8
      end do
    end do 

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef V_POLE
        state%m_lat(i,j) = (state%gd(i,j) + state%gd(i,j-1)) * 0.5_r8
#else
        state%m_lat(i,j) = (state%gd(i,j) + state%gd(i,j+1)) * 0.5_r8
#endif
      end do
    end do

  end subroutine calc_m_lon_m_lat

  subroutine calc_m_vtx(state)
    
    type(state_type), intent(inout) :: state
    type(mesh_type), pointer :: mesh

    integer i, j
    real(r8) pole
 
    mesh => state%mesh

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%m_vtx(i,j) = (state%gd(i,j  ) + state%gd(i+1,j  ) +&
#ifdef V_POLE
                            state%gd(i,j-1) + state%gd(i+1,j-1)) * 0.25_r8
#else
                            state%gd(i,j+1) + state%gd(i+1,j+1)) * 0.25_r8
#endif
      end do
    end do
#ifdef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%half_lat_start_idx
      pole = 0.0_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole + state%gd(i,j)
      end do
      pole = pole / mesh%num_half_lon
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%m_vtx(i,j) = pole
      end do
    endif
    if(mesh%has_north_pole()) then
      j = mesh%half_lat_end_idx
      pole = 0.0_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole + state%gd(i,j-1)
      end do
      pole = pole / mesh%num_half_lon
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%m_vtx(i,j) = pole
      end do
    endif
#endif

  end subroutine calc_m_vtx
 
  subroutine calc_mf_lon_n_mf_lat_n(state)
  
    type(state_type), intent(inout) :: state
    type(mesh_type), pointer :: mesh
    integer i, j
    
    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%mf_lon_n(i,j) = state%m_lon(i,j) * state%u(i,j)
      end do
    end do
    call parallel_fill_halo(state%mesh, state%mf_lon_n)
#ifndef V_POLE
    !! calculate the mass flux on u-grid using Gauss elimination method
    if (mesh%has_south_pole()) then
      ! call calc_mf_lon_n_south_pole(state)
    end if 
    if (mesh%has_north_pole()) then
      ! call calc_mf_lon_n_north_pole(state)
    end if
#endif

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%mf_lat_n(i,j) = state%m_lat(i,j) * state%v(i,j)
      end do
    end do
    call parallel_fill_halo(state%mesh, state%mf_lat_n)

  end subroutine calc_mf_lon_n_mf_lat_n

  subroutine calc_mf_lon_n_south_pole(state)
    
    type(state_type), intent(inout) :: state
    type(mesh_type), pointer ::  mesh
    real(r8) pole_const
    integer  lda, ldb
    integer  nrhs, info
    integer  ipiv(state%mesh%num_half_lon)
    real(r8) coef(state%mesh%num_half_lon, state%mesh%num_half_lon), rhs(state%mesh%num_half_lon)
    integer i, j, k
    
    mesh => state%mesh

    lda = mesh%num_half_lon
    ldb = mesh%num_half_lon
    nrhs = 1

    j = mesh%full_lat_start_idx
    pole_const = 0.0_r8
    do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
      pole_const = pole_const + state%mf_lat_n(i,j) * mesh%half_cos_lat(j)
    end do
    pole_const = pole_const / mesh%num_full_lon * mesh%dlon / mesh%dlat
    
    do k = 1, mesh%num_half_lon
   	  do i = 1, mesh%num_half_lon
        if( i == k) then 
          coef(i,k) = -1.0_r8
        else if(i-k == 1) then 
          coef(i,k) = 1.0_r8
        else
          coef(i,k) = 0.0_r8
        end if 
      end do
    end do 
    coef(1,mesh%num_half_lon)=1.0_r8

    do i = 1, mesh%num_half_lon
      rhs(i) = pole_const - state%mf_lat_n(i,j) * mesh%half_cos_lat(j) * mesh%dlon * 2.0 / mesh%dlat
    end do
    call DGESV(state%mesh%num_half_lon, nrhs, coef, lda, ipiv, rhs, ldb, info)

    ! write(*,*) "Solution:"
    ! write(*,'(f8.3)') rhs
    ! write(*,*) "INFO=", info
    
    do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
      state%mf_lon_n(i,j) = rhs(i)
    end do
    call parallel_fill_halo(state%mesh, state%mf_lon_n(:,j))
    
  end subroutine calc_mf_lon_n_south_pole

  subroutine calc_mf_lon_n_north_pole(state)
    
    type(state_type), intent(inout) :: state
    type(mesh_type), pointer :: mesh
    real(r8) pole_const
    integer  lda, ldb
    integer  nrhs, info
    integer  ipiv(state%mesh%num_half_lon)
    real(r8) coef(state%mesh%num_half_lon, state%mesh%num_half_lon), rhs(state%mesh%num_half_lon)
    integer i, j, k
    
    mesh => state%mesh

    lda = mesh%num_half_lon
    ldb = mesh%num_half_lon
    nrhs = 1

    j = mesh%full_lat_end_idx
    pole_const = 0.0_r8
    do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
      pole_const = pole_const - state%mf_lat_n(i,j-1) * mesh%half_cos_lat(j-1)
    end do
    pole_const = pole_const / mesh%num_full_lon * mesh%dlon / mesh%dlat
    
    do k = 1, mesh%num_half_lon
   	  do i = 1, mesh%num_half_lon
        if( i == k) then 
          coef(i,k) = -1.0_r8
        else if(i-j == 1) then 
          coef(i,k) = 1.0_r8
        else
          coef(i,k) = 0.0_r8
        end if 
      end do
    end do 
    coef(1,mesh%num_half_lon)=1.0_r8

    do i = 1, mesh%num_half_lon
      rhs(i) = pole_const + state%mf_lat_n(i,j-1) * mesh%half_cos_lat(j-1) * mesh%dlon * 2.0 / mesh%dlat
    end do
    call DGESV(state%mesh%num_half_lon, nrhs, coef, lda, ipiv, rhs, ldb, info)

    ! write(*,*) "Solution:"
    ! write(*,'(f8.3)') rhs
    ! write(*,*) "INFO=", info
    ! stop

    do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
      state%mf_lon_n(i, j) = rhs(i)
    end do
    call parallel_fill_halo(state%mesh, state%mf_lon_n(:,j))

  end subroutine calc_mf_lon_n_north_pole

  subroutine calc_qhu_qhv(state, tend, dt)
    
    type(state_type), intent(inout) :: state
    type(tend_type ), intent(inout) :: tend
    real(r8)        , intent(in   ) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

    call calc_pv_on_edge_midpoint(state)

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        tend%qhv(i,j) = state%pv_lon(i,j) * ((state%mf_lat_n(i,j  ) + state%mf_lat_n(i+1,j  )) * mesh%half_cos_lat(j  ) + &
                                             (state%mf_lat_n(i,j-1) + state%mf_lat_n(i+1,j-1)) * mesh%half_cos_lat(j-1)) * 0.25_r8 / &
                        mesh%full_cos_lat(j)
      end do
    end do

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
    	do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef V_POLE
        tend%qhu(i,j) = state%pv_lat(i,j) * (state%mf_lon_n(i-1,j  ) + state%mf_lon_n(i,j  ) +&
                                             state%mf_lon_n(i-1,j-1) + state%mf_lon_n(i,j-1)) * 0.25_r8
#else
        tend%qhu(i,j) = state%pv_lat(i,j) * (state%mf_lon_n(i-1,j  ) + state%mf_lon_n(i,j  ) +&
                                             state%mf_lon_n(i-1,j+1) + state%mf_lon_n(i,j+1)) * 0.25_r8
#endif
      end do
    end do
  end subroutine calc_qhu_qhv

  subroutine calc_qhu_qhv_2(state, tend, dt)
    
    type(state_type), intent(inout) :: state
    type(tend_type ), intent(inout) :: tend
    real(r8)        , intent(in   ) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

#ifndef V_POLE
    ! mirror potential vorticity
    j = mesh%half_lat_end_idx
    state%pv(:,j+1) = state%pv(:,j)
 
    j = mesh%half_lat_start_idx
    state%pv(:,j-1) = state%pv(:,j)
#endif

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
#ifdef V_POLE
        tend%qhv(i,j) = ((state%mf_lat_n(i  ,j  ) * mesh%half_cos_lat(j  ) * (state%pv(i-1,j  ) + state%pv(i,j  ) + state%pv(i,j+1)) / 3.0_r8 + &
                          state%mf_lat_n(i+1,j  ) * mesh%half_cos_lat(j  ) * (state%pv(i+1,j  ) + state%pv(i,j  ) + state%pv(i,j+1)) / 3.0_r8 + &
                          state%mf_lat_n(i  ,j+1) * mesh%half_cos_lat(j+1) * (state%pv(i-1,j+1) + state%pv(i,j+1) + state%pv(i,j  )) / 3.0_r8 + &
                          state%mf_lat_n(i+1,j+1) * mesh%half_cos_lat(j+1) * (state%pv(i+1,j+1) + state%pv(i,j+1) + state%pv(i,j  )) / 3.0_r8)) * 0.25_r8 / &
#else
        tend%qhv(i,j) = ((state%mf_lat_n(i  ,j  ) * mesh%half_cos_lat(j  ) * (state%pv(i-1,j  ) + state%pv(i,j  ) + state%pv(i,j-1)) / 3.0_r8 + &
                          state%mf_lat_n(i+1,j  ) * mesh%half_cos_lat(j  ) * (state%pv(i+1,j  ) + state%pv(i,j  ) + state%pv(i,j-1)) / 3.0_r8 + &  
                          state%mf_lat_n(i  ,j-1) * mesh%half_cos_lat(j-1) * (state%pv(i-1,j-1) + state%pv(i,j-1) + state%pv(i,j  )) / 3.0_r8 + &
                          state%mf_lat_n(i+1,j-1) * mesh%half_cos_lat(j-1) * (state%pv(i+1,j-1) + state%pv(i,j-1) + state%pv(i,j  )) / 3.0_r8)) * 0.25_r8 / &
#endif
                        mesh%full_cos_lat(j)
      end do
    end do

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
    	do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef V_POLE
        tend%qhu(i,j) = (state%mf_lon_n(i-1,j-1) * (state%pv(i-1,j-1) + state%pv(i-1,j) + state%pv(i,j)) / 3.0_r8 +&
                         state%mf_lon_n(i  ,j-1) * (state%pv(i  ,j-1) + state%pv(i-1,j) + state%pv(i,j)) / 3.0_r8 +&
                         state%mf_lon_n(i-1,j  ) * (state%pv(i-1,j+1) + state%pv(i-1,j) + state%pv(i,j)) / 3.0_r8 +&
                         state%mf_lon_n(i  ,j  ) * (state%pv(i  ,j+1) + state%pv(i-1,j) + state%pv(i,j)) / 3.0_r8) * 0.25_r8
#else
        tend%qhu(i,j) = (state%mf_lon_n(i-1,j+1) * (state%pv(i-1,j+1) + state%pv(i-1,j) + state%pv(i,j)) / 3.0_r8 +&
                         state%mf_lon_n(i  ,j+1) * (state%pv(i  ,j+1) + state%pv(i-1,j) + state%pv(i,j)) / 3.0_r8 +&
                         state%mf_lon_n(i-1,j  ) * (state%pv(i-1,j-1) + state%pv(i-1,j) + state%pv(i,j)) / 3.0_r8 +&
                         state%mf_lon_n(i  ,j  ) * (state%pv(i  ,j-1) + state%pv(i-1,j) + state%pv(i,j)) / 3.0_r8) * 0.25_r8
#endif
      end do
    end do

  end subroutine calc_qhu_qhv_2

  subroutine calc_dkedlon_dkedlat(state, tend, dt)
    
    type(state_type), intent(inout) :: state
    type(tend_type ), intent(inout) :: tend
    real(r8)        , intent(in   ) :: dt 

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        tend%dkedlon(i,j) = (state%ke(i+1,j) - state%ke(i,j)) / (radius * mesh%full_cos_lat(j) * mesh%dlon)
      end do
    end do

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef V_POLE
        tend%dkedlat(i,j) = (state%ke(i,j) - state%ke(i,j-1)) / (radius * mesh%dlat)
#else
        tend%dkedlat(i,j) = (state%ke(i,j+1) - state%ke(i,j)) / (radius * mesh%dlat)
#endif
      end do
    end do     

  end subroutine calc_dkedlon_dkedlat

  subroutine calc_dpedlon_dpedlat(static, state, tend, dt)
    
    type(static_type), intent(in   ) :: static
    type(state_type ), intent(inout) :: state
    type(tend_type  ), intent(inout) :: tend
    real(r8)         , intent(in   ) :: dt 

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        tend%dpedlon(i,j) = (state%gd(i+1,j) -   state%gd(i,j) +&
                           static%ghs(i+1,j) - static%ghs(i,j)) /&
                              (radius * mesh%full_cos_lat(j) * mesh%dlon)
      end do
    end do

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef V_POLE
        tend%dpedlat(i,j) = (state%gd(i,j) -   state%gd(i,j-1) +&
                           static%ghs(i,j) - static%ghs(i,j-1)) /&
                            (radius * mesh%dlat)
#else
        tend%dpedlat(i,j) = (state%gd(i,j+1) -   state%gd(i,j) +&
                           static%ghs(i,j+1) - static%ghs(i,j)) /&
                            (radius * mesh%dlat)
#endif
      end do
    end do     

  end subroutine calc_dpedlon_dpedlat

  subroutine calc_dmfdlon_dmfdlat(state, tend, dt)
  
    type(state_type), intent(inout) :: state
    type(tend_type ), intent(inout) :: tend
    real(r8)        , intent(in   ) :: dt 

    type(mesh_type), pointer :: mesh
    integer i, j
    real(r8) pole

    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        tend%dmfdlon(i,j) = (state%mf_lon_n(i,j) - state%mf_lon_n(i-1,j)) /&
                           (radius * mesh%full_cos_lat(j) * mesh%dlon)
      end do
    end do

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef V_POLE
        tend%dmfdlat(i,j) = (state%mf_lat_n(i,j+1) * mesh%half_cos_lat(j+1) - &
                             state%mf_lat_n(i,j  ) * mesh%half_cos_lat(j  )) /&
                             (radius * mesh%full_cos_lat(j) * mesh%dlat)
#else
        tend%dmfdlat(i,j) = (state%mf_lat_n(i,j  ) * mesh%half_cos_lat(j  ) - &
                             state%mf_lat_n(i,j-1) * mesh%half_cos_lat(j-1)) /&
                             (radius * mesh%full_cos_lat(j) * mesh%dlat)
#endif
      end do
    end do
#ifndef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_start_idx
      pole = 0.0_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole + state%mf_lat_n(i,j)
      end do
      pole = pole * 4.0_r8 / mesh%num_full_lon / (radius * mesh%dlat)
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        tend%dmfdlat(i,j) = pole
      end do
    end if

    if (mesh%has_north_pole()) then
      j = mesh%full_lat_end_idx
      pole = 0.0_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole + state%mf_lat_n(i,j-1)
      end do
      pole = -pole * 4.0_r8 / mesh%num_full_lon / (radius * mesh%dlat)
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        tend%dmfdlat(i,j) = pole
      end do
    end if
#endif
  end subroutine calc_dmfdlon_dmfdlat

subroutine test03()

  INTEGER,parameter :: N=5
  integer:: LDA, LDB
  INTEGER :: NRHS
  INTEGER :: INFO
  INTEGER :: IPIV(N)
  REAL(r8) :: A(N,N), B(N)
  
  integer i,j

  LDA=N;LDB=N
  NRHS=1
  
  ! A=reshape((/2.0,4.0,&
  !             3.0,-1.5/),(/2,2/))
  ! A = reshape([2.0, 4.0 , 1.0,&
  !              3.0, -1.5, 4.0,&
  !              4.0, -2.0, -5.0], (/3,3/))
  ! B = reshape((/20.0,-5.0,-6.0/),(/3,1/))

  do i = 1, n 
    do j = 1, n
      if( i == j) then 
        a(i,j) = -1.0
    	else if(i-j == 1) then 
        a(i,j) = 1.0
      else
        a(i,j) = 0.0
      end if 
    end do
  end do 
  a(1,n)=1
  do j=1,n
  write(*,100)  (a(i,j),i=1,n)
  end do

100 format(5f10.2)

  b(1) = 1.0
  b(2:n-1) = 0.0
  b(n) = 1.0

  call DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
  
  write(*,*) "Solution:"
  write(*,'(f8.3)') B
  write(*,*) "INFO=", INFO

end subroutine test03

end module operators_mod