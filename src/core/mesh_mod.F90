module mesh_mod
  
  use flogger
  use const_mod
  use namelist_mod

  implicit none

  private

  public mesh_type
  public mesh

  type mesh_type
    integer halo_width
    integer num_full_lon
    integer num_half_lon
    integer num_full_lat
    integer num_half_lat
    integer full_lon_start_idx
    integer full_lon_end_idx
    integer full_lat_start_idx
    integer full_lat_end_idx
    integer full_lat_start_idx_no_pole
    integer full_lat_end_idx_no_pole
    integer half_lon_start_idx
    integer half_lon_end_idx
    integer half_lat_start_idx
    integer half_lat_end_idx
    integer half_lat_start_idx_no_pole
    integer half_lat_end_idx_no_pole
    integer full_lon_lb
    integer full_lon_ub
    integer half_lon_lb
    integer half_lon_ub
    integer full_lat_lb
    integer full_lat_ub
    integer half_lat_lb
    integer half_lat_ub
    real(r8) start_lon
    real(r8) end_lon
    real(r8) start_lat
    real(r8) end_lat
    real(r8) dlon
    real(r8) dlat
    real(r8), allocatable :: full_dlon(:)
    real(r8), allocatable :: full_dlat(:)
    real(r8), allocatable :: full_lon(:)
    real(r8), allocatable :: half_lon(:)
    real(r8), allocatable :: full_lat(:)
    real(r8), allocatable :: half_lat(:)
    real(r8), allocatable :: full_cos_lon(:)
    real(r8), allocatable :: half_cos_lon(:)
    real(r8), allocatable :: full_sin_lon(:)
    real(r8), allocatable :: half_sin_lon(:)
    real(r8), allocatable :: full_cos_lat(:)
    real(r8), allocatable :: half_cos_lat(:)
    real(r8), allocatable :: full_sin_lat(:)
    real(r8), allocatable :: half_sin_lat(:)
    ! For output
    real(r8), allocatable :: full_lon_deg(:)
    real(r8), allocatable :: half_lon_deg(:)
    real(r8), allocatable :: full_lat_deg(:)
    real(r8), allocatable :: half_lat_deg(:)
    ! Coriolis parameter
    real(r8), allocatable :: full_f(:)
    real(r8), allocatable :: half_f(:)
  contains
    procedure :: init => mesh_init
    procedure :: has_south_pole => mesh_has_south_pole
    procedure :: has_north_pole => mesh_has_north_pole
    ! procedure :: is_south_pole => mesh_is_south_pole
    ! procedure :: is_north_pole => mesh_is_north_pole
    final :: mesh_final
  end type mesh_type
  
  type(mesh_type), target :: mesh

contains
  subroutine mesh_init(this, num_lon, num_lat)
    
    class(mesh_type), intent(inout)  :: this
    integer         , intent(in   )  :: num_lon
    integer         , intent(in   )  :: num_lat

    integer i, j 

    this%num_full_lon = num_lon
    this%num_half_lon = num_lon
#ifdef V_POLE
    this%num_full_lat = num_lat - 1
    this%num_half_lat = num_lat
#else
    this%num_full_lat = num_lat
    this%num_half_lat = num_lat - 1
#endif
    this%full_lon_start_idx = 1
    this%full_lon_end_idx = this%num_full_lon
    this%full_lat_start_idx = 1
    this%full_lat_end_idx = this%num_full_lat
    this%half_lon_start_idx = 1
    this%half_lon_end_idx = this%num_half_lon
    this%half_lat_start_idx = 1
    this%half_lat_end_idx = this%num_half_lat

    this%halo_width = 2
    this%start_lon  = 0.0_r8
    this%end_lon    = pi2
    this%start_lat  = -pi05
    this%end_lat    = pi05
#ifdef V_POLE
    this%full_lat_start_idx_no_pole = this%full_lat_start_idx
    this%full_lat_end_idx_no_pole   = this%full_lat_end_idx
    this%half_lat_start_idx_no_pole = merge(this%half_lat_start_idx + 1, this%half_lat_start_idx, this%has_south_pole())
    this%half_lat_end_idx_no_pole   = merge(this%half_lat_end_idx   - 1, this%half_lat_end_idx  , this%has_north_pole())
#else
    this%full_lat_start_idx_no_pole = merge(this%full_lat_start_idx + 1, this%full_lat_start_idx, this%has_south_pole())
    this%full_lat_end_idx_no_pole   = merge(this%full_lat_end_idx   - 1, this%full_lat_end_idx  , this%has_north_pole())
    this%half_lat_start_idx_no_pole = this%half_lat_start_idx
    this%half_lat_end_idx_no_pole   = this%half_lat_end_idx
#endif
    this%full_lon_lb = this%full_lon_start_idx - this%halo_width
    this%full_lon_ub = this%full_lon_end_idx   + this%halo_width
    this%full_lat_lb = this%full_lat_start_idx - 1
    this%full_lat_ub = this%full_lat_end_idx   + 1
    this%half_lon_lb = this%half_lon_start_idx - this%halo_width
    this%half_lon_ub = this%half_lon_end_idx   + this%halo_width
    this%half_lat_lb = this%half_lat_start_idx - 1
    this%half_lat_ub = this%half_lat_end_idx   + 1

    allocate(this%full_lon      (this%full_lon_lb:this%full_lon_ub))
    allocate(this%half_lon      (this%half_lon_lb:this%half_lon_ub))
    allocate(this%full_lat      (this%full_lat_lb:this%full_lat_ub))
    allocate(this%half_lat      (this%half_lat_lb:this%half_lat_ub))
    allocate(this%full_cos_lon  (this%full_lon_lb:this%full_lon_ub))
    allocate(this%half_cos_lon  (this%half_lon_lb:this%half_lon_ub))
    allocate(this%full_sin_lon  (this%full_lon_lb:this%full_lon_ub))
    allocate(this%half_sin_lon  (this%half_lon_lb:this%half_lon_ub))
    allocate(this%full_cos_lat  (this%full_lat_lb:this%full_lat_ub))
    allocate(this%half_cos_lat  (this%half_lat_lb:this%half_lat_ub))
    allocate(this%full_sin_lat  (this%full_lat_lb:this%full_lat_ub))
    allocate(this%half_sin_lat  (this%half_lat_lb:this%half_lat_ub))
    allocate(this%full_lon_deg  (this%full_lon_lb:this%full_lon_ub))
    allocate(this%half_lon_deg  (this%half_lon_lb:this%half_lon_ub))
    allocate(this%full_lat_deg  (this%full_lat_lb:this%full_lat_ub))
    allocate(this%half_lat_deg  (this%half_lat_lb:this%half_lat_ub))
    allocate(this%full_f        (this%full_lat_lb:this%full_lat_ub))
    allocate(this%half_f        (this%half_lat_lb:this%half_lat_ub))
    allocate(this%full_dlon     (this%full_lat_lb:this%full_lat_ub))
    allocate(this%full_dlat     (this%full_lat_lb:this%full_lat_ub))

    this%dlon = (this%end_lon - this%start_lon) / this%num_full_lon
    do i = this%full_lon_lb, this%full_lon_ub
      this%full_lon(i) = this%start_lon + (i - 1) * this%dlon
      this%half_lon(i) = this%full_lon(i) + 0.5_r8 * this%dlon
      this%full_lon_deg(i) = this%full_lon(i) * deg
      this%half_lon_deg(i) = this%half_lon(i) * deg
    end do

#ifdef V_POLE
    this%dlat = (this%end_lat - this%start_lat) / this%num_full_lat
    do j = this%half_lat_lb, this%half_lat_ub
      this%half_lat(j) = this%start_lat + (j - 1) * this%dlat
      if (abs(this%half_lat(j)) < 1.0e-14) this%half_lat(j) = 0.0_r8
      this%half_lat_deg(j) = this%half_lat(j) * deg
    end do
    this%half_lat(this%num_half_lat) = this%end_lat
    this%half_lat_deg(this%num_half_lat) = this%end_lat * deg

    do j = this%full_lat_lb, this%full_lat_ub
      this%full_lat(j) = this%half_lat(j) + 0.5_r8 * this%dlat
      this%full_lat_deg(j) = this%full_lat(j) * deg
    end do
#else
    this%dlat = (this%end_lat - this%start_lat) / this%num_half_lat
    do j = this%full_lat_lb, this%full_lat_ub
      this%full_lat(j) = this%start_lat + (j - 1) * this%dlat
      this%full_lat_deg(j) = this%full_lat(j) * deg
    end do
    this%full_lat(this%num_full_lat) = this%end_lat
    this%full_lat_deg(this%num_full_lat) = this%end_lat * deg

    do j = this%half_lat_lb, this%half_lat_ub
      this%half_lat(j) = this%full_lat(j) + 0.5_r8 * this%dlat
      this%half_lat_deg(j) = this%half_lat(j) * deg
    end do
#endif

    do i = this%full_lon_lb, this%full_lon_ub
      this%full_cos_lon(i) = cos(this%full_lon(i))
      this%full_sin_lon(i) = sin(this%full_lon(i))
    end do

    do i = this%half_lon_lb, this%half_lon_ub
      this%half_cos_lon(i) = cos(this%half_lon(i))
      this%half_sin_lon(i) = sin(this%half_lon(i))
    end do

    do j = this%full_lat_lb, this%full_lat_ub
      if (this%full_lat(j) >= -pi05 .and. this%full_lat(j) <= pi05) then
        this%full_cos_lat(j) = cos(this%full_lat(j))
        this%full_sin_lat(j) = sin(this%full_lat(j))
      end if
    end do

    do j = this%half_lat_lb, this%half_lat_ub
      if (this%half_lat(j) >= -pi05 .and. this%half_lat(j) <= pi05) then
        this%half_cos_lat(j) = cos(this%half_lat(j))
        this%half_sin_lat(j) = sin(this%half_lat(j))
      end if
    end do

#ifdef V_POLE
    if (this%has_south_pole()) then
      this%half_cos_lat(this%half_lat_start_idx) =  0.0_r8
      this%half_sin_lat(this%half_lat_start_idx) = -1.0_r8
    end if
    if (this%has_north_pole()) then
      this%half_cos_lat(this%half_lat_end_idx) = 0.0_r8
      this%half_sin_lat(this%half_lat_end_idx) = 1.0_r8
    end if
#else
    if (this%has_south_pole()) then
      this%full_cos_lat(this%full_lat_start_idx) = 0.0_r8
      this%full_sin_lat(this%full_lat_start_idx) = -1.0_r8
    end if
    
    if (this%has_north_pole()) then
      this%full_cos_lat(this%full_lat_end_idx) = 0.0_r8
      this%full_sin_lat(this%full_lat_end_idx) = 1.0_r8
    end if
#endif

    do j = this%full_lat_start_idx, this%full_lat_end_idx
      this%full_f(j) = 2.0_r8 * omega * this%full_sin_lat(j)
    end do
    do j = this%half_lat_start_idx, this%half_lat_end_idx
      this%half_f(j) = 2.0_r8 * omega * this%half_sin_lat(j)
    end do

    do j = this%full_lat_start_idx, this%full_lat_end_idx
      this%full_dlon(j) = radius * this%full_cos_lat(j) * this%dlon
    end do
    do j = this%full_lat_start_idx, this%full_lat_end_idx
      this%full_dlat(j) = radius * this%full_cos_lat(j) * this%dlat
    end do

#ifndef V_POLE    
    call reset_cos_lat_at_poles(mesh)
#endif
    call log_notice('Mesh module is initialized.')

  end subroutine mesh_init

  logical function mesh_has_south_pole(this) result(res)
    
    class(mesh_type), intent(in) :: this

    res = this%start_lat == -0.5_r8 * pi

  end function mesh_has_south_pole

  logical function mesh_has_north_pole(this) result(res)

    class(mesh_type), intent(in) :: this
    
    res = this%end_lat == 0.5_r8 * pi
      
  end function mesh_has_north_pole

  subroutine mesh_final(this)
   
    type(mesh_type), intent(inout) :: this

    if (allocated(this%full_lon           )) deallocate(this%full_lon     )
    if (allocated(this%full_lat           )) deallocate(this%full_lat     )
    if (allocated(this%half_lon           )) deallocate(this%half_lon     )
    if (allocated(this%half_lat           )) deallocate(this%half_lat     )
    if (allocated(this%full_cos_lon       )) deallocate(this%full_cos_lon )
    if (allocated(this%half_cos_lon       )) deallocate(this%half_cos_lon )
    if (allocated(this%full_sin_lon       )) deallocate(this%full_sin_lon )
    if (allocated(this%half_sin_lon       )) deallocate(this%half_sin_lon )
    if (allocated(this%full_cos_lat       )) deallocate(this%full_cos_lat )
    if (allocated(this%half_cos_lat       )) deallocate(this%half_cos_lat )
    if (allocated(this%full_sin_lat       )) deallocate(this%full_sin_lat )
    if (allocated(this%half_sin_lat       )) deallocate(this%half_sin_lat )
    if (allocated(this%full_lon_deg       )) deallocate(this%full_lon_deg )
    if (allocated(this%half_lon_deg       )) deallocate(this%half_lon_deg )
    if (allocated(this%full_lat_deg       )) deallocate(this%full_lat_deg )
    if (allocated(this%half_lat_deg       )) deallocate(this%half_lat_deg )
    if (allocated(this%full_f             )) deallocate(this%full_f       )
    if (allocated(this%half_f             )) deallocate(this%half_f       )
    if (allocated(this%full_dlon          )) deallocate(this%full_dlon    )
    if (allocated(this%full_dlat          )) deallocate(this%full_dlat    )
  end subroutine mesh_final
  
  subroutine reset_cos_lat_at_poles(this)
  
    type(mesh_type), intent(inout) :: this
    integer j 

    j = this%full_lat_start_idx
    this%full_cos_lat(j) = this%half_cos_lat(this%half_lat_start_idx_no_pole) * 0.25_r8

    j = this%full_lat_end_idx
    this%full_cos_lat(j) = this%half_cos_lat(this%half_lat_end_idx_no_pole) * 0.25_r8

  end subroutine reset_cos_lat_at_poles
end module mesh_mod