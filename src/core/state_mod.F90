module state_mod
 
  use flogger
  use const_mod
  use mesh_mod
  use namelist_mod
  use allocator_mod

  implicit none

  private

  public state_type
  public states
  public state_init_root

  type state_type
    type(mesh_type), pointer :: mesh => null()
    ! Prognostic variables
    real(r8), allocatable, dimension(:,:) :: u
    real(r8), allocatable, dimension(:,:) :: v
    real(r8), allocatable, dimension(:,:) :: gd
    ! Diagnostic variables
    real(r8), allocatable, dimension(:,:) :: pv
    real(r8), allocatable, dimension(:,:) :: m_vtx
    real(r8), allocatable, dimension(:,:) :: m_lon
    real(r8), allocatable, dimension(:,:) :: m_lat
    real(r8), allocatable, dimension(:,:) :: mf_lon_n
    real(r8), allocatable, dimension(:,:) :: mf_lat_n
    real(r8), allocatable, dimension(:,:) :: pv_lon
    real(r8), allocatable, dimension(:,:) :: pv_lat
    real(r8), allocatable, dimension(:,:) :: ke
    real(r8) vor_sp
    real(r8) vor_np
    real(r8) total_m
    real(r8) total_ke
    real(r8) total_e 
    real(r8) total_pe
    real(r8) total_pv
  contains
    procedure :: init => state_init
    procedure :: clear => state_clear
    final :: state_final
  end type state_type

  type(state_type), allocatable, target :: states(:)

contains
  
  subroutine state_init_root()
  
    integer i 

    if (.not. allocated(states)) then
      select case (trim(time_scheme))
      case ('pc2')
        allocate(states(3))
      case ('rk3')
        allocate(states(4))
      case ('rk4')
        allocate(states(5))
      case default
        call log_error('Wrong time scheme!')
      end select
      do i = lbound(states, 1), ubound(states, 1)
        call states(i)%init(mesh)
      end do
      call log_notice('state module is initialized')
    end if
  end subroutine state_init_root

  subroutine state_init(this, mesh)
  
    class(state_type), intent(inout)   :: this
    type(mesh_type  ), intent(in   ), target :: mesh 

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%u        , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%v        , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%gd       , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pv       , half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%m_vtx    , half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%m_lon    , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%m_lat    , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%mf_lon_n , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%mf_lat_n , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%pv_lon   , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pv_lat   , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%ke       , full_lon=.true., full_lat=.true.)

  end subroutine state_init

  subroutine state_clear(this)
    
    class(state_type), intent(inout) :: this

    if (allocated(this%u        ))  deallocate(this%u       )
    if (allocated(this%v        )) deallocate(this%v        )
    if (allocated(this%gd       )) deallocate(this%gd       )
    if (allocated(this%pv       )) deallocate(this%pv       )
    if (allocated(this%m_vtx    )) deallocate(this%m_vtx    )
    if (allocated(this%m_lon    )) deallocate(this%m_lon    )
    if (allocated(this%m_lat    )) deallocate(this%m_lat    )
    if (allocated(this%mf_lon_n )) deallocate(this%mf_lon_n )
    if (allocated(this%mf_lat_n )) deallocate(this%mf_lat_n )
    if (allocated(this%pv_lon   )) deallocate(this%pv_lon   )
    if (allocated(this%pv_lat   )) deallocate(this%pv_lat   )
    if (allocated(this%ke       )) deallocate(this%ke       )

  end subroutine state_clear

  subroutine state_final(this)
  
    type(state_type), intent(inout) :: this

    call this%clear()

  end subroutine state_final

end module state_mod