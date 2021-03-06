module static_mod
  
  use flogger
  use const_mod
  use mesh_mod
  use allocator_mod

  implicit none

  private
  public static_type
  public static
  public static_init_root

  type static_type
    type(mesh_type), pointer :: mesh => null()
    real(r8), allocatable :: ghs(:,:)
  contains
    procedure :: init => static_init
    procedure :: clear => static_clear
    final :: static_final
  end type static_type

  type(static_type), target :: static

contains
  
  subroutine static_init_root()
    
    call static%init(mesh)
    call log_notice('static module is initialized.')

  end subroutine static_init_root
  
  subroutine static_init(this, mesh)
    
    class(static_type), intent(inout)         :: this 
    type(mesh_type   ), intent(in   ), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%ghs, full_lon=.true., full_lat=.true.)

  end subroutine static_init

  subroutine static_clear(this)
    
    class(static_type), intent(inout) :: this

    if (allocated(this%ghs)) deallocate(this%ghs)

  end subroutine static_clear

  subroutine static_final(this)
    
    type(static_type), intent(inout) :: this

    call this%clear()

  end subroutine static_final

end module static_mod