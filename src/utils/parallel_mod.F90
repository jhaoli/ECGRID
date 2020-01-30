module parallel_mod

  use mesh_mod

  implicit none

  private
  
  public parallel_fill_halo
  public parallel_final

  interface parallel_fill_halo
    module procedure parallel_fill_halo_1d_r8_1
    module procedure parallel_fill_halo_1d_r8_2
    module procedure parallel_fill_halo_2d_r8
  end interface parallel_fill_halo


contains

  subroutine parallel_final()

    integer ierr

    ! call mpi_finalize(ierr)

  end subroutine parallel_final
  
  subroutine parallel_fill_halo_1d_r8_1(halo_width, array, left_halo, right_halo)

    integer, intent(in   )           :: halo_width
    real(8), intent(inout)           :: array(:)
    logical, intent(in   ), optional :: left_halo
    logical, intent(in   ), optional :: right_halo

    integer i, m, n

    if (merge(left_halo, .true., present(left_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * halo_width
      do i = 1, halo_width
        array(m+i) = array(n+i)
      end do
    end if

    if (merge(right_halo, .true., present(right_halo))) then
      m = ubound(array, 1) - halo_width
      n = lbound(array, 1) + halo_width - 1
      do i = 1, halo_width
        array(m+i) = array(n+i)
      end do
    end if

  end subroutine parallel_fill_halo_1d_r8_1

  subroutine parallel_fill_halo_1d_r8_2(mesh, array, left_halo, right_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(:)
    logical        , intent(in   ), optional :: left_halo
    logical        , intent(in   ), optional :: right_halo

    integer i, m, n

    if (merge(left_halo, .true., present(left_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * mesh%halo_width
      do i = 1, mesh%halo_width
        array(m+i) = array(n+i)
      end do
    end if

    if (merge(right_halo, .true., present(right_halo))) then
      m = ubound(array, 1) - mesh%halo_width
      n = lbound(array, 1) + mesh%halo_width - 1
      do i = 1, mesh%halo_width
        array(m+i) = array(n+i)
      end do
    end if

  end subroutine parallel_fill_halo_1d_r8_2

  subroutine parallel_fill_halo_2d_r8(mesh, array, left_halo, right_halo, top_halo, bottom_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(:,:)
    logical        , intent(in   ), optional :: left_halo
    logical        , intent(in   ), optional :: right_halo
    logical        , intent(in   ), optional :: top_halo
    logical        , intent(in   ), optional :: bottom_halo

    integer i, j, m, n

    if (merge(left_halo, .true., present(left_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * mesh%halo_width
      do j = lbound(array, 2), ubound(array, 2)
        do i = 1, mesh%halo_width
          array(m+i,j) = array(n+i,j)
        end do
      end do
    end if

    if (merge(right_halo, .true., present(right_halo))) then
      m = ubound(array, 1) - mesh%halo_width
      n = lbound(array, 1) + mesh%halo_width - 1
      do j = lbound(array, 2), ubound(array, 2)
        do i = 1, mesh%halo_width
          array(m+i,j) = array(n+i,j)
        end do
      end do
    end if

  end subroutine parallel_fill_halo_2d_r8

end module parallel_mod