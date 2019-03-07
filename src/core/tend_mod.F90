module tend_mod

  use const_mod
  use mesh_mod

  implicit none

  private

  public tend_init
  public tend_final
  public tend_type
  public tend

  type tend_on_cell_type
    real(real_kind), allocatable :: gd(:)
  end type tend_on_cell_type

  type tend_on_edge_type
    real(real_kind), allocatable :: gd   (:)
    real(real_kind), allocatable :: u    (:)
    real(real_kind), allocatable :: iap_u(:)
  end type tend_on_edge_type

  type tend_on_vertex_type
    real(real_kind), allocatable :: gd(:)
  end type tend_on_vertex_type

  type tend_type
    type(tend_on_cell_type)   cell
    type(tend_on_edge_type)   edge
    type(tend_on_vertex_type) vertex
  end type tend_type

  type(tend_type), allocatable :: tend(:)

contains

  subroutine tend_init()

    integer i

    if (.not. allocated(tend)) then
      allocate(tend(0:2))
      do i = lbound(tend, 1), ubound(tend, 1)
        call allocate_tend(tend(i))
      end do
    end if

  end subroutine tend_init

  subroutine tend_final()

    integer i

    if (allocated(tend)) then
      do i = lbound(tend, 1), ubound(tend, 1)
        call deallocate_tend(tend(i))
      end do
    end if

  end subroutine tend_final

  subroutine allocate_tend(tend)

    type(tend_type), intent(out) :: tend

    if (.not. allocated(tend%cell%gd)) allocate(tend%cell%gd(nCells))

    if (.not. allocated(tend%edge%gd   )) allocate(tend%edge%gd   (nEdges))
    if (.not. allocated(tend%edge%u    )) allocate(tend%edge%u    (nEdges))
    if (.not. allocated(tend%edge%iap_u)) allocate(tend%edge%iap_u(nEdges))

    if (.not. allocated(tend%vertex%gd)) allocate(tend%vertex%gd(nVertices))

  end subroutine allocate_tend

  subroutine deallocate_tend(tend)

    type(tend_type), intent(inout) :: tend

    if (allocated(tend%cell%gd)) deallocate(tend%cell%gd)

    if (allocated(tend%edge%gd   )) deallocate(tend%edge%gd   )
    if (allocated(tend%edge%u    )) deallocate(tend%edge%u    )
    if (allocated(tend%edge%iap_u)) deallocate(tend%edge%iap_u)

    if (allocated(tend%vertex%gd)) deallocate(tend%vertex%gd)

  end subroutine deallocate_tend

end module tend_mod