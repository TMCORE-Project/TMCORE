module state_mod

  use const_mod
  use mesh_mod

  implicit none

  private

  public state_init
  public state_final
  public state_type
  public state

  type state_on_cell_type
    real(real_kind), allocatable :: gd    (:)
    real(real_kind), allocatable :: u     (:)
    real(real_kind), allocatable :: v     (:)
    real(real_kind), allocatable :: vor   (:)
    real(real_kind), allocatable :: div   (:)
    real(real_kind), allocatable :: ke    (:)
    real(real_kind), allocatable :: pv    (:)
    real(real_kind), allocatable :: iap_gh(:)
  end type state_on_cell_type

  type state_on_edge_type
    real(real_kind), allocatable :: gd    (:)
    real(real_kind), allocatable :: u     (:) ! Normal  velocity
    real(real_kind), allocatable :: v     (:) ! Tangent velocity
    real(real_kind), allocatable :: pv    (:)
    real(real_kind), allocatable :: pv_flx(:)
    real(real_kind), allocatable :: iap_u (:)
  end type state_on_edge_type

  type state_on_vertex_type
    real(real_kind), allocatable :: gd (:)
    real(real_kind), allocatable :: vor(:)
    real(real_kind), allocatable :: pv (:)
  end type state_on_vertex_type

  type state_type
    type(state_on_cell_type)   cell
    type(state_on_edge_type)   edge
    type(state_on_vertex_type) vertex
    real(real_kind) total_mass
    real(real_kind) total_energy
    real(real_kind) total_potential_enstropy
    real(real_kind) total_angular_momentum
    real(real_kind) total_absolute_vorticity
  end type state_type

  type(state_type), allocatable :: state(:)

contains

  subroutine state_init()

    integer i

    if (.not. allocated(state)) then
      allocate(state(-4:2))
      do i = lbound(state, 1), ubound(state, 1)
        call allocate_state(state(i))
      end do
    end if

  end subroutine state_init

  subroutine state_final()

    integer i

    if (allocated(state)) then
      do i = lbound(state, 1), ubound(state, 1)
        call deallocate_state(state(i))
      end do
      deallocate(state)
    end if

  end subroutine state_final

  subroutine allocate_state(state)

    type(state_type), intent(out) :: state

    ! Only consider serial case, leave parallel case for later work.
    if (.not. allocated(state%cell%gd    )) allocate(state%cell%gd    (nCells))
    if (.not. allocated(state%cell%u     )) allocate(state%cell%u     (nCells))
    if (.not. allocated(state%cell%v     )) allocate(state%cell%v     (nCells))
    if (.not. allocated(state%cell%vor   )) allocate(state%cell%vor   (nCells))
    if (.not. allocated(state%cell%div   )) allocate(state%cell%div   (nCells))
    if (.not. allocated(state%cell%ke    )) allocate(state%cell%ke    (nCells))
    if (.not. allocated(state%cell%pv    )) allocate(state%cell%pv    (nCells))

    if (.not. allocated(state%edge%gd    )) allocate(state%edge%gd    (nEdges))
    if (.not. allocated(state%edge%u     )) allocate(state%edge%u     (nEdges))
    if (.not. allocated(state%edge%v     )) allocate(state%edge%v     (nEdges))
    if (.not. allocated(state%edge%pv    )) allocate(state%edge%pv    (nEdges))
    if (.not. allocated(state%edge%pv_flx)) allocate(state%edge%pv_flx(nEdges))
    if (.not. allocated(state%edge%iap_u )) allocate(state%edge%iap_u (nEdges))

    if (.not. allocated(state%vertex%gd  )) allocate(state%vertex%gd  (nVertices))
    if (.not. allocated(state%vertex%vor )) allocate(state%vertex%vor (nVertices))
    if (.not. allocated(state%vertex%pv  )) allocate(state%vertex%pv  (nVertices))

  end subroutine allocate_state

  subroutine deallocate_state(state)

    type(state_type), intent(inout) :: state

    if (allocated(state%cell%gd    )) deallocate(state%cell%gd    )
    if (allocated(state%cell%u     )) deallocate(state%cell%u     )
    if (allocated(state%cell%v     )) deallocate(state%cell%v     )
    if (allocated(state%cell%vor   )) deallocate(state%cell%vor   )
    if (allocated(state%cell%div   )) deallocate(state%cell%div   )
    if (allocated(state%cell%ke    )) deallocate(state%cell%ke    )
    if (allocated(state%cell%pv    )) deallocate(state%cell%pv    )

    if (allocated(state%edge%gd    )) deallocate(state%edge%gd    )
    if (allocated(state%edge%u     )) deallocate(state%edge%u     )
    if (allocated(state%edge%v     )) deallocate(state%edge%v     )
    if (allocated(state%edge%pv    )) deallocate(state%edge%pv    )
    if (allocated(state%edge%pv_flx)) deallocate(state%edge%pv_flx)
    if (allocated(state%edge%iap_u )) deallocate(state%edge%iap_u )

    if (allocated(state%vertex%gd  )) deallocate(state%vertex%gd  )
    if (allocated(state%vertex%vor )) deallocate(state%vertex%vor )
    if (allocated(state%vertex%pv  )) deallocate(state%vertex%pv  )

  end subroutine deallocate_state

end module state_mod