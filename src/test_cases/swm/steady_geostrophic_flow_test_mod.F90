module steady_geostrophic_flow_test_mod

  use params_mod
  use mesh_mod
  use static_mod
  use state_mod

  implicit none

  private

  public steady_geostrophic_flow_test_set_initial_condition

  real(real_kind), parameter :: alpha = 0.d0
  real(real_kind), parameter :: u0 = 2 * pi * radius / (12 * 86400)
  real(real_kind), parameter :: gd0 = 2.94e4 ! m2 s-2

contains

  subroutine steady_geostrophic_flow_test_set_initial_condition()

    real (kind=real_kind), allocatable :: psiVertex(:)

    integer                            :: iCell, iEdge, iVertex
    
    write(6, *) '[Notice]: Use steady geostrophic flow initial condition.'
    
    !
    ! Initialize wind field
    !
    allocate(psiVertex(nVertices))
    do iVertex = lbound(state(1)%vertex%pv, 1), ubound(state(1)%vertex%pv, 1)
       psiVertex(iVertex) = -radius * u0 * ( dsin(latVertex(iVertex)) * dcos(alpha)                         &
                                           - dcos(lonVertex(iVertex)) * dcos(latVertex(iVertex)) * dsin(alpha))
    end do
    
    do iEdge = lbound(state(1)%edge%u, 1), ubound(state(1)%edge%u, 1)
       state(1)%edge%u(iEdge) = -( psiVertex(verticesOnEdge(2,iEdge)) &
                                 - psiVertex(verticesOnEdge(1,iEdge)) ) / dvEdge(iEdge)
    end do
    deallocate(psiVertex)

    !
    ! Initialize height field (actually, fluid thickness field)
    !
    do iCell = lbound(static%cell%ghs, 1), ubound(static%cell%ghs, 1)
       state(1)%cell%gd(iCell) = gd0 - (radius * omega * u0 + 0.5d0 * u0**2.d0) * &
                                       (-dcos(lonCell(iCell)) * dcos(latCell(iCell)) * dsin(alpha) &
                                        +dsin(latCell(iCell)) * dcos(alpha))**2.d0
    end do
    
    static%cell%ghs = 0.d0

  end subroutine steady_geostrophic_flow_test_set_initial_condition

end module steady_geostrophic_flow_test_mod
