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

    real (kind=RKIND), parameter   :: u0    = 2.d0 * pi * a
    real (kind=RKIND), allocatable :: psiVertex(:)

    integer                        :: iCell, iEdge, iVtx
    
    !
    ! Initialize wind field
    !
    allocate(psiVertex(nVertices))
    do iVtx = 1, nVertices
       psiVertex(iVtx) = -a * u0 * ( dsin(latVertex(iVtx)) * dcos(alpha)                         &
                                    -dcos(lonVertex(iVtx)) * dcos(latVertex(iVtx)) * dsin(alpha))
    end do
    
    do iEdge = 1,nEdges
       u(iEdge) = -( psiVertex(verticesOnEdge(2,iEdge)) &
                    -psiVertex(verticesOnEdge(1,iEdge)) ) / dvEdge(iEdge)
    end do
    deallocate(psiVertex)

    !
    ! Initialize height field (actually, fluid thickness field)
    !
    do iCell = 1, nCells
       wh(iCell) = gh0 - (a * omega * u0 + 0.5d0 * u0**2.d0) * &
                   (-dcos(lonCell(iCell)) * dcos(latCell(iCell)) * dsin(alpha) &
                    +dsin(latCell(iCell)) * dcos(alpha))**2.d0
    end do
    
    wh_s = 0.d0

  end subroutine steady_geostrophic_flow_test_set_initial_condition

end module steady_geostrophic_flow_test_mod
