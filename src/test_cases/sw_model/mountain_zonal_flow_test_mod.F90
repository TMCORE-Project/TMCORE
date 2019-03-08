module mountain_zonal_flow_test_mod

  use params_mod
  use mesh_mod
  use static_mod
  use state_mod

  implicit none

  private

  public mountain_zonal_flow_test_set_initial_condition

  real(real_kind), parameter :: alpha = 0.0
  real(real_kind), parameter :: u0 = 20.0
  real(real_kind), parameter :: gd0 = 5960.0 * g ! m2 s-2
  real(real_kind), parameter :: lon0 = pi * 1.5
  real(real_kind), parameter :: lat0 = pi / 6.0
  real(real_kind), parameter :: ghs0 = 2000.0 * g
  real(real_kind), parameter :: R = pi / 9.0

contains

  subroutine mountain_zonal_flow_test_set_initial_condition()

    real(real_kind) cos_lat, sin_lat, cos_lon, sin_lon, cos_alpha, sin_alpha, dlon, d
    real(real_kind) psi_vertex(nVertices)
    integer iCell, iEdge, iVertex

    write(6, *) '[Notice]: Use mountain zonal flow initial condition.'

    cos_alpha = cos(alpha)
    sin_alpha = sin(alpha)

    do iCell = lbound(static%cell%ghs, 1), ubound(static%cell%ghs, 1)
      dlon = abs(lonCell(iCell) - lon0)
      dlon = min(dlon, 2 * pi - dlon)
      d = min(R, sqrt(dlon**2 + (latCell(iCell) - lat0)**2))
      static%cell%ghs(iCell) = ghs0 * (1.0d0 - d / R)
    end do

    do iVertex = lbound(state(1)%vertex%pv, 1), ubound(state(1)%vertex%pv, 1)
      cos_lat = cos(latVertex(iVertex))
      sin_lat = sin(latVertex(iVertex))
      cos_lon = cos(lonVertex(iVertex))
      sin_lon = sin(lonVertex(iVertex))
      psi_vertex(iVertex) = - radius * u0 * (sin_lat * cos_alpha - cos_lon * cos_lat * sin_alpha)
    end do

    do iEdge = lbound(state(1)%edge%u, 1), ubound(state(1)%edge%u, 1)
      state(1)%edge%u(iEdge) = - (psi_vertex(verticesOnEdge(2,iEdge)) - psi_vertex(verticesOnEdge(1,iEdge))) / dvEdge(iEdge)
    end do

    do iCell = lbound(state(1)%cell%gd, 1), ubound(state(1)%cell%gd, 1)
      state(1)%cell%gd(iCell) = gd0 - (radius * omega * u0 + 0.5d0 * u0**2) * ( &
        - cos(lonCell(iCell)) * cos(latCell(iCell)) * sin_alpha + &
          sin(latCell(iCell)) * cos_alpha &
      )**2 - static%cell%ghs(iCell)
    end do

  end subroutine mountain_zonal_flow_test_set_initial_condition

end module mountain_zonal_flow_test_mod