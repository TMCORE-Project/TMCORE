module cosine_bell_rotation_test_mod

  use params_mod
  use mesh_mod
  use static_mod
  use state_mod

  implicit none

  private

  public cosine_bell_rotation_test_set_initial_condition
  public cosine_bell_rotation_test_update_wind

  real(real_kind), parameter :: h0 = 1000.0d0 ! m
  real(real_kind), parameter :: u0 = 2.0d0 * pi * radius / 12.0d0 / 86400.0d0
  real(real_kind), parameter :: lon0 = 1.5d0 * pi
  real(real_kind), parameter :: lat0 = 0.0d0
  real(real_kind), parameter :: R0 = radius / 3.0d0
  real(real_kind), parameter :: alpha = 0.0d0

contains

  subroutine cosine_bell_rotation_test_set_initial_condition()

    real(real_kind) sin_lat0, cos_lat0, r
    integer iCell

    sin_lat0 = sin(lat0)
    cos_lat0 = cos(lat0)

    do iCell = lbound(state(1)%cell%gd, 1), ubound(state(1)%cell%gd, 1)
      r = radius * acos(sin_lat0 * sin(latCell(iCell)) + cos_lat0 * cos(latCell(iCell)) * cos(lonCell(iCell) - lon0))
      if (r < R0) then
        state(1)%cell%gd(iCell) = 0.5d0 * R0 * (1.0d0 + cos(pi * r / R0))
      else
        state(1)%cell%gd(iCell) = 0.0d0
      end if
    end do

    call cosine_bell_rotation_test_update_wind(0.0d0, state(1))

  end subroutine cosine_bell_rotation_test_set_initial_condition

  subroutine cosine_bell_rotation_test_update_wind(seconds, state)

    real(real_kind),  intent(in)    :: seconds
    type(state_type), intent(inout) :: state

    real(real_kind) cos_lat, sin_lat, cos_lon, cos_alpha, sin_alpha, dlon, d
    real(real_kind) psi_vertex(nVertices)
    integer iEdge, iVertex

    cos_alpha = cos(alpha)
    sin_alpha = sin(alpha)

    do iVertex = lbound(state%vertex%pv, 1), ubound(state%vertex%pv, 1)
      cos_lat = cos(latVertex(iVertex))
      sin_lat = sin(latVertex(iVertex))
      cos_lon = cos(lonVertex(iVertex))
      psi_vertex(iVertex) = - radius * u0 * (sin_lat * cos_alpha - cos_lon * cos_lat * sin_alpha)
    end do

    do iEdge = lbound(state%edge%u, 1), ubound(state%edge%u, 1)
      state%edge%u(iEdge) = - (psi_vertex(verticesOnEdge(2,iEdge)) - psi_vertex(verticesOnEdge(1,iEdge))) / dvEdge(iEdge)
    end do

  end subroutine cosine_bell_rotation_test_update_wind

end module cosine_bell_rotation_test_mod