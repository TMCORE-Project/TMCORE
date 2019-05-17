module deformation_flow_test_mod

  use params_mod
  use log_mod
  use string_mod
  use sphere_geometry_mod
  use mesh_mod
  use static_mod
  use state_mod

  implicit none

  private

  public deformation_flow_test_set_initial_condition
  public deformation_flow_test_update_wind

  real(real_kind), parameter :: T0 = 12.0d0 * 86400.0d0
  real(real_kind), parameter :: u0 = 2.0d0 * pi * radius / T0
  real(real_kind), parameter :: R0 = radius / 3.0d0
  real(real_kind), parameter :: alpha = 0.5d0 * pi

contains

  subroutine deformation_flow_test_set_initial_condition()

    real(real_kind) lon1, lat1, lon2, lat2, r1, r2, r, b, c
    real(real_kind) x1, y1, z1, x2, y2, z2, d1, d2
    integer iCell

    static%cell%ghs = 0.0d0

    lon1 = pi * 5.0d0 / 6.0d0
    lat1 = 0.0d0
    lon2 = pi * 7.0d0 / 6.0d0
    lat2 = 0.0d0
    call cartesian_transform(lon1, lat1, x1, y1, z1)
    call cartesian_transform(lon2, lat2, x2, y2, z2)

    select case (string_split(subcase_name, 2, ':'))
    case ('slotted_cylinders')
      r = 0.5d0 * radius
      b = 0.1d0 * g
      c = 1.0d0 * g
      do iCell = lbound(state(1)%cell%gd, 1), ubound(state(1)%cell%gd, 1)
        r1 = calc_distance(lon1, lat1, lonCell(iCell), latCell(iCell))
        r2 = calc_distance(lon2, lat2, lonCell(iCell), latCell(iCell))
        if ((r1 <= r .and. abs(lonCell(iCell) - lon1) >= r / 6.0d0 / radius) .or. (r2 <= r .and. abs(lonCell(iCell) - lon2) >= r / 6.0d0 / radius)) then
          state(1)%cell%gd(iCell) = c
        else if (r1 <= r .and. abs(lonCell(iCell) - lon1) < r / 6.0d0 / radius .and. latCell(iCell) - lat1 < -5.0d0 / 12.0d0 * r / radius) then
          state(1)%cell%gd(iCell) = c
        else if (r2 <= r .and. abs(lonCell(iCell) - lon2) < r / 6.0d0 / radius .and. latCell(iCell) - lat2 >  5.0d0 / 12.0d0 * r / radius) then
          state(1)%cell%gd(iCell) = c
        else
          state(1)%cell%gd(iCell) = b
        end if
      end do
    case ('cosine_bells')

    case ('gaussian_hills')
      do iCell = lbound(state(1)%cell%gd, 1), ubound(state(1)%cell%gd, 1)
        d1 = ((x1 - xCell(iCell))**2 + (y1 - yCell(iCell))**2 + (z1 - zCell(iCell))**2) / radius**2
        d2 = ((x2 - xCell(iCell))**2 + (y2 - yCell(iCell))**2 + (z2 - zCell(iCell))**2) / radius**2
        state(1)%cell%gd(iCell) = 0.95 * (exp(-5.0 * d1) + exp(-5.0 * d2))
      end do
    case default
      call log_error('Unknown subcase_name ' // trim(subcase_name) // '!')
    end select

    call deformation_flow_test_update_wind(0.0d0, state(1))

  end subroutine deformation_flow_test_set_initial_condition

  subroutine deformation_flow_test_update_wind(seconds, state)

    real(real_kind),  intent(in)    :: seconds
    type(state_type), intent(inout) :: state

    real(real_kind) c, k, cos_T
    real(real_kind) psi_vertex(nVertices)
    integer iEdge, iVertex

    select case (string_split(subcase_name, 1, ':'))
    case ('case1')

    case ('case2')

    case ('case3')

    case ('case4')
      c = pi * 2.0d0 * seconds / T0
      k = 10.0d0 * radius / T0
      cos_T = cos(pi * seconds / T0)
      do iVertex = lbound(state%vertex%pv, 1), ubound(state%vertex%pv, 1)
        psi_vertex(iVertex) = k * sin(lonVertex(iVertex) - c)**2 * cos(latVertex(iVertex))**2 * cos_T - u0 * sin(latVertex(iVertex))
      end do
      psi_vertex = psi_vertex * radius
      do iEdge = lbound(state%edge%u, 1), ubound(state%edge%u, 1)
        state%edge%u(iEdge) = - (psi_vertex(verticesOnEdge(2,iEdge)) - psi_vertex(verticesOnEdge(1,iEdge))) / dvEdge(iEdge)
      end do
    case default
      call log_error('Unknown subcase_name ' // trim(subcase_name) // '!')
    end select

  end subroutine deformation_flow_test_update_wind

end module deformation_flow_test_mod
