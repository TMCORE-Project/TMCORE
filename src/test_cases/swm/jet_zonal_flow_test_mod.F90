module jet_zonal_flow_test_mod

  use params_mod, r8 => real_kind
  use mesh_mod
  use state_mod
  use static_mod

  implicit none

  private

  public jet_zonal_flow_test_set_initial_condition

  real(r8), parameter :: u_max = 80.0_r8
  real(r8), parameter :: lat0 = pi / 7.0_r8
  real(r8), parameter :: lat1 = pi / 2.0_r8 - lat0
  real(r8), parameter :: en = exp(-4.0_r8 / (lat1 - lat0)**2_r8)
  real(r8), parameter :: gh0 = g * 1.0e4_r8
  real(r8), parameter :: ghd = g * 120_r8
  real(r8), parameter :: lat2 = pi / 4.0_r8
  real(r8), parameter :: alpha = 1.0_r8 / 3.0_r8
  real(r8), parameter :: beta = 1.0_r8 / 15.0_r8

contains

  subroutine jet_zonal_flow_test_set_initial_condition()

    real(r8), allocatable :: psiVertex(:)
    integer iCell, iEdge, iVertex
    integer neval, ierr
    real(r8) abserr

    call log_notice('Use jet zonal flow initial condition.')

    static%cell%ghs = 0.0_r8

    allocate(psiVertex(nVertices))
    do iVertex = lbound(state(1)%vertex%pv, 1), ubound(state(1)%vertex%pv, 1)
      if (latVertex(iVertex) <= lat0 .or. latVertex(iVertex) >= lat1) then
        psiVertex(iVertex) = 0.0_r8
      else
        psiVertex(iVertex) = u_max / en * exp(1.0_r8 / (latVertex(iVertex) - lat0) / (latVertex(iVertex) - lat1)) * &
          (- 1.0_r8 / (latVertex(iVertex) - lat0)    / (latVertex(iVertex) - lat1)**2 &
           - 1.0_r8 / (latVertex(iVertex) - lat0)**2 / (latVertex(iVertex) - lat1))
      end if
    end do

    do iEdge = lbound(state(1)%edge%u, 1), ubound(state(1)%edge%u, 1)
      state(1)%edge%u(iEdge) = -( psiVertex(verticesOnEdge(2,iEdge)) &
                                - psiVertex(verticesOnEdge(1,iEdge)) ) / dvEdge(iEdge)
    end do
    deallocate(psiVertex)

    do iCell = lbound(static%cell%ghs, 1), ubound(static%cell%ghs, 1)
      call qags(gh_integrand, -0.5*pi, latCell(iCell), 1.0e-10, 1.0e-3, state(1)%cell%gd(iCell), abserr, neval, ierr)
      state(1)%cell%gd(iCell) = gh0 - state(1)%cell%gd(iCell)
      state(1)%cell%gd(iCell) = state(1)%cell%gd(iCell) + ghd * cos(latCell(iCell)) * &
        exp(-((lonCell(iCell) - pi) / alpha)**2) * exp(-((lat2 - latCell(iCell)) / beta)**2)
    end do

  end subroutine jet_zonal_flow_test_set_initial_condition

  real function gh_integrand(lat) result(res)

    real, intent(in) :: lat

    real u, f

    u = u_function(lat)
    f = 2 * omega * sin(lat)
    res = radius * u * (f + tan(lat) / radius * u)

  end function gh_integrand

  real function u_function(lat) result(res)

    real, intent(in) :: lat

    if (lat <= lat0 .or. lat >= lat1) then
      res = 0.0
    else
      res = u_max / en * exp(1 / (lat - lat0) / (lat - lat1))
    end if

  end function u_function

end module jet_zonal_flow_test_mod
