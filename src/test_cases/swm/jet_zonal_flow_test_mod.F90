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

    integer iCell, iEdge, iVertex
    integer neval, ierr
    real(r8) abserr
    real(r8) u, v
    real(r8) angle1,angle2
    real(r8) lon

    call log_notice('Use jet zonal flow initial condition.')

    static%cell%ghs = 0.0_r8

    v = 0._r8
    do iEdge = 1,nEdges
      u      = u_function(latEdge(iEdge))
      state(1)%edge%u(iEdge) = cos(angleEdge(iEdge))*sqrt(u**2+v**2)
    enddo

    do iCell = lbound(static%cell%ghs, 1), ubound(static%cell%ghs, 1)
      call qags(gh_integrand, -0.5_r8*pi, latCell(iCell), 1.0e-10, 1.0e-3, state(1)%cell%gd(iCell), abserr, neval, ierr)
      
      state(1)%cell%gd(iCell) = gh0 - state(1)%cell%gd(iCell)
      
      if(lonCell(iCell)>pi)then
        lon = lonCell(iCell) - 2._r8 * pi
      else
        lon = lonCell(iCell)
      endif
      
      state(1)%cell%gd(iCell) = state(1)%cell%gd(iCell) + ghd * cos(latCell(iCell)) * exp(-(lon / alpha)**2) * exp(-((lat2 - latCell(iCell)) / beta)**2)
    end do

  end subroutine jet_zonal_flow_test_set_initial_condition

  real function gh_integrand(lat) result(res)

    real(r8), intent(in) :: lat

    real(r8) u, f

    u = u_function(lat)
    f = 2 * omega * sin(lat)
    res = radius * u * (f + tan(lat) / radius * u)

  end function gh_integrand

  real function u_function(lat) result(res)

    real(r8), intent(in) :: lat

    if (lat <= lat0 .or. lat >= lat1) then
      res = 0._r8
    else
      res = u_max / en * exp(1 / (lat - lat0) / (lat - lat1))
    end if

  end function u_function

end module jet_zonal_flow_test_mod
