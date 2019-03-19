module rossby_haurwitz_wave_test_mod

  use params_mod
  use mesh_mod
  use static_mod
  use state_mod

  implicit none

  private

  public rossby_haurwitz_wave_test_set_initial_condition

  real (kind=real_kind), parameter :: h0 = 8000.d0
  real (kind=real_kind), parameter :: w = 7.848e-6
  real (kind=real_kind), parameter :: K = 7.848e-6
  real (kind=real_kind), parameter :: R = 4.d0

contains

  subroutine rossby_haurwitz_wave_test_set_initial_condition()
    implicit none
    
    integer :: iCell, iEdge, iVertex
    real (kind=real_kind), dimension(nVertices) :: psiVertex
    real (kind=real_kind), dimension(nCells   ) :: A, B, C
    
    write(6, *) '[Notice]: Use rossby haurwitz wave initial condition.'
    
    !
    ! Initialize wind field
    !
    do iVertex = lbound(state(1)%vertex%pv, 1), ubound(state(1)%vertex%pv, 1)
      psiVertex(iVertex) = - radius * radius * w *   dsin(latVertex(iVertex))       &
                           + radius * radius * K * ( dcos(latVertex(iVertex))**R)   &
                                                   * dsin(latVertex(iVertex))       &
                                                   * dcos(R * lonVertex(iVertex))
    end do
    
    do iEdge = lbound(state(1)%edge%u, 1), ubound(state(1)%edge%u, 1)
      state(1)%edge%u(iEdge) = -( psiVertex(verticesOnEdge(2,iEdge)) &
                                - psiVertex(verticesOnEdge(1,iEdge)) &
                                ) / dvEdge(iEdge)
    end do

    A = 0.5d0 * w * (2.d0 * omega + w) * cos(latCell)**2.d0 + &
         0.25d0 * K**2.d0 * dcos(latCell)**(2.d0*R) * ((R+1.d0)*cos(latCell)**2.d0 + 2.d0*R**2.d0 - R - 2.d0 - 2.d0*R**2.d0 * cos(latCell)**(-2.d0))
    B = (2.0*(omega + w)*K / ((R+1.0)*(R+2.0))) * cos(latCell)**R * ((R**2.0 + 2.0*R + 2.0) - ((R+1.0)*cos(latCell))**2.0)
    C = 0.25 * K**2.0 * cos(latCell)**(2.0*R) * ((R+1.0)*dcos(latCell)**2.0 - R - 2.0)
    
    !
    ! Initialize height field (actually, fluid thickness field)
    !
    do iCell = lbound(state(1)%cell%gd, 1), ubound(state(1)%cell%gd, 1)
      state(1)%cell%gd(iCell) = ( g * h0                                                     &
                                + radius * radius * A(iCell)                                 &
                                + radius * radius * B(iCell) * dcos(    R * lonCell(iCell))  &
                                + radius * radius * C(iCell) * dcos(2.0*R * lonCell(iCell)))
    end do
    
    static%cell%ghs = 0.d0
    !
    ! Initialize vorticity on the cells.
    !
    state(1)%cell%pv  = ( fCell + 2 * w * sin(latCell) - K * sin(latCell) * cos(latCell)**R * (R**2 + 3 * R + 2) * cos(R * lonCell) ) / state(1)%cell%gd

  end subroutine rossby_haurwitz_wave_test_set_initial_condition

end module rossby_haurwitz_wave_test_mod