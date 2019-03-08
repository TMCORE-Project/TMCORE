module tmcore_common_mod

  use params_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use log_mod
  use mesh_mod
  use state_mod
  use tend_mod
  use operators_mod

  implicit none

  private

  public time_integrate
  public calc_kinetic_energy
  public calc_tangent_wind
  public calc_tangent_vor_flux
  public calc_pv_on_vertex
  public calc_pv_on_edge
  public spatial_operators_interface
  public update_state_interface

  interface
    subroutine spatial_operators_interface(state, tend)
      import state_type
      import tend_type
      type(state_type), intent(inout) :: state
      type(tend_type),  intent(inout) :: tend
    end subroutine spatial_operators_interface

    subroutine update_state_interface(dt, tend, old_state, new_state)
      import state_type
      import tend_type
      import real_kind
      real(real_kind),  intent(in)    :: dt
      type(tend_type),  intent(in)    :: tend
      type(state_type), intent(in)    :: old_state
      type(state_type), intent(inout) :: new_state
    end subroutine update_state_interface
  end interface

contains

  subroutine time_integrate(spatial_operators, update_state)

    procedure(spatial_operators_interface), pointer :: spatial_operators
    procedure(update_state_interface), pointer :: update_state

    select case (time_scheme)
    case (1)
      call time_integrate_CRK3(old, new, spatial_operators, update_state)
    case (2)
      call time_integrate_CRK4(old, new, spatial_operators, update_state)
    case (3)
      call time_integrate_PC2(old, new, spatial_operators, update_state)
    case default
      stop 'Unknown integration scheme'
    end select

  end subroutine time_integrate

  subroutine time_integrate_CRK3(old, new, spatial_operators, update_state)

    integer, intent(in) :: old
    integer, intent(in) :: new
    procedure(spatial_operators_interface), pointer :: spatial_operators
    procedure(update_state_interface), pointer :: update_state
    
    integer             :: one   = -1, &
                           two   = -2, &
                           three = -3
    
    real(real_kind) phi3_norm2, R1R2, R1R3, R2R3, beta

    state(one)%cell%gd    = state(old)%cell%gd
    state(one)%edge%iap_u = state(old)%edge%iap_u
    
    call spatial_operators(state(one), tend(one))
    call update_state(0.5d0 * dt, tend(one), state(one), state(two))

    call spatial_operators(state(two), tend(two))
    call update_state(-dt    , tend(one), state(one  ), state(three))
    call update_state(2.d0*dt, tend(two), state(three), state(three))

    call spatial_operators(state(three), tend(three))


    tend(new)%cell%gd    = (tend(one)%cell%gd    + 4.d0*tend(two)%cell%gd    + tend(three)%cell%gd   )/6.d0
    tend(new)%edge%iap_u = (tend(one)%edge%iap_u + 4.d0*tend(two)%edge%iap_u + tend(three)%edge%iap_u)/6.d0
    
    if (conserve_energy) then
      phi3_norm2 = inner_product(tend(new  ),tend(new  ))
      R1R2       = inner_product(tend(one  ),tend(two  ))
      R1R3       = inner_product(tend(one  ),tend(three))
      R2R3       = inner_product(tend(two  ),tend(three))
      
      beta       = (2.d0*R1R2 - R1R3 + 2.d0*R2R3)/(3.d0*phi3_norm2)
      call log_add_diag('beta', beta)
    else
      beta = 1.0d0
    end if
    
    call update_state(beta * dt, tend(new), state(old), state(new))

  end subroutine time_integrate_CRK3

  subroutine time_integrate_CRK4(old, new, spatial_operators, update_state)

    integer, intent(in) :: old
    integer, intent(in) :: new
    procedure(spatial_operators_interface), pointer :: spatial_operators
    procedure(update_state_interface), pointer :: update_state

    integer             :: one   = -1,&
                           two   = -2,&
                           three = -3,&
                           four  = -4
    real(real_kind) phi4_norm2, R1R2, R2R3, R3R4, beta

    state(one)%cell%gd    = state(old)%cell%gd
    state(one)%edge%iap_u = state(old)%edge%iap_u
    
    call spatial_operators(state(one), tend(one))
    call update_state(0.5d0 * dt, tend(one), state(one), state(two))

    call spatial_operators(state(two), tend(two))
    call update_state(0.5d0 * dt, tend(two), state(one), state(three))

    call spatial_operators(state(three), tend(three))
    call update_state(dt, tend(three), state(one), state(four))
    
    call spatial_operators(state(four), tend(four))

    tend(new)%cell%gd    = (tend(one)%cell%gd    + 2.d0*tend(two)%cell%gd    + 2.d0*tend(three)%cell%gd    + tend(four)%cell%gd   )/6.d0
    tend(new)%edge%iap_u = (tend(one)%edge%iap_u + 2.d0*tend(two)%edge%iap_u + 2.d0*tend(three)%edge%iap_u + tend(four)%edge%iap_u)/6.d0

    if (conserve_energy) then
      phi4_norm2 = inner_product(tend(new  ),tend(new  ))
      R1R2       = inner_product(tend(one  ),tend(two  ))
      R2R3       = inner_product(tend(two  ),tend(three))
      R3R4       = inner_product(tend(three),tend(four ))
      
      beta       = (R1R2 + R2R3 + R3R4)/(3.d0*phi4_norm2)
      call log_add_diag('beta', beta)
    else
      beta = 1.0d0
    end if
    
    call update_state(beta * dt, tend(new), state(old), state(new))

  end subroutine time_integrate_CRK4

  subroutine time_integrate_PC2(old, new, spatial_operators, update_state)

    integer, intent(in) :: old
    integer, intent(in) :: new
    procedure(spatial_operators_interface), pointer :: spatial_operators
    procedure(update_state_interface), pointer :: update_state

    real(real_kind) ip1, ip2, beta

    call spatial_operators(state(old), tend(old))
    call update_state(0.5d0 * dt, tend(old), state(old), state(new))

    call spatial_operators(state(new), tend(old))
    call update_state(0.5d0 * dt, tend(old), state(old), state(new))

    call spatial_operators(state(new), tend(new))

    if (conserve_energy) then
      ip1 = inner_product(tend(old), tend(new))
      ip2 = inner_product(tend(new), tend(new))
      beta = merge(ip1 / ip2, 1.0d0, ip1 /= 0.0d0 .and. ip2 /= 0.0d0)
      call log_add_diag('beta', beta)
    else
      beta = 1.0d0
    end if
    
    call update_state(beta * dt, tend(new), state(old), state(new))

  end subroutine time_integrate_PC2

  ! Degree of freedom to keep energy conservation
  subroutine calc_kinetic_energy(u_edge, ke_cell)

    real(real_kind), intent(in)  :: u_edge (:)
    real(real_kind), intent(out) :: ke_cell(:)

    integer iCell

    do iCell = lbound(ke_cell, 1), ubound(ke_cell, 1)
      ke_cell(iCell) = sum(                                                       &
                            u_edge  (edgesOnCell(1:nEdgesOnCell(iCell),iCell))**2 &
                          * areaEdge(edgesOnCell(1:nEdgesOnCell(iCell),iCell))    &
                          ) * 0.25d0
    end do
    ke_cell = ke_cell / areaCell

  end subroutine calc_kinetic_energy

  subroutine calc_tangent_wind(u_edge, v_edge)

    real(real_kind), intent(in)  :: u_edge(:)
    real(real_kind), intent(out) :: v_edge(:)

    integer iEdge

    do iEdge = lbound(v_edge, 1), ubound(v_edge, 1)
      v_edge(iEdge) = sum(                                                         &
                           u_edge       (edgesOnEdge(1:nEdgesOnEdge(iEdge),iEdge)) &
                         * weightsOnEdge            (1:nEdgesOnEdge(iEdge),iEdge)  &
                         )
    end do

  end subroutine calc_tangent_wind

  subroutine calc_tangent_vor_flux(u_edge, gd_edge, pv_edge, tangent_vor_flux)

    real(real_kind), intent(in)  :: u_edge          (:)
    real(real_kind), intent(in)  :: gd_edge         (:)
    real(real_kind), intent(in)  :: pv_edge         (:)
    real(real_kind), intent(out) :: tangent_vor_flux(:)

    integer iEdge, iEdgePrime, i

    tangent_vor_flux = 0.0d0
    do iEdge = lbound(u_edge, 1), ubound(u_edge, 1)
      do i = 1, nEdgesOnEdge(iEdge)
        iEdgePrime = edgesOnEdge(i,iEdge)

        tangent_vor_flux(iEdge) = tangent_vor_flux(iEdge) + weightsOnEdge(i,iEdge)                         &
                                                          * u_edge (iEdgePrime)                            &
                                                          * gd_edge(iEdgePrime)                            &
                                                          * 0.5d0 * (pv_edge(iEdge) + pv_edge(iEdgePrime))
      end do
    end do

  end subroutine calc_tangent_vor_flux

  subroutine calc_pv_on_vertex(vor_vertex, gd_vertex, pv_vertex)

    real(real_kind), intent(in)  :: vor_vertex(:)
    real(real_kind), intent(in)  :: gd_vertex (:)
    real(real_kind), intent(out) :: pv_vertex (:)

    pv_vertex = (fVertex + vor_vertex) / gd_vertex

  end subroutine calc_pv_on_vertex

  ! Degree of freedom to keep energy conservation
  subroutine calc_pv_on_edge(u_edge, v_edge, gd_edge, gd_tend_vertex, pv_vertex, pv_cell, pv_edge)

    real(real_kind), intent(in)  :: u_edge        (:)
    real(real_kind), intent(in)  :: v_edge        (:)
    real(real_kind), intent(in)  :: gd_edge       (:)
    real(real_kind), intent(in)  :: gd_tend_vertex(:)
    real(real_kind), intent(in)  :: pv_vertex     (:)
    real(real_kind), intent(in)  :: pv_cell       (:)
    real(real_kind), intent(out) :: pv_edge       (:)

    select case (pv_scheme)
    case (1)
      call calc_pv_on_edge_APVM(u_edge, v_edge, pv_vertex, pv_cell, pv_edge)
    case (2)
      call calc_pv_on_edge_CLUST()
    case (3)
      call calc_pv_on_edge_LUST()
    case (4)
      call calc_pv_on_edge_conservative_APVM(u_edge, v_edge, gd_edge, gd_tend_vertex, pv_vertex, pv_cell, pv_edge)
    case (5)
      call calc_pv_on_edge_order2()
    case (6)
      call calc_pv_on_edge_smoothed_order2()
    case default
      call log_error('Unknow PV scheme, please choose from 1 (APVM), 2 (CLUST), 3 (LUST), 4 (conservative APVM), 5 (order2), 6 (smoothed order2)!')
    end select

  end subroutine calc_pv_on_edge

  subroutine calc_pv_on_edge_APVM(u_edge, v_edge, pv_vertex, pv_cell, pv_edge)

    real(real_kind), intent(in)  :: u_edge   (:)
    real(real_kind), intent(in)  :: v_edge   (:)
    real(real_kind), intent(in)  :: pv_vertex(:)
    real(real_kind), intent(in)  :: pv_cell  (:)
    real(real_kind), intent(out) :: pv_edge  (:)

    integer iEdge

    call scalar_v2e_interp_operator(pv_vertex, pv_edge)

    if (apvm_weight /= 0.0d0) then
      do iEdge = lbound(pv_edge, 1), ubound(pv_edge, 1)
        pv_edge(iEdge) = pv_edge(iEdge) - apvm_weight * dt * ( v_edge(iEdge) * (pv_vertex(verticesOnEdge(2,iEdge)) - pv_vertex(verticesOnEdge(1,iEdge))) / dvEdge(iEdge) &
                                                             + u_edge(iEdge) * (pv_cell  (cellsOnEdge   (2,iEdge)) - pv_cell  (cellsOnEdge   (1,iEdge))) / dcEdge(iEdge) &
                                                             )
      end do
    end if

  end subroutine calc_pv_on_edge_APVM

  subroutine calc_pv_on_edge_conservative_APVM(u_edge, v_edge, gd_edge, gd_tend_vertex, pv_vertex, pv_cell, pv_edge)

    real(real_kind), intent(in)  :: u_edge        (:)
    real(real_kind), intent(in)  :: v_edge        (:)
    real(real_kind), intent(in)  :: gd_edge       (:)
    real(real_kind), intent(in)  :: gd_tend_vertex(:)
    real(real_kind), intent(in)  :: pv_vertex     (:)
    real(real_kind), intent(in)  :: pv_cell       (:)
    real(real_kind), intent(out) :: pv_edge       (:)

    real(real_kind) dpv_edge     (lbound(pv_edge, 1):ubound(pv_edge, 1))
    real(real_kind)  pv_flux_edge(lbound(pv_edge, 1):ubound(pv_edge, 1))
    real(real_kind) dpv_flux_edge(lbound(pv_edge, 1):ubound(pv_edge, 1))
    real(real_kind)  vor_tend_vertex(lbound(pv_vertex, 1):ubound(pv_vertex, 1))
    real(real_kind) dvor_tend_vertex(lbound(pv_vertex, 1):ubound(pv_vertex, 1))

    integer iEdge
    real(real_kind) eps

    call scalar_v2e_interp_operator(pv_vertex, pv_edge)

    do iEdge = lbound(pv_edge, 1), ubound(pv_edge, 1)
      dpv_edge(iEdge) = ( v_edge(iEdge) * (pv_vertex(verticesOnEdge(2,iEdge)) - pv_vertex(verticesOnEdge(1,iEdge))) / dvEdge(iEdge) &
                        + u_edge(iEdge) * (pv_cell  (cellsOnEdge   (2,iEdge)) - pv_cell  (cellsOnEdge   (1,iEdge))) / dcEdge(iEdge) &
                        )
    end do

    call calc_tangent_vor_flux(u_edge, gd_edge,  pv_edge,  pv_flux_edge)
    call calc_tangent_vor_flux(u_edge, gd_edge, dpv_edge, dpv_flux_edge)

    call calc_vor_tend_on_vertex( pv_flux_edge,  vor_tend_vertex)
    call calc_vor_tend_on_vertex(dpv_flux_edge, dvor_tend_vertex)

    eps = -( sum(pv_vertex**2 * gd_tend_vertex  * areaTriangle)         &
           - sum(pv_vertex    * vor_tend_vertex * areaTriangle) * 2.0d0 &
           ) / (2.0d0 * dt * sum(pv_vertex * dvor_tend_vertex * areaTriangle))

    pv_edge = pv_edge - eps * dt * dpv_edge

  end subroutine calc_pv_on_edge_conservative_APVM

  subroutine calc_pv_on_edge_CLUST()

    ! Under construction

  end subroutine calc_pv_on_edge_CLUST

  subroutine calc_pv_on_edge_LUST()

    ! Under construction

  end subroutine calc_pv_on_edge_LUST

  subroutine calc_pv_on_edge_order2()

    ! Under construction

  end subroutine calc_pv_on_edge_order2

  subroutine calc_pv_on_edge_smoothed_order2()

    ! Under construction

  end subroutine calc_pv_on_edge_smoothed_order2

  subroutine calc_vor_tend_on_vertex(u_tend_edge, vor_tend_vertex)

    real(real_kind), intent(in)  :: u_tend_edge    (:)
    real(real_kind), intent(out) :: vor_tend_vertex(:)

    call curl_operator(u_tend_edge, vor_tend_vertex)

  end subroutine calc_vor_tend_on_vertex

end module tmcore_common_mod
