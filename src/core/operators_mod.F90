module operators_mod

  use params_mod
  use log_mod
  use mesh_mod
  use static_mod
  use state_mod
  use tend_mod
  use poly_fit_mod

  implicit none

  private

  public div_operator
  public grad_operator
  public curl_operator
  public scalar_c2e_interp_operator
  public scalar_c2v_interp_operator
  public scalar_v2c_interp_operator
  public scalar_v2e_interp_operator
  public iap_swm_operator
  public inverse_iap_swm_operator
  public inner_product
  public calc_kinetic_energy
  public calc_tangent_wind
  public calc_tangent_vor_flux
  public calc_pv_on_vertex
  public calc_pv_on_edge
  public add_upwind_correction_on_cell

  interface inner_product
    module procedure inner_product_state
    module procedure inner_product_tend
    module procedure inner_product_state_tend
  end interface inner_product

contains

  subroutine div_operator(f_edge, div)

    real(real_kind), intent(in)  :: f_edge(:)
    real(real_kind), intent(out) :: div   (:)

    integer iCell

    do iCell = lbound(div, 1), ubound(div, 1)
      div(iCell) = -sum( nSignEdge            (1:nEdgesOnCell(iCell),iCell)  &
                       * f_edge   (edgesOnCell(1:nEdgesOnCell(iCell),iCell)) &
                       * dvEdge   (edgesOnCell(1:nEdgesOnCell(iCell),iCell)) &
                       )
    end do
    div = div / areaCell

  end subroutine div_operator

  subroutine grad_operator(f_cell, grad)

    real(real_kind), intent(in)  :: f_cell(:)
    real(real_kind), intent(out) :: grad  (:)

    integer iCell, iEdge, i

    do iCell = lbound(f_cell, 1), ubound(f_cell, 1)
      do i = 1, nEdgesOnCell(iCell)
        iEdge = edgesOnCell(i,iCell)
        
        grad(iEdge) = -sum( nSignEdge         (1:2,iEdge)  &
                          * f_cell(cellsOnEdge(1:2,iEdge)) &
                          ) / dcEdge(iEdge)
      end do
    end do

  end subroutine grad_operator

  subroutine curl_operator(f_edge, curl)

    real(real_kind), intent(in)  :: f_edge(:)
    real(real_kind), intent(out) :: curl  (:)

    integer iVertex

    do iVertex = lbound(curl, 1), ubound(curl, 1)
      curl(iVertex) = -sum( tSignEdge           (:,iVertex)  &
                          * f_edge(edgesOnVertex(:,iVertex)) &
                          * dcEdge(edgesOnVertex(:,iVertex)) &
                          )
    end do
    curl = curl / areaTriangle

  end subroutine curl_operator

  subroutine scalar_c2e_interp_operator(f_cell, f_edge)

    real(real_kind), intent(in ) :: f_cell(:)
    real(real_kind), intent(out) :: f_edge(:)

    integer iEdge

    do iEdge = lbound(f_edge, 1), ubound(f_edge, 1)
      f_edge(iEdge) = 0.5d0 * (f_cell(cellsOnEdge(1,iEdge)) + f_cell(cellsOnEdge(2,iEdge)))
    end do

  end subroutine scalar_c2e_interp_operator

  subroutine scalar_c2v_interp_operator(f_cell, f_vertex)

    real(real_kind), intent(in)  :: f_cell  (:)
    real(real_kind), intent(out) :: f_vertex(:)

    integer iVertex, iCell, i

    f_vertex = 0.0d0
    do iVertex = lbound(f_vertex, 1), ubound(f_vertex, 1)
      do i = 1, nCellsOnVertex(iVertex)
        iCell = cellsOnVertex(i,iVertex)

        f_vertex(iVertex) = f_vertex(iVertex) + f_cell(iCell) * kiteAreasOnVertex(i,iVertex)
      end do
    end do
    f_vertex = f_vertex / areaTriangle

  end subroutine scalar_c2v_interp_operator

  subroutine scalar_v2c_interp_operator(f_vertex, f_cell)

    real(real_kind), intent(in)  :: f_vertex(:)
    real(real_kind), intent(out) :: f_cell  (:)

    integer iVertex, iCell, i

    f_cell = 0.0d0
    do iVertex = lbound(f_vertex, 1), ubound(f_vertex, 1)
      do i = 1, nCellsOnVertex(iVertex)
        iCell = cellsOnVertex(i,iVertex)

        f_cell(iCell) = f_cell(iCell) + f_vertex(iVertex) * kiteAreasOnVertex(i,iVertex) / areaCell(iCell)
      end do
    end do

  end subroutine scalar_v2c_interp_operator

  subroutine scalar_v2e_interp_operator(f_vertex, f_edge)

    real(real_kind), intent(in)  :: f_vertex(:)
    real(real_kind), intent(out) :: f_edge  (:)

    integer iEdge

    do iEdge = lbound(f_edge, 1), ubound(f_edge, 1)
      f_edge(iEdge) = 0.5d0 * ( f_vertex(verticesOnEdge(1,iEdge)) &
                              + f_vertex(verticesOnEdge(2,iEdge)) &
                              )
    end do

  end subroutine scalar_v2e_interp_operator

  subroutine iap_swm_operator(gd_edge, u_edge, iap_u_edge)

    real(real_kind), intent(in)  :: gd_edge   (:)
    real(real_kind), intent(in)  :: u_edge    (:)
    real(real_kind), intent(out) :: iap_u_edge(:)

    iap_u_edge = sqrt(gd_edge) * u_edge

  end subroutine iap_swm_operator

  subroutine inverse_iap_swm_operator(gd_edge, iap_u_edge, u_edge)

    real(real_kind), intent(in)  :: gd_edge   (:)
    real(real_kind), intent(in)  :: iap_u_edge(:)
    real(real_kind), intent(out) :: u_edge    (:)

    u_edge = iap_u_edge / sqrt(gd_edge)

  end subroutine inverse_iap_swm_operator

  real(real_kind) function inner_product_state(state1, state2) result(res)

    type(state_type), intent(in) :: state1
    type(state_type), intent(in) :: state2

    res = (sum(state1%edge%iap_u * state2%edge%iap_u * areaEdge) + sum(state1%cell%gd * state2%cell%gd * areaCell)) / totalArea

  end function inner_product_state

  real(real_kind) function inner_product_tend(tend1, tend2) result(res)

    type(tend_type), intent(in) :: tend1
    type(tend_type), intent(in) :: tend2

    res = (sum(tend1%edge%iap_u * tend2%edge%iap_u * areaEdge) + sum(tend1%cell%gd * tend2%cell%gd * areaCell)) / totalArea

  end function inner_product_tend

  real(real_kind) function inner_product_state_tend(state, tend) result(res)

    type(state_type), intent(in) :: state
    type(tend_type),  intent(in) :: tend

    res = (sum(state%edge%iap_u * tend%edge%iap_u * areaEdge)  + sum((state%cell%gd + static%cell%ghs) * tend%cell%gd * areaCell)) /totalArea

  end function inner_product_state_tend

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
      call calc_pv_on_edge_conservative_APVM(u_edge, v_edge, gd_edge, gd_tend_vertex, pv_vertex, pv_cell, pv_edge)
    case (3)
      call calc_pv_on_edge_high_order_upwind(u_edge, v_edge, gd_edge, gd_tend_vertex, pv_vertex, pv_edge)
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

  subroutine calc_pv_on_edge_high_order_upwind(u_edge, v_edge, gd_edge, gd_tend_vertex, pv_vertex, pv_edge)
    real(real_kind), intent(in)  :: u_edge        (:)
    real(real_kind), intent(in)  :: v_edge        (:)
    real(real_kind), intent(in)  :: gd_edge       (:)
    real(real_kind), intent(in)  :: gd_tend_vertex(:)
    real(real_kind), intent(in)  :: pv_vertex     (:)
    real(real_kind), intent(out) :: pv_edge       (:)
    
    real(real_kind) :: coef2 = 0.5d0
    real(real_kind) :: coef4 = 0.25d0

    real(real_kind) d2fdx2_vertex1, d2fdx2_vertex2 ! 2nd order derivatives
    real(real_kind) d4fdx4_vertex1, d4fdx4_vertex2 ! 4th order derivatives
    
    integer i, iEdge, iVertex1, iVertex2
    
    pv_edge = 0.5d0 * ( pv_vertex(verticesOnEdge(1,:)) + pv_vertex(verticesOnEdge(2,:)) )

    coef2 = 0.0d0
    if (interp_pv_order == 3 ) coef2 =  1.d0
    if (interp_pv_order >= 4 ) coef2 =  0.0d0
    if (interp_pv_order == 6 ) coef4 =  0.0d0

    if (interp_pv_order == 3 .or. interp_pv_order == 4) then
      do iEdge = lbound(u_edge, 1), ubound(u_edge, 1)
        iVertex1 = verticesOnEdge(1,iEdge)
        iVertex2 = verticesOnEdge(2,iEdge)

        ! Compute 2nd order derivatives.
        d2fdx2_vertex1 = sum( deriv2OnVertex(1:nFit2Vertices(iVertex1)-1,1,iEdge) * (pv_vertex(fit2Vertices(1:nFit2Vertices(iVertex1)-1,iVertex1)) - pv_vertex(iVertex1)) )
        d2fdx2_vertex2 = sum( deriv2OnVertex(1:nFit2Vertices(iVertex2)-1,2,iEdge) * (pv_vertex(fit2Vertices(1:nFit2Vertices(iVertex2)-1,iVertex2)) - pv_vertex(iVertex2)) )
        
        pv_edge(iEdge) = pv_edge(iEdge) - dvEdge(iEdge)**2 * ( (d2fdx2_vertex1 + d2fdx2_vertex2) &
                        - sign(1.0d0, v_edge(iEdge)) * coef2 * (d2fdx2_vertex2 - d2fdx2_vertex1) ) / 12.0d0
      end do
    else if (interp_pv_order == 5 .or. interp_pv_order == 6) then
      do iEdge = lbound(u_edge, 1), ubound(u_edge, 1)
        iVertex1 = verticesOnEdge(1,iEdge)
        iVertex2 = verticesOnEdge(2,iEdge)

        ! Compute 2nd order derivatives.
        d2fdx2_vertex1 = sum( deriv2OnVertex(1:nFit2Vertices(iVertex1)-1,1,iEdge) * (pv_vertex(fit2Vertices(1:nFit2Vertices(iVertex1)-1,iVertex1)) - pv_vertex(iVertex1)) )
        d2fdx2_vertex2 = sum( deriv2OnVertex(1:nFit2Vertices(iVertex2)-1,2,iEdge) * (pv_vertex(fit2Vertices(1:nFit2Vertices(iVertex2)-1,iVertex2)) - pv_vertex(iVertex2)) )
        ! Compute 4th order derivatives.
        d4fdx4_vertex1 = sum( deriv4OnVertex(1:nFit4Vertices(iVertex1)-1,1,iEdge) * (pv_vertex(fit4Vertices(1:nFit4Vertices(iVertex1)-1,iVertex1)) - pv_vertex(iVertex1)) )
        d4fdx4_vertex2 = sum( deriv4OnVertex(1:nFit4Vertices(iVertex2)-1,2,iEdge) * (pv_vertex(fit4Vertices(1:nFit4Vertices(iVertex2)-1,iVertex2)) - pv_vertex(iVertex2)) )
        
        pv_edge(iEdge) = pv_edge(iEdge) - dvEdge(iEdge)**2 *  (d2fdx2_vertex1 + d2fdx2_vertex2)  / 12.0d0 &
                                        + dvEdge(iEdge)**4 * ((d4fdx4_vertex1 + d4fdx4_vertex2)           &
                       - sign(1.0d0, v_edge(iEdge)) * coef4 * (d4fdx4_vertex2 - d4fdx4_vertex1)) / 60.d0
      end do
    end if

  end subroutine calc_pv_on_edge_high_order_upwind

  subroutine calc_vor_tend_on_vertex(u_tend_edge, vor_tend_vertex)

    real(real_kind), intent(in)  :: u_tend_edge    (:)
    real(real_kind), intent(out) :: vor_tend_vertex(:)

    call curl_operator(u_tend_edge, vor_tend_vertex)

  end subroutine calc_vor_tend_on_vertex
  
  subroutine add_upwind_correction_on_cell(f_cell, u_edge, f_edge)

    real(real_kind), intent(in   ) :: f_cell(:)
    real(real_kind), intent(in   ) :: u_edge(:)
    real(real_kind), intent(inout) :: f_edge(:)

    real(real_kind) :: coef2 = 0.5d0
    real(real_kind) :: coef4 = 0.25d0
    
    real(real_kind) d2fdx2_cell1, d2fdx2_cell2 ! 2nd order derivatives
    real(real_kind) d4fdx4_cell1, d4fdx4_cell2 ! 4th order derivatives
  
    integer i, iEdge, iCell1, iCell2

    coef2 = 0.0d0
    if (adv_order == 3               ) coef2 =  1.0d0
    if (adv_order == 3 .and. adv_mono) coef2 = 0.25d0
    if (adv_order >= 4               ) coef2 =  0.0d0
    if (adv_order == 6               ) coef4 =  0.0d0
    
    if (adv_order == 3 .or. adv_order == 4) then
      do iEdge = lbound(f_edge, 1), ubound(f_edge, 1)
        iCell1 = cellsOnEdge(1,iEdge)
        iCell2 = cellsOnEdge(2,iEdge)

        ! Compute 2nd order derivatives.
        d2fdx2_cell1 = sum( deriv2OnCell(1:nFit2Cells(iCell1)-1,1,iEdge) * (f_cell(fit2Cells(1:nFit2Cells(iCell1)-1,iCell1)) - f_cell(iCell1)) )
        d2fdx2_cell2 = sum( deriv2OnCell(1:nFit2Cells(iCell2)-1,2,iEdge) * (f_cell(fit2Cells(1:nFit2Cells(iCell2)-1,iCell2)) - f_cell(iCell2)) )
        
        f_edge(iEdge) = f_edge(iEdge) - dcEdge(iEdge)**2 * ( (d2fdx2_cell1 + d2fdx2_cell2) - sign(1.0d0, u_edge(iEdge)) * coef2 * (d2fdx2_cell2 - d2fdx2_cell1) ) / 12.0d0
      end do
    else if (adv_order == 5 .or. adv_order == 6) then
      do iEdge = lbound(f_edge, 1), ubound(f_edge, 1)
        iCell1 = cellsOnEdge(1,iEdge)
        iCell2 = cellsOnEdge(2,iEdge)

        ! Compute 2nd order derivatives.
        d2fdx2_cell1 = sum( deriv2OnCell(1:nFit2Cells(iCell1)-1,1,iEdge) * (f_cell(fit2Cells(1:nFit2Cells(iCell1)-1,iCell1)) - f_cell(iCell1)) )
        d2fdx2_cell2 = sum( deriv2OnCell(1:nFit2Cells(iCell2)-1,2,iEdge) * (f_cell(fit2Cells(1:nFit2Cells(iCell2)-1,iCell2)) - f_cell(iCell2)) )
        ! Compute 4th order derivatives.
        d4fdx4_cell1 = sum( deriv4OnCell(1:nFit4Cells(iCell1)-1,1,iEdge) * (f_cell(fit4Cells(1:nFit4Cells(iCell1)-1,iCell1)) - f_cell(iCell1)) )
        d4fdx4_cell2 = sum( deriv4OnCell(1:nFit4Cells(iCell2)-1,2,iEdge) * (f_cell(fit4Cells(1:nFit4Cells(iCell2)-1,iCell2)) - f_cell(iCell2)) )
        
        f_edge(iEdge) = f_edge(iEdge) - dcEdge(iEdge)**2 * (d2fdx2_cell1 + d2fdx2_cell2) / 12.0d0 &
                                      + dcEdge(iEdge)**4 * ((d4fdx4_cell1 + d4fdx4_cell2) &
                                                           - sign(1.0d0, u_edge(iEdge)) * coef4 * (d4fdx4_cell2 - d4fdx4_cell1)) / 60.0d0
      end do
    end if

  end subroutine add_upwind_correction_on_cell

end module operators_mod
