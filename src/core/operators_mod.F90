module operators_mod

  use params_mod
  use mesh_mod
  use static_mod
  use state_mod
  use tend_mod

  implicit none

  private

  public div_operator
  public grad_operator
  public curl_operator
  public scalar_c2e_interp_operator
  public scalar_c2v_interp_operator
  public scalar_v2c_interp_operator
  public scalar_v2e_interp_operator
  public iap_sw_operator
  public inverse_iap_sw_operator
  public inner_product

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
      div(iCell) = -sum(                                                     &
                         nSignEdge            (1:nEdgesOnCell(iCell),iCell)  &
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

    real(real_kind), intent(in)  :: f_cell(:)
    real(real_kind), intent(out) :: f_edge(:)

    integer iEdge

    do iEdge = lbound(f_edge, 1), ubound(f_edge, 1)
      f_edge(iEdge) = 0.5d0 * (                              &
                                f_cell(cellsOnEdge(1,iEdge)) &
                              + f_cell(cellsOnEdge(2,iEdge)) &
                              )
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

  subroutine iap_sw_operator(gd_edge, u_edge, iap_u_edge)

    real(real_kind), intent(in)  :: gd_edge   (:)
    real(real_kind), intent(in)  :: u_edge    (:)
    real(real_kind), intent(out) :: iap_u_edge(:)

    iap_u_edge = sqrt(gd_edge) * u_edge

  end subroutine iap_sw_operator

  subroutine inverse_iap_sw_operator(gd_edge, iap_u_edge, u_edge)

    real(real_kind), intent(in)  :: gd_edge   (:)
    real(real_kind), intent(in)  :: iap_u_edge(:)
    real(real_kind), intent(out) :: u_edge    (:)

    u_edge = iap_u_edge / sqrt(gd_edge)

  end subroutine inverse_iap_sw_operator

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

end module operators_mod
