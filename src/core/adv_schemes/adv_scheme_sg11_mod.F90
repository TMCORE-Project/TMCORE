module adv_scheme_sg11_mod

  use params_mod
  use timer_mod
  use mesh_mod
  use sphere_geometry_mod
  use math_mod
  use operators_mod
  use poly_fit_mod

  implicit none

  private

  public adv_scheme_sg11_init
  public adv_scheme_sg11_final
  public adv_scheme_sg11_run

contains

  subroutine adv_scheme_sg11_init()

    call poly_fit_init()

  end subroutine adv_scheme_sg11_init

  subroutine adv_scheme_sg11_final()

    call poly_fit_final()

  end subroutine adv_scheme_sg11_final

  subroutine adv_scheme_sg11_run(f_cell, u_edge, f_edge, tend_f_cell)

    real(real_kind), intent(in ) :: f_cell     (:)
    real(real_kind), intent(in ) :: u_edge     (:)
    real(real_kind), intent(out) :: f_edge     (:)
    real(real_kind), intent(out) :: tend_f_cell(:)

    real(real_kind) flux(lbound(u_edge, 1):ubound(u_edge, 1))

    call scalar_c2e_interp_operator(f_cell, f_edge)
    call add_upwind_correction(f_cell, u_edge, f_edge)

    flux = u_edge * f_edge
    call div_operator(flux, tend_f_cell)

  end subroutine adv_scheme_sg11_run

  subroutine add_upwind_correction(f_cell, u_edge, f_edge)


    real(real_kind), intent(in ) :: f_cell(:)
    real(real_kind), intent(in ) :: u_edge(:)
    real(real_kind), intent(out) :: f_edge(:)

    real(real_kind) :: coef2 = 0.5d0
    real(real_kind) :: coef4 = 0.25d0
    
    real(real_kind) d2fdx2_cell1, d2fdx2_cell2 ! 2nd order derivatives
    real(real_kind) d4fdx4_cell1, d4fdx4_cell2 ! 4th order derivatives
    
    integer i, iEdge, iCell1, icell2

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
        d2fdx2_cell1 = 0.0d0
        do i = 1, nFitCellsOnCell(2,iCell1) - 1
          d2fdx2_cell1 = d2fdx2_cell1 + derivOnCell(i,1,2,iEdge) * (f_cell(fitCellsOnCell(i,2,iCell1)) - f_cell(iCell1))
        end do
        d2fdx2_cell2 = 0.0d0
        do i = 1, nFitCellsOnCell(2,iCell2) - 1
          d2fdx2_cell2 = d2fdx2_cell2 + derivOnCell(i,2,2,iEdge) * (f_cell(fitCellsOnCell(i,2,iCell2)) - f_cell(iCell2))
        end do
        
        f_edge(iEdge) = f_edge(iEdge) - dcEdge(iEdge)**2 * ( &
          (d2fdx2_cell1 + d2fdx2_cell2) - sign(1.0d0, u_edge(iEdge)) * coef2 * (d2fdx2_cell2 - d2fdx2_cell1) &
        ) / 12.0d0
      end do
    else if (adv_order == 5 .or. adv_order == 6) then
      do iEdge = lbound(f_edge, 1), ubound(f_edge, 1)
        iCell1 = cellsOnEdge(1,iEdge)
        iCell2 = cellsOnEdge(2,iEdge)

        ! Compute 2nd order derivatives.
        d2fdx2_cell1 = 0.0d0
        do i = 1, nFitCellsOnCell(2,iCell1) - 1
          d2fdx2_cell1 = d2fdx2_cell1 + derivOnCell(i,1,2,iEdge) * (f_cell(fitCellsOnCell(i,2,iCell1)) - f_cell(iCell1))
        end do
        d2fdx2_cell2 = 0.0d0
        do i = 1, nFitCellsOnCell(2,iCell2) - 1
          d2fdx2_cell2 = d2fdx2_cell2 + derivOnCell(i,2,2,iEdge) * (f_cell(fitCellsOnCell(i,2,iCell2)) - f_cell(iCell2))
        end do
        
        ! Compute 4th order derivatives.
        d4fdx4_cell1 = 0.0d0
        do i = 1, nFitCellsOnCell(4,iCell1) - 1
          d4fdx4_cell1 = d4fdx4_cell1 + derivOnCell(i,1,4,iEdge) * (f_cell(fitCellsOnCell(i,4,iCell1)) - f_cell(iCell1))
        end do
        d4fdx4_cell2 = 0.0d0
        do i = 1, nFitCellsOnCell(4,iCell2) - 1
          d4fdx4_cell2 = d4fdx4_cell2 + derivOnCell(i,2,4,iEdge) * (f_cell(fitCellsOnCell(i,4,iCell2)) - f_cell(iCell2))
        end do

        f_edge(iEdge) = f_edge(iEdge) - dcEdge(iEdge)**2 * (d2fdx2_cell1 + d2fdx2_cell2) / 12.0d0 &
          + dcEdge(iEdge)**4 * ((d4fdx4_cell1 + d4fdx4_cell2) - sign(1.0d0, u_edge(iEdge)) * coef4 * (d4fdx4_cell2 - d4fdx4_cell1)) / 60.0d0
      end do
    end if

  end subroutine add_upwind_correction

end module adv_scheme_sg11_mod