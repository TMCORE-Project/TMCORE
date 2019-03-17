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
    call add_upwind_correction_on_cell(f_cell, u_edge, f_edge)

    flux = u_edge * f_edge
    call div_operator(flux, tend_f_cell)

  end subroutine adv_scheme_sg11_run

end module adv_scheme_sg11_mod