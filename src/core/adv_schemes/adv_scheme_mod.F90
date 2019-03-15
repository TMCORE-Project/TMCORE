module adv_scheme_mod

  use params_mod
  use adv_scheme_sg11_mod

  implicit none

contains

  subroutine adv_scheme_init()

    call adv_scheme_sg11_init()

  end subroutine adv_scheme_init

  subroutine adv_scheme_final()

    call adv_scheme_sg11_final()

  end subroutine adv_scheme_final

  subroutine adv_scheme_run(f_cell, u_edge, f_edge, tend_f_cell)

    real(real_kind), intent(in ) :: f_cell     (:)
    real(real_kind), intent(in ) :: u_edge     (:)
    real(real_kind), intent(out) :: f_edge     (:)
    real(real_kind), intent(out) :: tend_f_cell(:)

    call adv_scheme_sg11_run(f_cell, u_edge, f_edge, tend_f_cell)

  end subroutine adv_scheme_run

end module adv_scheme_mod