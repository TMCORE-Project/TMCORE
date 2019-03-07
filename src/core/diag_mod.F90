module diag_mod

  use params_mod
  use time_mod
  use log_mod
  use static_mod
  use state_mod

  private

  public diag_init
  public diag_run
  public calc_total_mass_interface
  public calc_total_energy_interface
  public calc_total_potential_enstropy_interface

  interface
    subroutine calc_total_mass_interface(state)
      import state_type
      type(state_type), intent(inout) :: state
    end subroutine calc_total_mass_interface

    subroutine calc_total_energy_interface(state, static)
      import state_type
      import static_type
      type(state_type),  intent(inout) :: state
      type(static_type), intent(in) :: static
    end subroutine calc_total_energy_interface

    subroutine calc_total_potential_enstropy_interface(state)
      import state_type
      type(state_type), intent(inout) :: state
    end subroutine calc_total_potential_enstropy_interface
  end interface

  procedure(calc_total_mass_interface), pointer :: calc_total_mass_ptr
  procedure(calc_total_energy_interface), pointer :: calc_total_energy_ptr
  procedure(calc_total_potential_enstropy_interface), pointer :: calc_total_potential_enstropy_ptr

contains

  subroutine diag_init(calc_total_mass, calc_total_energy, calc_total_potential_enstropy)

    procedure(calc_total_mass_interface), pointer :: calc_total_mass
    procedure(calc_total_energy_interface), pointer :: calc_total_energy
    procedure(calc_total_potential_enstropy_interface), pointer :: calc_total_potential_enstropy

    calc_total_mass_ptr => calc_total_mass
    calc_total_energy_ptr => calc_total_energy
    calc_total_potential_enstropy_ptr => calc_total_potential_enstropy

  end subroutine diag_init

  subroutine diag_run(state, static)

    type(state_type),  intent(inout) :: state
    type(static_type), intent(in) :: static

    call calc_total_mass_ptr(state)
    call calc_total_energy_ptr(state, static)
    call calc_total_potential_enstropy_ptr(state)

    call log_add_diag('total_mass', state%total_mass)
    call log_add_diag('total_energy', state%total_energy)
    call log_add_diag('total_potential_enstropy', state%total_potential_enstropy)

  end subroutine diag_run

end module diag_mod