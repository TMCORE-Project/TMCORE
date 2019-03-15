module tmcore_adv_mod

  use params_mod
  use log_mod
  use mesh_mod
  use operators_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use timer_mod
  use static_mod
  use state_mod
  use tend_mod
  use operators_mod
  use diag_mod
  use history_mod
  use time_scheme_mod
  use adv_scheme_mod

  implicit none

  private

  public tmcore_adv_init
  public tmcore_adv_final
  public tmcore_adv_run

  procedure(update_wind_interface), pointer :: update_wind_ptr

contains

  subroutine tmcore_adv_init(namelist_file_path, update_wind)

    character(*), intent(in) :: namelist_file_path
    procedure(update_wind_interface) update_wind

    update_wind_ptr => update_wind

    call params_parse_namelist(namelist_file_path)
    call log_init()
    call time_init()
    call timer_init()
    call mesh_init()
    call adv_scheme_init()
    call static_init()
    call state_init()
    call tend_init()
    call diag_init(calc_total_mass)
    call history_init()

  end subroutine tmcore_adv_init

  subroutine tmcore_adv_final()

    call adv_scheme_final()
    call mesh_final()
    call static_final()
    call state_final()
    call tend_final()
    call history_final()

  end subroutine tmcore_adv_final

  subroutine tmcore_adv_run()
  
    call diag_run(state(old), static)
    call history_write(state(old), static)
    call log_step()

    do while (.not. time_is_finished())
      call time_integrate(spatial_operators, update_state)
      call time_advance()
      if (time_is_alerted('hist0.output')) call history_write(state(old), static)
      call log_step()
    end do

  end subroutine tmcore_adv_run

  subroutine spatial_operators(state, tend)

    type(state_type), intent(inout) :: state
    type(tend_type),  intent(inout) :: tend

    call adv_scheme_run(state%cell%gd, state%edge%u, state%edge%gd, tend%cell%gd)

  end subroutine spatial_operators

  subroutine update_state(dt, tend, old_state, new_state)

    real(real_kind),  intent(in)    :: dt
    type(tend_type),  intent(in)    :: tend
    type(state_type), intent(in)    :: old_state
    type(state_type), intent(inout) :: new_state

    new_state%cell%gd = old_state%cell%gd + dt * tend%cell%gd
    call update_wind_ptr(time_elapsed_seconds() + dt, new_state)

  end subroutine update_state

  subroutine calc_total_mass(state)

    type(state_type), intent(inout) :: state

    state%total_mass = sum(state%cell%gd * areaCell)

  end subroutine calc_total_mass

end module tmcore_adv_mod
