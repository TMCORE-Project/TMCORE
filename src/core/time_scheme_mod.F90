module time_scheme_mod

  use params_mod
  use operators_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use log_mod
  use state_mod
  use tend_mod

  implicit none

  private

  public time_integrate
  public spatial_operators_interface
  public update_state_interface
  public update_wind_interface

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

    subroutine update_wind_interface(seconds, state)
      import state_type
      import real_kind
      real(real_kind),  intent(in)    :: seconds
      type(state_type), intent(inout) :: state
    end subroutine update_wind_interface
  end interface

contains

  subroutine time_integrate(spatial_operators, update_state)

    procedure(spatial_operators_interface) spatial_operators
    procedure(update_state_interface)      update_state

    select case (time_scheme)
    case (1)
      call time_integrate_CRK3(old, new, spatial_operators, update_state)
    case (2)
      call time_integrate_CRK4(old, new, spatial_operators, update_state)
    case (3)
      call time_integrate_PC2 (old, new, spatial_operators, update_state)
    case default
      stop 'Unknown integration scheme'
    end select

  end subroutine time_integrate

  subroutine time_integrate_CRK3(old, new, spatial_operators, update_state)

    integer, intent(in) :: old
    integer, intent(in) :: new
    procedure(spatial_operators_interface) spatial_operators
    procedure(update_state_interface)      update_state

    integer, parameter :: one   = -1
    integer, parameter :: two   = -2
    integer, parameter :: three = -3
    
    real(real_kind) R4R4, R1R2, R1R3, R2R3, beta

    state(one)%cell%gd    = state(old)%cell%gd
    state(one)%edge%iap_u = state(old)%edge%iap_u

    call spatial_operators(state(one), tend(one))
    call update_state(0.5d0 * dt, tend(one), state(one  ), state(two  ))

    call spatial_operators(state(two), tend(two))
    call update_state(      - dt, tend(one), state(one  ), state(three))
    call update_state(2.0d0 * dt, tend(two), state(three), state(three))

    call spatial_operators(state(three), tend(three))

    tend(new)%cell%pv    = (tend(one)%cell%pv    + 4.0d0 * tend(two)%cell%pv    + tend(three)%cell%pv   ) / 6.0d0
    tend(new)%cell%gd    = (tend(one)%cell%gd    + 4.0d0 * tend(two)%cell%gd    + tend(three)%cell%gd   ) / 6.0d0
    tend(new)%edge%iap_u = (tend(one)%edge%iap_u + 4.0d0 * tend(two)%edge%iap_u + tend(three)%edge%iap_u) / 6.0d0
    
    if (conserve_energy) then
      R4R4 = inner_product(tend(new), tend(new  ))
      R1R2 = inner_product(tend(one), tend(two  ))
      R1R3 = inner_product(tend(one), tend(three))
      R2R3 = inner_product(tend(two), tend(three))
      
      beta = (2.0d0 * R1R2 - R1R3 + 2.0d0 * R2R3) / (3.0d0 * R4R4)
      call log_add_diag('beta', beta)
    else
      beta = 1.0d0
    end if
    
    call update_state(beta * dt, tend(new), state(old), state(new))

  end subroutine time_integrate_CRK3

  subroutine time_integrate_CRK4(old, new, spatial_operators, update_state)

    integer, intent(in) :: old
    integer, intent(in) :: new
    procedure(spatial_operators_interface) spatial_operators
    procedure(update_state_interface)      update_state

    integer, parameter :: one   = -1
    integer, parameter :: two   = -2
    integer, parameter :: three = -3
    integer, parameter :: four  = -4
    real(real_kind) R4R4, R1R2, R2R3, R3R4, beta

    state(one)%cell%gd    = state(old)%cell%gd
    state(one)%cell%pv    = state(old)%cell%pv
    state(one)%edge%u     = state(old)%edge%u
    state(one)%edge%iap_u = state(old)%edge%iap_u

    call spatial_operators(state(one), tend(one))
    call update_state(0.5d0 * dt, tend(one), state(one), state(two))

    call spatial_operators(state(two), tend(two))
    call update_state(0.5d0 * dt, tend(two), state(one), state(three))

    call spatial_operators(state(three), tend(three))
    call update_state(dt, tend(three), state(one), state(four))
    
    call spatial_operators(state(four), tend(four))

    tend(new)%cell%pv    = (tend(one)%cell%pv    + 2.0d0 * tend(two)%cell%pv    + 2.0d0 * tend(three)%cell%pv    + tend(four)%cell%pv   ) / 6.0d0
    tend(new)%cell%gd    = (tend(one)%cell%gd    + 2.0d0 * tend(two)%cell%gd    + 2.0d0 * tend(three)%cell%gd    + tend(four)%cell%gd   ) / 6.0d0
    tend(new)%edge%iap_u = (tend(one)%edge%iap_u + 2.0d0 * tend(two)%edge%iap_u + 2.0d0 * tend(three)%edge%iap_u + tend(four)%edge%iap_u) / 6.0d0

    if (conserve_energy) then
      R4R4 = inner_product(tend(new  ), tend(new  ))
      R1R2 = inner_product(tend(one  ), tend(two  ))
      R2R3 = inner_product(tend(two  ), tend(three))
      R3R4 = inner_product(tend(three), tend(four ))
      
      beta = (R1R2 + R2R3 + R3R4) / (3.d0 * R4R4)
      call log_add_diag('beta', beta)
    else
      beta = 1.0d0
    end if
    
    call update_state(beta * dt, tend(new), state(old), state(new))

  end subroutine time_integrate_CRK4

  subroutine time_integrate_PC2(old, new, spatial_operators, update_state)

    integer, intent(in) :: old
    integer, intent(in) :: new
    procedure(spatial_operators_interface) spatial_operators
    procedure(update_state_interface)      update_state

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

end module time_scheme_mod