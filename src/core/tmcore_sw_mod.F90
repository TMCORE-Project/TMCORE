module tmcore_sw_mod

  use params_mod
  use log_mod
  use mesh_mod
  use advection_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use static_mod
  use state_mod
  use tend_mod
  use operators_mod
  use diag_mod
  use history_mod
  use time_scheme_mod

  implicit none

  private

  public tmcore_sw_init
  public tmcore_sw_final
  public tmcore_sw_run

  procedure(spatial_operators_interface), pointer :: spatial_operators_ptr
  procedure(update_state_interface), pointer :: update_state_ptr
  procedure(calc_total_mass_interface), pointer :: calc_total_mass_ptr
  procedure(calc_total_energy_interface), pointer :: calc_total_energy_ptr
  procedure(calc_total_potential_enstropy_interface), pointer :: calc_total_potential_enstropy_ptr
  procedure(calc_total_absolute_vorticity_interface), pointer :: calc_total_absolute_vorticity_ptr

contains

  subroutine tmcore_sw_init(namelist_file_path)

    character(*), intent(in) :: namelist_file_path

    spatial_operators_ptr => spatial_operators
    update_state_ptr => update_state
    calc_total_mass_ptr => calc_total_mass
    calc_total_energy_ptr => calc_total_energy
    calc_total_potential_enstropy_ptr => calc_total_potential_enstropy
    calc_total_absolute_vorticity_ptr => calc_total_absolute_vorticity

    call params_parse_namelist(namelist_file_path)
    call log_init()
    call time_init()
    call mesh_init()
    call advecion_init()
    call static_init()
    call state_init()
    call tend_init()
    call diag_init(calc_total_mass_ptr, calc_total_energy_ptr, calc_total_potential_enstropy_ptr, calc_total_absolute_vorticity_ptr)
    call history_init()

  end subroutine tmcore_sw_init

  subroutine tmcore_sw_final()

    call mesh_final()
    call static_final()
    call state_final()
    call tend_final()
    call history_final()

  end subroutine tmcore_sw_final

  subroutine tmcore_sw_run()
  
    call scalar_c2e_interp_operator(state(old)%cell  %gd   , state(old)%edge  %gd    , adv_order            , state(old)%edge%u, adv_monotonic)
    call iap_sw_operator           (state(old)%edge  %gd   , state(old)%edge  %u     , state(old)%edge%iap_u                                  )
    call scalar_c2v_interp_operator(state(old)%cell  %gd   , state(old)%vertex%gd                                                             )
    call curl_operator             (state(old)%edge  %u    , state(old)%vertex%vor                                                            )
    call calc_pv_on_vertex         (state(old)%vertex%vor  , state(old)%vertex%gd    , state(old)%vertex%pv                                   )

    call diag_run(state(old), static)
    call history_write(state(old), static)
    call log_step()

    do while (.not. time_is_finished())
      call time_integrate(spatial_operators_ptr, update_state_ptr)
      call time_advance()
      
      call scalar_c2e_interp_operator(state(old)%cell  %gd , state(old)%edge  %gd   , adv_order            , state(old)%edge%u, adv_monotonic)
      call inverse_iap_sw_operator   (state(old)%edge  %gd , state(old)%edge  %iap_u, state(old)%edge  %u                                    )
      call scalar_c2v_interp_operator(state(old)%cell  %gd , state(old)%vertex%gd                                                            )
      call curl_operator             (state(old)%edge  %u  , state(old)%vertex%vor                                                           )
      call calc_pv_on_vertex         (state(old)%vertex%vor, state(old)%vertex%gd   , state(old)%vertex%pv                                   )
      call diag_run(state(old), static)
      
      if (time_is_alerted('hist0.output')) call history_write(state(old), static)
      call log_step()
    end do

  end subroutine tmcore_sw_run

  subroutine spatial_operators(state, tend)

    type(state_type), intent(inout) :: state
    type(tend_type),  intent(inout) :: tend

    call scalar_c2e_interp_operator(state%cell  %gd     , state%edge  %gd   , adv_order        , state%edge%u, adv_monotonic)
    call inverse_iap_sw_operator   (state%edge  %gd     , state%edge  %iap_u, state%edge%u   )
    call calc_gd_tend_on_cell      (state%edge  %u      , state%edge  %gd   , tend %cell%gd  )
    call scalar_c2e_interp_operator(tend %cell  %gd     , tend %edge  %gd   , adv_order        , state%edge%u, adv_monotonic)
    call scalar_c2v_interp_operator(tend %cell  %gd     , tend %vertex%gd                    )
    call scalar_c2v_interp_operator(state%cell  %gd     , state%vertex%gd                    )
    call curl_operator             (state%edge  %u      , state%vertex%vor                   )
    call calc_tangent_wind         (state%edge  %u      , state%edge  %v                     )
    call calc_kinetic_energy       (state%edge  %u      , state%cell  %ke                    )
    call calc_pv_on_vertex         (state%vertex%vor    , state%vertex%gd   , state%vertex%pv)
    call scalar_v2c_interp_operator(state%vertex%pv     , state%cell  %pv                    )
    call calc_pv_on_edge           (state%edge  %u      , state%edge  %v    , state%edge  %gd  , tend %vertex%gd     , state %vertex%pv  , state%cell%pv   , state%edge%pv)
    call calc_tangent_vor_flux     (state%edge  %u      , state%edge  %gd   , state%edge  %pv  , state%edge  %pv_flx                                                      )
    call calc_u_tend_on_edge       (state%edge  %u      , state%cell  %ke   , state%cell  %gd  , state%edge  %gd     , static%cell  %ghs ,                                &
                                    state%edge  %pv_flx , state%vertex%pv   , tend %cell  %gd  , tend %edge  %gd     , tend  %edge  %u   , tend%edge%iap_u                )

  end subroutine spatial_operators

  subroutine calc_gd_tend_on_cell(u_edge, gd_edge, gd_tend_cell)

    real(real_kind), intent(in)  :: u_edge      (:)
    real(real_kind), intent(in)  :: gd_edge     (:)
    real(real_kind), intent(out) :: gd_tend_cell(:)

    real(real_kind) flux(lbound(u_edge, 1):ubound(u_edge, 1))
    integer iCell

    flux = u_edge * gd_edge
    
    call div_operator(flux, gd_tend_cell)

  end subroutine calc_gd_tend_on_cell

  subroutine calc_u_tend_on_edge(u_edge, ke_cell, gd_cell, gd_edge, ghs_cell, pv_flx_edge, pv_vertex, gd_tend_cell, gd_tend_edge, u_tend_edge, iap_u_tend_edge)

    real(real_kind), intent(in)  :: u_edge         (:)
    real(real_kind), intent(in)  :: ke_cell        (:)
    real(real_kind), intent(in)  :: gd_cell        (:)
    real(real_kind), intent(in)  :: gd_edge        (:)
    real(real_kind), intent(in)  :: ghs_cell       (:)
    real(real_kind), intent(in)  :: pv_flx_edge    (:)
    real(real_kind), intent(in)  :: pv_vertex      (:)
    real(real_kind), intent(in)  :: gd_tend_cell   (:)
    real(real_kind), intent(in)  :: gd_tend_edge   (:)
    real(real_kind), intent(out) :: u_tend_edge    (:)
    real(real_kind), intent(out) :: iap_u_tend_edge(:)

    real(real_kind) iap_gd_edge(lbound(u_tend_edge, 1):ubound(u_tend_edge, 1))
    real(real_kind) dkedx      (lbound(u_tend_edge, 1):ubound(u_tend_edge, 1))
    real(real_kind) dghdx      (lbound(u_tend_edge, 1):ubound(u_tend_edge, 1))

    iap_gd_edge = sqrt(gd_edge)
    dkedx = ( ke_cell(cellsOnEdge(2,:)) -  ke_cell(cellsOnEdge(1,:))) / dcEdge
    dghdx = ( gd_cell(cellsOnEdge(2,:)) -  gd_cell(cellsOnEdge(1,:)) + &
             ghs_cell(cellsOnEdge(2,:)) - ghs_cell(cellsOnEdge(1,:))) / dcEdge

    u_tend_edge = pv_flx_edge - dkedx - dghdx

    iap_u_tend_edge = iap_gd_edge * u_tend_edge + 0.5d0 * u_edge / iap_gd_edge * gd_tend_edge

  end subroutine calc_u_tend_on_edge

  subroutine update_state(dt, tend, old_state, new_state)

    real(real_kind),  intent(in)    :: dt
    type(tend_type),  intent(in)    :: tend
    type(state_type), intent(in)    :: old_state
    type(state_type), intent(inout) :: new_state

    new_state%edge%iap_u = old_state%edge%iap_u + dt * tend%edge%iap_u
    new_state%cell%gd    = old_state%cell%gd    + dt * tend%cell%gd

  end subroutine update_state

  subroutine calc_total_mass(state)

    type(state_type), intent(inout) :: state

    state%total_mass = sum(state%cell%gd * areaCell)

  end subroutine calc_total_mass

  subroutine calc_total_energy(state, static)

    type(state_type),  intent(inout) :: state
    type(static_type), intent(in   ) :: static

    state%total_energy = sum(state%edge%iap_u**2 * areaEdge) + sum((state%cell%gd + static%cell%ghs)**2 * areaCell)

  end subroutine calc_total_energy

  subroutine calc_total_potential_enstropy(state)

    type(state_type), intent(inout) :: state

    state%total_potential_enstropy = sum(state%vertex%gd * state%vertex%pv**2 * areaTriangle)

  end subroutine calc_total_potential_enstropy

  subroutine calc_total_angular_momentum(state)

    type(state_type), intent(inout) :: state

    state%total_angular_momentum = sum(state%edge%gd * (state%edge%u * cos(angleEdge) + radius * omega * cos(latEdge)) * cos(latEdge) * areaEdge) / sum(areaEdge)

  end subroutine calc_total_angular_momentum

  subroutine calc_total_absolute_vorticity(state)

    type(state_type), intent(inout) :: state

    state%total_absolute_vorticity = sum((state%vertex%vor + fVertex) * areaTriangle) / sum(areaTriangle)

  end subroutine calc_total_absolute_vorticity
end module tmcore_sw_mod
