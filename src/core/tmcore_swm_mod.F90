module tmcore_swm_mod

  use params_mod
  use log_mod
  use mesh_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use static_mod
  use state_mod
  use tend_mod
  use operators_mod
  use diag_mod
  use history_mod
  use time_scheme_mod
  use adv_scheme_mod
  use poly_fit_mod

  implicit none

  private

  public tmcore_swm_init
  public tmcore_swm_final
  public tmcore_swm_run

contains

  subroutine tmcore_swm_init(namelist_file_path)

    character(*), intent(in) :: namelist_file_path

    call params_parse_namelist(namelist_file_path)
    call log_init()
    call time_init()
    call mesh_init()
    call adv_scheme_init()
    call static_init()
    call state_init()
    call tend_init()
    call diag_init(calc_total_mass, calc_total_energy, calc_total_potential_enstropy, calc_total_absolute_vorticity)
    call history_init()

  end subroutine tmcore_swm_init

  subroutine tmcore_swm_final()

    call adv_scheme_final()
    call mesh_final()
    call static_final()
    call state_final()
    call tend_final()
    call history_final()

  end subroutine tmcore_swm_final

  subroutine tmcore_swm_run()
  
    call scalar_c2e_interp_operator   (state(old)%cell  %gd , state(old)%edge  %gd                          )
    call add_upwind_correction_on_cell(state(old)%cell  %gd , state(old)%edge  %iap_u, state(old)%edge%gd   )
    call iap_swm_operator             (state(old)%edge  %gd , state(old)%edge  %u    , state(old)%edge%iap_u)
    call scalar_c2v_interp_operator   (state(old)%cell  %gd , state(old)%vertex%gd                          )
    call curl_operator                (state(old)%edge  %u  , state(old)%vertex%vor                         )
    call calc_pv_on_vertex            (state(old)%vertex%vor, state(old)%vertex%gd   , state(old)%vertex%pv )

    call diag_run(state(old), static)
    call history_write(state(old), static)
    call log_step()

    do while (.not. time_is_finished())
      call time_integrate(spatial_operators, update_state)
      call time_advance()
      call scalar_c2e_interp_operator   (state(old)%cell  %gd , state(old)%edge  %gd                         )
      call add_upwind_correction_on_cell(state(old)%cell  %gd , state(old)%edge  %iap_u, state(old)%edge  %gd) ! iap_u has the same sign as u
      call inverse_iap_swm_operator     (state(old)%edge  %gd , state(old)%edge  %iap_u, state(old)%edge  %u )
      call scalar_c2v_interp_operator   (state(old)%cell  %gd , state(old)%vertex%gd                         )
      call curl_operator                (state(old)%edge  %u  , state(old)%vertex%vor                        )
      call calc_pv_on_vertex            (state(old)%vertex%vor, state(old)%vertex%gd   , state(old)%vertex%pv)
      call diag_run(state(old), static)
      if (time_is_alerted('hist0.output')) call history_write(state(old), static)
      call log_step()
    end do

  end subroutine tmcore_swm_run

  subroutine spatial_operators(state, tend)

    type(state_type), intent(inout) :: state
    type(tend_type),  intent(inout) :: tend

    ! Inverse IAP transformation
    call scalar_c2e_interp_operator   (state%cell  %gd     , state%edge  %gd                    )
    call add_upwind_correction_on_cell(state%cell  %gd     , state%edge  %iap_u, state%edge%gd  )
    call inverse_iap_swm_operator     (state%edge  %gd     , state%edge  %iap_u, state%edge%u   )
    
    ! Calculate tend of geopotential depth (gd)
    call calc_gd_tend_on_cell         (state%edge  %u      , state%edge  %gd   , tend %cell%gd  )
    
    ! Interpolate tend gd to edge and vertex
    call scalar_c2e_interp_operator   (tend %cell  %gd     , tend %edge  %gd                    )
    call add_upwind_correction_on_cell(tend %cell  %gd     , state%edge  %iap_u, tend %edge%gd  )
    call scalar_c2v_interp_operator   (tend %cell  %gd     , tend %vertex%gd                    )
    call scalar_c2v_interp_operator   (state%cell  %gd     , state%vertex%gd                    )
    
    ! Calculate potential vorticity to edge
    call curl_operator                (state%edge  %u      , state%vertex%vor                   )
    call calc_tangent_wind            (state%edge  %u      , state%edge  %v                     )
    call calc_kinetic_energy          (state%edge  %u      , state%cell  %ke                    )
    call calc_pv_on_vertex            (state%vertex%vor    , state%vertex%gd   , state%vertex%pv)
    call scalar_v2c_interp_operator   (state%vertex%pv     , state%cell  %pv                    )
    call calc_pv_on_edge              (state%edge  %u      , state%edge  %v    , state%edge  %gd  , tend %vertex%gd     , state %vertex%pv  , state%cell%pv   , state%edge%pv)
    
    ! Calculate tend of u and iap_u
    call calc_tangent_vor_flux        (state%edge  %u      , state%edge  %gd   , state%edge  %pv  , state%edge  %pv_flx                                                      )
    call calc_u_tend_on_edge          (state%edge  %u      , state%cell  %ke   , state%cell  %gd  , state%edge  %gd     , static%cell  %ghs ,                                &
                                       state%edge  %pv_flx , tend %edge  %gd   , tend %edge  %u   , tend %edge  %iap_u                                                       )

  end subroutine spatial_operators

  subroutine calc_gd_tend_on_cell(u_edge, gd_edge, gd_tend_cell)

    real(real_kind), intent(in)  :: u_edge      (:)
    real(real_kind), intent(in)  :: gd_edge     (:)
    real(real_kind), intent(out) :: gd_tend_cell(:)

    real(real_kind) flux(lbound(u_edge, 1):ubound(u_edge, 1))

    flux = u_edge * gd_edge
    
    call div_operator(flux, gd_tend_cell)

  end subroutine calc_gd_tend_on_cell

  subroutine calc_vorticity_tend_on_cell(u_edge, vor_edge, vor_tend_cell)

  real(real_kind), intent(in)  :: u_edge      (:)
  real(real_kind), intent(in)  :: vor_edge     (:)
  real(real_kind), intent(out) :: vor_tend_cell(:)

  real(real_kind) flux(lbound(u_edge, 1):ubound(u_edge, 1))

  flux = u_edge * vor_edge
  
  call div_operator(flux, vor_tend_cell)

  end subroutine calc_vorticity_tend_on_cell
  
  subroutine calc_u_tend_on_edge(u_edge, ke_cell, gd_cell, gd_edge, ghs_cell, pv_flx_edge, gd_tend_edge, u_tend_edge, iap_u_tend_edge)

    real(real_kind), intent(in)  :: u_edge         (:)
    real(real_kind), intent(in)  :: ke_cell        (:)
    real(real_kind), intent(in)  :: gd_cell        (:)
    real(real_kind), intent(in)  :: gd_edge        (:)
    real(real_kind), intent(in)  :: ghs_cell       (:)
    real(real_kind), intent(in)  :: pv_flx_edge    (:)
    real(real_kind), intent(in)  :: gd_tend_edge   (:)
    real(real_kind), intent(out) :: u_tend_edge    (:)
    real(real_kind), intent(out) :: iap_u_tend_edge(:)

    real(real_kind) iap_gd_edge(lbound(u_tend_edge, 1):ubound(u_tend_edge, 1))
    
    real(real_kind) E          (lbound(u_tend_edge, 1):ubound(u_tend_edge, 1))
    real(real_kind) dEdx       (lbound(u_tend_edge, 1):ubound(u_tend_edge, 1))
    
    real(real_kind) d3fdx3_cell1, &
                    d3fdx3_cell2
    
    integer iCell1, iCell2, iEdge
      
    E = ke_cell + gd_cell - ghs_cell
    
    iap_gd_edge = sqrt(gd_edge)
    
    dEdx  = ( E(cellsOnEdge(2,:)) - E(cellsOnEdge(1,:))) / dcEdge
    
    if (flux_4th_order_correct)then
      do iEdge = lbound(u_tend_edge, 1), ubound(u_tend_edge, 1)
        iCell1 = cellsOnEdge(1,iEdge)
        iCell2 = cellsOnEdge(2,iEdge)
        
        ! Compute 4th order derives
        d3fdx3_cell1 = sum( deriv3OnCell(1:nFit3Cells(iCell1)-1,1,iEdge) * (E(fit3Cells(1:nFit3Cells(iCell1)-1,iCell1)) - E(iCell1)) )
        d3fdx3_cell2 = sum( deriv3OnCell(1:nFit3Cells(iCell2)-1,2,iEdge) * (E(fit3Cells(1:nFit3Cells(iCell2)-1,iCell2)) - E(iCell2)) )
      
        dEdx(iEdge)  = dEdx(iEdge) - ( d3fdx3_cell2 - d3fdx3_cell1 ) * dcEdge(iEdge) / 24.d0
        
      end do     ! do iEdge
    end if

    u_tend_edge = pv_flx_edge - dEdx

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

end module tmcore_swm_mod
