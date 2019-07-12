module test_cases_mod

  use params_mod
  use log_mod
  use state_mod
  use mountain_zonal_flow_test_mod
  use steady_geostrophic_flow_test_mod
  use rossby_haurwitz_wave_test_mod
  use linearized_rossby_wave_test_mod
  use jet_zonal_flow_test_mod
  use cosine_bell_rotation_test_mod
  use deformation_flow_test_mod
  
  implicit none

contains

  subroutine test_cases_set_initial_condition()

    select case (case_name)
    case ('mountain_zonal_flow')
      call mountain_zonal_flow_test_set_initial_condition()
    case ('steady_geostrophic_flow')
      call steady_geostrophic_flow_test_set_initial_condition()
    case ('rossby_haurwitz_wave')
      call rossby_haurwitz_wave_test_set_initial_condition()
    case ('linearized_rossby_wave')
      call linearized_rossby_wave_test_set_initial_condition()
    case ('jet_zonal_flow')
      call jet_zonal_flow_test_set_initial_condition()
    case ('cosine_bell_rotation')
      call cosine_bell_rotation_test_set_initial_condition()
    case ('deformation_flow')
      call deformation_flow_test_set_initial_condition()
    case default
      call log_error('Unknown test case ' // trim(case_name) // '!')
    end select

  end subroutine test_cases_set_initial_condition

  subroutine test_cases_update_wind(seconds, state)

    real(real_kind),  intent(in)    :: seconds
    type(state_type), intent(inout) :: state

    select case (case_name)
    case ('cosine_bell_rotation')
      call cosine_bell_rotation_test_update_wind(seconds, state)
    case ('deformation_flow')
      call deformation_flow_test_update_wind(seconds, state)
    case default
      call log_error('Unknown test case ' // trim(case_name) // '!')
    end select

  end subroutine test_cases_update_wind

end module test_cases_mod
