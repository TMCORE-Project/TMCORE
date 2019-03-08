module test_cases_mod

  use params_mod
  use log_mod
  use mountain_zonal_flow_test_mod

  implicit none

contains

  subroutine test_cases_set_initial_condition()

    select case (case_name)
    case ('mountain_zonal_flow')
      call mountain_zonal_flow_test_set_initial_condition()
    case default
      call log_error('Unknown test case ' // trim(case_name) // '!')
    end select

  end subroutine test_cases_set_initial_condition

end module test_cases_mod
