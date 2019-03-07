program test_tmcore_sw

  use params_mod
  use log_mod
  use tmcore_sw_mod
  use mountain_zonal_flow_test_mod

  implicit none

  integer start_clock, end_clock

  call system_clock(start_clock)

  call tmcore_sw_init('./namelist.tmcore_sw')

  select case (case_name)
  case ('mountain_zonal_flow')
    call mountain_zonal_flow_test_set_initial_condition()
  case default
    call log_error('Unknown test case ' // trim(case_name) // '!')
  end select

  call tmcore_sw_run()

  call tmcore_sw_final()

  call system_clock(end_clock)

end program test_tmcore_sw