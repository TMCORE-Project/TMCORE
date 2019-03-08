program test_tmcore_sw

  use params_mod
  use log_mod
  use tmcore_sw_mod
  use test_cases_mod

  implicit none

  integer start_clock, end_clock

  call system_clock(start_clock)

  call tmcore_sw_init('./namelist.tmcore_sw')

  call test_cases_set_initial_condition()

  call tmcore_sw_run()

  call tmcore_sw_final()

  call system_clock(end_clock)

end program test_tmcore_sw