program test_tmcore_swm

  use params_mod
  use log_mod
  use tmcore_swm_mod
  use test_cases_mod

  implicit none

  integer start_clock, end_clock

  call system_clock(start_clock)

  call tmcore_swm_init('./namelist.tmcore_swm')

  call test_cases_set_initial_condition()

  call tmcore_swm_run()

  call tmcore_swm_final()

  call system_clock(end_clock)

end program test_tmcore_swm