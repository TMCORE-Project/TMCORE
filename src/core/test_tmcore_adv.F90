program test_tmcore_adv

  use params_mod
  use log_mod
  use tmcore_adv_mod
  use test_cases_mod

  implicit none

  integer start_clock, end_clock

  call system_clock(start_clock)

  call tmcore_adv_init('./namelist.tmcore_adv', test_cases_update_wind)

  call test_cases_set_initial_condition()

  call tmcore_adv_run()

  call tmcore_adv_final()

  call system_clock(end_clock)

end program test_tmcore_adv