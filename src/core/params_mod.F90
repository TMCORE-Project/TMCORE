module params_mod

  use const_mod

  implicit none

  integer run_days
  integer run_hours
  integer run_minutes
  integer run_seconds
  integer :: start_time(5) = [0, 0, 0, 0, 0]
  integer :: end_time  (5) = [0, 0, 0, 0, 0]
  real(real_kind) dt

  namelist /time_control/ &
    run_days,             &
    run_hours,            &
    run_minutes,          &
    run_seconds,          &
    start_time,           &
    end_time,             &
    dt

  character(max_name_len) history_periods(1)
  character(max_name_len) :: time_units = 'days'
  character(max_file_path_len) :: output_file_prefix = 'N/A'
  character(max_name_len) :: frames_per_file = 'N/A'

  namelist /io_control/   &
    history_periods,      &
    time_units,           &
    output_file_prefix,   &
    frames_per_file

  integer time_scheme                 ! Time integration scheme
  integer energy_scheme               ! Total energy conservation scheme, 1: tau_n = 2.d0*(phi,F)/(phi,phi); 2: tau_n = beta_n*dt
  integer pv_scheme                   ! 1: APVM; 2: CLUST 4: conservative_APVM
  logical :: conserve_energy = .true.
  real(real_kind) :: apvm_weight = 0.5
  real(real_kind) :: clust_weight = 0.25
  character(max_file_path_len) mesh_file_path

  namelist /tmcore/       &
    time_scheme,          &
    energy_scheme,        &
    pv_scheme,            &
    conserve_energy,      &
    apvm_weight,          &
    clust_weight,         &
    mesh_file_path

  character(max_name_len) case_name

  namelist /test_case/     &
    case_name

contains

  subroutine params_parse_namelist(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path)
    read(10, nml=time_control)
    read(10, nml=io_control)
    read(10, nml=tmcore)
    read(10, nml=test_case)
    close(10)

  end subroutine params_parse_namelist

end module params_mod
