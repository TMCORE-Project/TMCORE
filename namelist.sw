&sw_model
    run_days                          = 0
    run_hours                         = 2
    run_minutes                       = 0
    run_seconds                       = 0
    history_interval                  = 3600
    config_dt                         = 60
    config_test_case                  = 6
    config_time_integration           = 'CRK4'
	config_energy_conservation_scheme = 2,       !choose energy conservation scheme from 1: tau_n = 2.d0*(phi,F)/(phi,phi) or 2: tau_n = beta_n*dt
    config_energy_conservation        = .true.
    config_enstrophy_conservation     = .false.
	config_PV_scheme                  = 'APVM'	 ! PV interpolate scheme choose from 'order2', 'order2_smooth','APVM', 'APVM_Conservation', 'LUST' or 'CLUST', defalut is 'APVM'
    config_apvm_upwinding             = 0.0d0
	config_CLUST_upwinding            = 0.25d0
/
&io
    config_output_file            = ''   ! 'output_case2_2562_600s_NEC_PC3.nc', if config_output_file is '', the model will output with defalut name
/
&mesh
    config_mesh_file = 'x1.163842.grid.nc'
/
