MODULE module_para
    
    implicit none
    integer         ,  parameter ::    IKIND  = 8
    integer         ,  parameter ::    RKIND  = 8
    real(kind=RKIND),  parameter ::    omega  = 7.292d-5		 ! angular velocity of the earth rotation
	real(kind=RKIND),  parameter ::    a      = 6371220d0		 ! radius of the earth
    real(kind=RKIND),  parameter ::    g      = 9.80616d0		 ! gravity acceleration
	real(kind=RKIND),  parameter ::    pi     = datan(1.d0)*4.d0 ! pi = datan(1d0)*4.d0 = 3.141592653589793238462643384279d0

    integer(kind=IKIND)::    run_days                      
    integer(kind=IKIND)::    run_hours                     
    integer(kind=IKIND)::    run_minutes                   
    integer(kind=IKIND)::    run_seconds                   
    integer(kind=IKIND)::    history_interval                      ! output interval in seconds
    real(kind=RKIND)   ::    config_dt                             ! integral time step
    integer            ::    config_test_case                      ! select test case
    character*200      ::    config_time_integration               ! integral scheme
    logical            ::    config_energy_conservation            ! if turn on energy conservation option
    integer            ::    config_energy_conservation_scheme     ! choose energy conservation scheme from 1:tau_n = 2.d0*(phi,F)/(phi,phi) or 2:tau_n = beta_n*dt
    logical            ::    config_enstrophy_conservation         ! if turn on enstrophy conservation option
    character*500      ::    config_PV_scheme                      ! PV interpolate scheme choose from 'APVM' or 'CLUST', defalut is 'APVM'
    real(kind=RKIND)   ::    config_apvm_upwinding                 ! APVM parameter, default is 0.5
    real(kind=RKIND)   ::    config_CLUST_upwinding                ! CLUST parameter, default is 0.25
    character*500      ::    config_mesh_file                      ! mesh file
    character*500      ::    config_output_file                    ! output netCDF file name
    
    namelist /sw_model/ run_days                          ,&
                        run_hours                         ,&
                        run_minutes                       ,&
                        run_seconds                       ,&
                        history_interval                  ,&
                        config_dt                         ,&
                        config_test_case                  ,&
                        config_time_integration           ,&
                        config_energy_conservation        ,&
                        config_energy_conservation_scheme ,&
                        config_enstrophy_conservation     ,&
                        config_PV_scheme                  ,&
                        config_apvm_upwinding             ,&
                        config_CLUST_upwinding
    namelist /io/       config_output_file
    namelist /mesh/     config_mesh_file
    
    contains
    
    subroutine read_namelist
    implicit none
    
    ! Set the initial value
    run_days                      = 100
    run_hours                     = 0    
    run_minutes                   = 0    
    run_seconds                   = 0    
    history_interval              = 3600   ! output interval in seconds    
    config_dt                     = 400    ! integral time step    
    config_test_case              = 2      ! select test case    
    config_time_integration       = 'CRK4' ! integral scheme    
    config_energy_conservation    = .true.
    config_enstrophy_conservation = .false.
    config_PV_scheme              = 'APVM'
    config_apvm_upwinding         = 0.5    
    config_CLUST_upwinding        = 0.25
    config_mesh_file              = 'x1.2562.grid.nc' ! mesh file    
    
    open(1,file = 'namelist.sw',status='old')
    read(1,nml=sw_model)
    read(1,nml=io)
    read(1,nml=mesh)
    close(1)
    end subroutine read_namelist
    
END MODULE module_para
    
