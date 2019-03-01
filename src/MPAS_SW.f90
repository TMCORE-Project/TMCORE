!*****************************************************************************!
!			 The Barotropic Model on Arbitrarily-Structured C-grids	          !
!                                  by Lilong Zhou                             !
!*****************************************************************************!	
!
PROGRAM MAIN

    use module_para      , only: RKIND,read_namelist,config_dt,config_output_file,history_interval,run_days,run_hours,run_minutes,run_seconds
    use module_array     , only: init_arrays,wu,wh,wh_s,u
    use mesh_info        , only: read_mesh,areaEdge,areaCell,sum_areaEdge,sum_areaCell
    use test_cases       , only: choose_case
    use output_netCDF    , only: write_netCDF
    use diagnostics_tools, only: totalMass               ,&
                                 totalAbsoluteVorticity  ,&
                                 totalEnergy             ,&
                                 totalPotentialEnstrophy ,&
                                 totalAngularMomentum

    implicit none
    
    integer*8 nt
    integer*8 i,nstep
    integer   output_idx,output_num
    integer   time_1,time_2
    
    ! Diagnostic
    real(kind=RKIND) dtn
    real(kind=RKIND) Mass0               ,&
                     Mass                ,&
                     AbsoluteVotricity0  ,&
                     AbsoluteVotricity   ,&
                     Energy0             ,&
                     Energy              ,&
                     PotentialEnstrophy0 ,&
                     PotentialEnstrophy  ,&
                     AngularMomentum0    ,&
                     AngularMomentum
    
    
    call SYSTEM_CLOCK(time_1)
    
    call read_namelist
    call read_mesh
    call init_arrays
    call choose_case
    
    nt         = run_days*86400.d0 + run_hours*3600.d0 + run_minutes*60.d0 + run_seconds
    nstep      = int(nt/config_dt)
    
    output_idx = 0
    output_num = nt/history_interval
    
    print*,'Output time ',output_idx,'/',output_num
    call write_netCDF(wu,wh,wh_s,output_idx)
    
    Mass0               = totalMass(wh)
    Energy0             = totalEnergy(wu,wh,wh_s)
    AbsoluteVotricity0  = totalAbsoluteVorticity(wu,wh)
    PotentialEnstrophy0 = totalPotentialEnstrophy(wu,wh)
    AngularMomentum0    = totalAngularMomentum(wu,wh)
    
    i=0
    do i = 1,nstep
        call integration(dtn)
        
        if(mod(int(i*config_dt,8),history_interval)==0)then
            output_idx = output_idx + 1
            call write_netCDF(wu,wh,wh_s,output_idx)
            
            Mass               = totalMass(wh)
            Energy             = totalEnergy(wu,wh,wh_s)
            AbsoluteVotricity  = totalAbsoluteVorticity(wu,wh)
            PotentialEnstrophy = totalPotentialEnstrophy(wu,wh)
            AngularMomentum     = totalAngularMomentum(wu,wh)
            
            print*,'Output time ',output_idx,'/',output_num,' dtn = ',dtn
            print*,' MCR  = ',(Mass-Mass0)/Mass0
            print*,' ECR  = ',(Energy-Energy0)/Energy0
            print*,' TAV  = ',AbsoluteVotricity
            print*,' PECR = ',(PotentialEnstrophy-PotentialEnstrophy0)/PotentialEnstrophy0
            print*,' AMCR = ',(AngularMomentum-AngularMomentum0)/AngularMomentum0
        endif
    enddo
    
    call SYSTEM_CLOCK(time_2)
    print*,'It took ',dble(time_2-time_1)/10000.0,' seconds to run this program'
end program main
    