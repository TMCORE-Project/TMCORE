SUBROUTINE CRK3(dt,dtn)
    use module_para      , only: RKIND,config_energy_conservation,config_energy_conservation_scheme
    use module_array     , only: u,wu,wh,wh_s
    use mesh_info        , only: nEdges,nCells
    use diagnostics_tools, only: inner_product
    implicit none
    
    real(kind=RKIND),intent(in ) :: dt
    real(kind=RKIND),intent(out) :: dtn
    real(kind=RKIND),dimension(nEdges) :: wu1,wu2,wu3
    real(kind=RKIND),dimension(nCells) :: wh1,wh2,wh3
    real(kind=RKIND),dimension(nEdges) :: Lwu1,Lwu2,Lwu3
    real(kind=RKIND),dimension(nCells) :: Lwh1,Lwh2,Lwh3
    real(kind=RKIND),dimension(nEdges) :: phi3_wu
    real(kind=RKIND),dimension(nCells) :: phi3_wh
    
    real(kind=RKIND)                   :: phi3_norm2,R1R2,R1R3,R2R3
    real(kind=RKIND)                   :: beta_n
    
    wu1 = wu
    wh1 = wh
    
    call L_operator(wu1,wh1,wh_s,Lwu1,Lwh1)
    
    wu2 = wu - 0.5d0*dt*Lwu1
    wh2 = wh - 0.5d0*dt*Lwh1
    
    call L_operator(wu2,wh2,wh_s,Lwu2,Lwh2)
    
    wu3 = wu + dt*Lwu1 - 2.d0*dt*Lwu2
    wh3 = wh + dt*Lwh1 - 2.d0*dt*Lwh2
    
    call L_operator(wu3,wh3,wh_s,Lwu3,Lwh3)
    
    phi3_wu = (Lwu1+4.d0*Lwu2+Lwu3)/6.d0
    phi3_wh = (Lwh1+4.d0*Lwh2+Lwh3)/6.d0
    
    if(config_energy_conservation)then
        if(config_energy_conservation_scheme==1)then
            dtn        = 2.d0*inner_product(phi3_wu,phi3_wh,wu,wh)/inner_product(phi3_wu,phi3_wh,phi3_wu,phi3_wh)
        elseif(config_energy_conservation_scheme==2)then
            phi3_norm2 = inner_product(phi3_wu,phi3_wh,phi3_wu,phi3_wh)
            R1R2       = inner_product(Lwu1   ,Lwh1   ,Lwu2   ,Lwh2   )
            R1R3       = inner_product(Lwu1   ,Lwh1   ,Lwu3   ,Lwh3   )
            R2R3       = inner_product(Lwu2   ,Lwh2   ,Lwu3   ,Lwh3   )
            
            beta_n     = (2.d0*R1R2 - R1R3 + 2.d0*R2R3)/(3.d0*phi3_norm2)
            dtn        = beta_n*dt
        endif
    else
        dtn        = dt
    endif
    !print*,beta_n*dt-2.d0*inner_product(phi3_wu,phi3_wh,wu,wh)/phi3_norm2
    
    wu = wu - dtn*phi3_wu
    wh = wh - dtn*phi3_wh
    
    ! Accuracy Check
    !print*,dt *inner_product(phi3_wu,phi3_wh,phi3_wu,phi3_wh)+2.d0*inner_product(phi3_wu,phi3_wh,wu,wh)
    !print*,dtn*inner_product(phi3_wu,phi3_wh,phi3_wu,phi3_wh)+2.d0*inner_product(phi3_wu,phi3_wh,wu,wh)
    !print*,( dtn*inner_product(phi3_wu,phi3_wh,phi3_wu,phi3_wh)+2.d0*inner_product(phi3_wu,phi3_wh,wu,wh)  &
    !       -(dt *inner_product(phi3_wu,phi3_wh,phi3_wu,phi3_wh)+2.d0*inner_product(phi3_wu,phi3_wh,wu,wh)))&
    !       /inner_product(phi3_wu,phi3_wh,phi3_wu,phi3_wh) + dt
    !print*,dtn
    !print*,inner_product(phi3_wu,phi3_wh,phi3_wu,phi3_wh)
END SUBROUTINE CRK3