SUBROUTINE CRK4(dt,dtn)
    use module_para      , only: RKIND,config_energy_conservation,config_energy_conservation_scheme
    use module_array     , only: u,wu,wh,wh_s
    use mesh_info        , only: nEdges,nCells
    use diagnostics_tools, only: inner_product
    implicit none
    
    real(kind=RKIND),intent(in ) :: dt
    real(kind=RKIND),intent(out) :: dtn
    real(kind=RKIND),dimension(nEdges) :: wu1,wu2,wu3,wu4
    real(kind=RKIND),dimension(nCells) :: wh1,wh2,wh3,wh4
    real(kind=RKIND),dimension(nEdges) :: Lwu1,Lwu2,Lwu3,Lwu4
    real(kind=RKIND),dimension(nCells) :: Lwh1,Lwh2,Lwh3,Lwh4
    real(kind=RKIND),dimension(nEdges) :: phi4_wu
    real(kind=RKIND),dimension(nCells) :: phi4_wh
    
    real(kind=RKIND)                   :: phi4_norm2,R1R2,R2R3,R3R4
    real(kind=RKIND)                   :: beta_n
    
    wu1 = wu
    wh1 = wh
    
    call L_operator(wu1,wh1,wh_s,Lwu1,Lwh1)
    
    wu2 = wu - 0.5d0*dt*Lwu1
    wh2 = wh - 0.5d0*dt*Lwh1
    
    call L_operator(wu2,wh2,wh_s,Lwu2,Lwh2)
    
    wu3 = wu - 0.5d0*dt*Lwu2
    wh3 = wh - 0.5d0*dt*Lwh2
    
    call L_operator(wu3,wh3,wh_s,Lwu3,Lwh3)
    
    wu4 = wu - dt*Lwu3
    wh4 = wh - dt*Lwh3
    
    call L_operator(wu4,wh4,wh_s,Lwu4,Lwh4)
    
    phi4_wu = (Lwu1+2.d0*Lwu2+2.d0*Lwu3+Lwu4)/6.d0
    phi4_wh = (Lwh1+2.d0*Lwh2+2.d0*Lwh3+Lwh4)/6.d0
    
    if(config_energy_conservation)then
        if(config_energy_conservation_scheme==1)then
            dtn        = 2.d0*inner_product(phi4_wu,phi4_wh,wu,wh)/inner_product(phi4_wu,phi4_wh,phi4_wu,phi4_wh)
        elseif(config_energy_conservation_scheme==2)then
            phi4_norm2 = inner_product(phi4_wu,phi4_wh,phi4_wu,phi4_wh)
            R1R2       = inner_product(Lwu1   ,Lwh1   ,Lwu2   ,Lwh2   )
            R2R3       = inner_product(Lwu2   ,Lwh2   ,Lwu3   ,Lwh3   )
            R3R4       = inner_product(Lwu3   ,Lwh3   ,Lwu4   ,Lwh4   )
            
            beta_n     = (R1R2 + R2R3 + R3R4)/(3.d0*phi4_norm2)
            dtn        = beta_n*dt
        endif
    else
        dtn        = dt
    endif
    !print*,beta_n*dt-2.d0*inner_product(phi4_wu,phi4_wh,wu,wh)/phi4_norm2
    
    wu = wu - dtn*phi4_wu
    wh = wh - dtn*phi4_wh
    
    ! Accuracy Check
    !print*,dt *inner_product(phi4_wu,phi4_wh,phi4_wu,phi4_wh)+2.d0*inner_product(phi4_wu,phi4_wh,wu,wh)
    !print*,dtn*inner_product(phi4_wu,phi4_wh,phi4_wu,phi4_wh)+2.d0*inner_product(phi4_wu,phi4_wh,wu,wh)
    !print*,( dtn*inner_product(phi4_wu,phi4_wh,phi4_wu,phi4_wh)+2.d0*inner_product(phi4_wu,phi4_wh,wu,wh)  &
    !       -(dt *inner_product(phi4_wu,phi4_wh,phi4_wu,phi4_wh)+2.d0*inner_product(phi4_wu,phi4_wh,wu,wh)))&
    !       /inner_product(phi4_wu,phi4_wh,phi4_wu,phi4_wh) + dt
    !print*,dtn
    !print*,inner_product(phi4_wu,phi4_wh,phi4_wu,phi4_wh)
END SUBROUTINE CRK4