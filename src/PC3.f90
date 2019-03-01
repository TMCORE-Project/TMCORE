SUBROUTINE PC3(dt,dtn)
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
    
    real(kind=RKIND)                   :: PC3_norm2
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
        
    if(config_energy_conservation)then
        PC3_norm2 = inner_product(Lwu3,Lwh3,Lwu3,Lwh3)
        if(config_energy_conservation_scheme==1)then
            dtn = 2.d0*inner_product(wu1,wh1,Lwu3,Lwh3)/PC3_norm2
        elseif(config_energy_conservation_scheme==2)then
            beta_n = inner_product(Lwu2,Lwh2,Lwu3,Lwh3)/PC3_norm2
            dtn    = beta_n*dt
        endif
    else
        dtn        = dt
    endif
    !print*,beta_n*dt-2.d0*inner_product(PC3_wu,PC3_wh,wu,wh)/PC3_norm2
    
    wu = wu - dtn*Lwu3
    wh = wh - dtn*Lwh3
    
    ! Accuracy Check
    !print*,dt *inner_product(PC3_wu,PC3_wh,PC3_wu,PC3_wh)+2.d0*inner_product(PC3_wu,PC3_wh,wu,wh)
    !print*,dtn*inner_product(PC3_wu,PC3_wh,PC3_wu,PC3_wh)+2.d0*inner_product(PC3_wu,PC3_wh,wu,wh)
    !print*,( dtn*inner_product(PC3_wu,PC3_wh,PC3_wu,PC3_wh)+2.d0*inner_product(PC3_wu,PC3_wh,wu,wh)  &
    !       -(dt *inner_product(PC3_wu,PC3_wh,PC3_wu,PC3_wh)+2.d0*inner_product(PC3_wu,PC3_wh,wu,wh)))&
    !       /inner_product(PC3_wu,PC3_wh,PC3_wu,PC3_wh) + dt
    !print*,dtn
    !print*,inner_product(PC3_wu,PC3_wh,PC3_wu,PC3_wh)
END SUBROUTINE PC3