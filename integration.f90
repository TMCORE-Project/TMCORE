SUBROUTINE integration(dtn)
    use module_para , only : RKIND,config_dt, config_time_integration
    implicit none
    real(kind=RKIND) , intent(out) :: dtn
    
    if(trim(adjustl(config_time_integration))=='CRK4')then
        call CRK4(config_dt,dtn)
    elseif(trim(adjustl(config_time_integration))=='CRK3')then
        call CRK3(config_dt,dtn)
    elseif(trim(adjustl(config_time_integration))=='Gill4')then
        call Gill4(config_dt,dtn)
    elseif(trim(adjustl(config_time_integration))=='PC3')then
        call PC3(config_dt,dtn)
    else
        stop 'Unknown integration scheme'
    endif
    
    !print*,dtn
    
END SUBROUTINE integration
    
