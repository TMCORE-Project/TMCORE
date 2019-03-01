MODULE output_netCDF
    use netcdf
    use module_para       , only: RKIND, g, a, omega                ,&
                                  config_output_file                ,&
                                  config_dt                         ,&
                                  config_energy_conservation        ,&
                                  config_energy_conservation_scheme ,&
                                  config_time_integration           ,&
                                  config_test_case                  ,&
                                  config_PV_scheme
    use mesh_info         , only: handle_err   ,&
                                  nCells       ,&
                                  nEdges       ,&
                                  nVertices    ,&
                                  lonCell      ,&
                                  latCell      ,&
                                  lonEdge      ,&
                                  latEdge      ,&
                                  lonVertex    ,&
                                  latVertex    ,&
                                  cellsOnEdge  ,&
                                  dvEdge       ,&
                                  dcEdge       ,&
                                  areaCell     ,&
                                  areaEdge     ,&
                                  areaTriangle ,&
                                  angleEdge
                                  
    use diagnostics_tools , only: IIAP              ,&
                                  compute_vorticity ,&
                                  compute_pv_vertex ,&
                                  interp_C2V
    
    implicit none
    include 'netcdf.inc'
    
    contains
    
    subroutine write_netCDF(wu,wh,wh_s,output_idx)
        implicit none
        integer                            ,intent(in):: output_idx
        real(kind=RKIND) ,dimension(nEdges),intent(in):: wu
        real(kind=RKIND) ,dimension(nCells),intent(in):: wh
        real(kind=RKIND) ,dimension(nCells),intent(in):: wh_s
        
        real(kind=RKIND) ,dimension(nEdges)           :: u
        real(kind=RKIND) ,dimension(nVertices)        :: vorticity
        real(kind=RKIND) ,dimension(nVertices)        :: wh_vertex
        real(kind=RKIND) ,dimension(nVertices)        :: pv_vertex
        
        real(kind=RKIND) ,dimension(nCells)           :: h
        real(kind=RKIND) ,dimension(nCells)           :: hs   ! hs : surface height, hs = wh_s/g
        real(kind=RKIND) ,dimension(nCells)           :: hphs ! hphs = h + hs
        
        character*200                                 :: nc_file
        
        integer       ncid,status
        integer       nCells_dimid    ,&
                      nEdges_dimid    ,&
                      nVertices_dimid ,&
                      TWO_dimid       ,&
                      output_num_dimid
        
        integer       lonCell_id      ,&
                      latCell_id      ,&
                      lonEdge_id      ,&
                      latEdge_id      ,&
                      lonVertex_id    ,&
                      latVertex_id    ,&
                      cellsOnEdge_id  ,&
                      areaCell_id     ,&
                      areaEdge_id     ,&
                      areaTriangle_id ,&
                      angleEdge_id    ,&
                      dvEdge_id       ,&
                      dcEdge_id       ,&
                      u_id            ,&
                      wu_id           ,&
                      h_id            ,&
                      hs_id           ,&
                      hphs_id         ,&
                      vorticity_id    ,&
                      wh_vertex_id    ,&
                      pv_vertex_id
        
        character*200 nCells_c,dt_c,EC_sign,EC_scheme_sign,case_num_c
        
        ! Set output file name
        if(config_output_file/='')then
            nc_file = config_output_file
        else
            write(nCells_c,'(i)')nCells
            write(dt_c    ,'(i)')int(config_dt)
            if(config_energy_conservation)then
                EC_sign = 'EC'
                write(EC_scheme_sign,'(i)') config_energy_conservation_scheme
            else
                EC_sign        = 'NEC'
                EC_scheme_sign = ''
            endif
            
            write(case_num_c,'(i)')config_test_case
            
            nc_file = 'output_'//'case'//trim(adjustl(case_num_c))//'_'//trim(adjustl(nCells_c))//'_'//trim(adjustl(dt_c))//'s_'//trim(adjustl(EC_sign))//trim(adjustl(EC_scheme_sign))&
                      //'_'//trim(adjustl(config_time_integration))//'_'//trim(adjustl(config_PV_scheme))//'.nc'
        endif
                
        ! Diagnostic fields
        call IIAP             (wu, wh       , u                    )
        call compute_vorticity(u , vorticity                       )
        call interp_C2V       (wh, wh_vertex                       )
        call compute_pv_vertex(u , wh_vertex, vorticity ,pv_vertex )
        
        h    = wh/g
        hs   = wh_s/g
        hphs = h + hs
        
        if(output_idx==0)then
            
            !print*,'nf90_create'
            status = nf90_create(path=trim(adjustl(nc_file)),cmode=NF90_CLOBBER,ncid=ncid)
            !if(status/=nf90_noerr) call handle_err(status)
            
            status = nf90_put_att(ncid,NF90_GLOBAL,'g'    ,g    )
            status = nf90_put_att(ncid,NF90_GLOBAL,'a'    ,a    )
            status = nf90_put_att(ncid,NF90_GLOBAL,'omega',omega)
            
            !print*,'nf90_def_dim'
            status = nf90_def_dim(ncid,'nCells'   ,nCells             ,nCells_dimid    )
            status = nf90_def_dim(ncid,'nEdges'   ,nEdges             ,nEdges_dimid    )
            status = nf90_def_dim(ncid,'nVertices',nVertices          ,nVertices_dimid )
            status = nf90_def_dim(ncid,'TWO'      ,2                  ,TWO_dimid       )
            status = nf90_def_dim(ncid,'nt'       ,NF90_UNLIMITED     ,output_num_dimid)
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_def_var'
            status = nf90_def_var(ncid,'lonCell'     ,NF90_DOUBLE,(/nCells_dimid                    /),lonCell_id     )
            status = nf90_def_var(ncid,'latCell'     ,NF90_DOUBLE,(/nCells_dimid                    /),latCell_id     )
            status = nf90_def_var(ncid,'lonEdge'     ,NF90_DOUBLE,(/nEdges_dimid                    /),lonEdge_id     )
            status = nf90_def_var(ncid,'latEdge'     ,NF90_DOUBLE,(/nEdges_dimid                    /),latEdge_id     )
            status = nf90_def_var(ncid,'lonVertex'   ,NF90_DOUBLE,(/nVertices_dimid                 /),lonVertex_id   )
            status = nf90_def_var(ncid,'latVertex'   ,NF90_DOUBLE,(/nVertices_dimid                 /),latVertex_id   )
            status = nf90_def_var(ncid,'areaCell'    ,NF90_DOUBLE,(/nCells_dimid                    /),areaCell_id    )
            status = nf90_def_var(ncid,'areaEdge'    ,NF90_DOUBLE,(/nEdges_dimid                    /),areaEdge_id    )
            status = nf90_def_var(ncid,'areaTriangle',NF90_DOUBLE,(/nVertices_dimid                 /),areaTriangle_id)
            status = nf90_def_var(ncid,'areaEdge'    ,NF90_DOUBLE,(/nEdges_dimid                    /),areaEdge_id    )
            status = nf90_def_var(ncid,'dvEdge'      ,NF90_DOUBLE,(/nEdges_dimid                    /),dvEdge_id      )
            status = nf90_def_var(ncid,'dcEdge'      ,NF90_DOUBLE,(/nEdges_dimid                    /),dcEdge_id      )
            status = nf90_def_var(ncid,'cellsOnEdge' ,NF90_DOUBLE,(/TWO_dimid      ,nEdges_dimid    /),cellsOnEdge_id )
            status = nf90_def_var(ncid,'u'           ,NF90_DOUBLE,(/nEdges_dimid   ,output_num_dimid/),u_id           )
            status = nf90_def_var(ncid,'U'           ,NF90_DOUBLE,(/nEdges_dimid   ,output_num_dimid/),wu_id          )
            status = nf90_def_var(ncid,'h'           ,NF90_DOUBLE,(/nCells_dimid   ,output_num_dimid/),h_id           )
            status = nf90_def_var(ncid,'hs'          ,NF90_DOUBLE,(/nCells_dimid   ,output_num_dimid/),hs_id          )
            status = nf90_def_var(ncid,'hphs'        ,NF90_DOUBLE,(/nCells_dimid   ,output_num_dimid/),hphs_id        )
            status = nf90_def_var(ncid,'vorticity'   ,NF90_DOUBLE,(/nVertices_dimid,output_num_dimid/),vorticity_id   )
            status = nf90_def_var(ncid,'wh_vertex'   ,NF90_DOUBLE,(/nVertices_dimid,output_num_dimid/),wh_vertex_id   )
            status = nf90_def_var(ncid,'pv_vertex'   ,NF90_DOUBLE,(/nVertices_dimid,output_num_dimid/),pv_vertex_id   )
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_put_att'
            !status = nf90_put_att()
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_enddef'
            status = nf90_enddef(ncid)
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_put_var'
            status = nf90_put_var(ncid,lonCell_id      ,lonCell                                          )
            status = nf90_put_var(ncid,latCell_id      ,latCEll                                          )
            status = nf90_put_var(ncid,lonEdge_id      ,lonEdge                                          )
            status = nf90_put_var(ncid,latEdge_id      ,latEdge                                          )
            status = nf90_put_var(ncid,lonVertex_id    ,lonVertex                                        )
            status = nf90_put_var(ncid,latVertex_id    ,latVertex                                        )
            status = nf90_put_var(ncid,areaCell_id     ,areaCell                                         )
            status = nf90_put_var(ncid,areaEdge_id     ,areaEdge                                         )
            status = nf90_put_var(ncid,areaTriangle_id ,areaTriangle                                     )
            status = nf90_put_var(ncid,areaEdge_id     ,areaEdge                                         )
            status = nf90_put_var(ncid,dvEdge_id       ,dvEdge                                           )
            status = nf90_put_var(ncid,dcEdge_id       ,dcEdge                                           )
            status = nf90_put_var(ncid,cellsOnEdge_id  ,cellsOnEdge                                      )
            status = nf90_put_var(ncid,u_id            ,u            ,start=(/1,1/),count=(/nEdges   ,1/))
            status = nf90_put_var(ncid,wu_id           ,wu           ,start=(/1,1/),count=(/nEdges   ,1/))
            status = nf90_put_var(ncid,h_id            ,h            ,start=(/1,1/),count=(/nCells   ,1/))
            status = nf90_put_var(ncid,hs_id           ,hs           ,start=(/1,1/),count=(/nCells   ,1/))
            status = nf90_put_var(ncid,hphs_id         ,hphs         ,start=(/1,1/),count=(/nCells   ,1/))
            status = nf90_put_var(ncid,vorticity_id    ,vorticity    ,start=(/1,1/),count=(/nVertices,1/))
            status = nf90_put_var(ncid,wh_vertex_id    ,wh_vertex    ,start=(/1,1/),count=(/nVertices,1/))
            status = nf90_put_var(ncid,pv_vertex_id    ,pv_vertex    ,start=(/1,1/),count=(/nVertices,1/))
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_close'
            status = nf90_close(ncid)
            !if(status/=nf90_noerr) call handle_err(status)
        else
            !print*,'nf90_open'
            status = nf90_open(trim(adjustl(nc_file)),NF90_WRITE,ncid)
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_inq_varid'
            status = nf90_inq_varid(ncid,'u'         ,u_id         )
            status = nf90_inq_varid(ncid,'U'         ,wu_id        )
            status = nf90_inq_varid(ncid,'h'         ,h_id         )
            status = nf90_inq_varid(ncid,'hs'        ,hs_id        )
            status = nf90_inq_varid(ncid,'hphs'      ,hphs_id      )
            status = nf90_inq_varid(ncid,'vorticity' ,vorticity_id )
            status = nf90_inq_varid(ncid,'wh_vertex' ,wh_vertex_id )
            status = nf90_inq_varid(ncid,'pv_vertex' ,pv_vertex_id )
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_put_var'
            status = nf90_put_var(ncid,u_id         ,u        ,start=(/1,output_idx+1/),count=(/nEdges   ,1/))
            status = nf90_put_var(ncid,wu_id        ,wu       ,start=(/1,output_idx+1/),count=(/nEdges   ,1/))
            status = nf90_put_var(ncid,h_id         ,h        ,start=(/1,output_idx+1/),count=(/nCells   ,1/))
            status = nf90_put_var(ncid,hs_id        ,hs       ,start=(/1,output_idx+1/),count=(/nCells   ,1/))
            status = nf90_put_var(ncid,hphs_id      ,hphs     ,start=(/1,output_idx+1/),count=(/nCells   ,1/))
            status = nf90_put_var(ncid,vorticity_id ,vorticity,start=(/1,output_idx+1/),count=(/nVertices,1/))
            status = nf90_put_var(ncid,wh_vertex_id ,wh_vertex,start=(/1,output_idx+1/),count=(/nVertices,1/))
            status = nf90_put_var(ncid,pv_vertex_id ,pv_vertex,start=(/1,output_idx+1/),count=(/nVertices,1/))
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_put_att'
            !status = nf90_put_att
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_close'
            status = nf90_close(ncid)
            !if(status/=nf90_noerr) call handle_err(status)
        endif
        
    end subroutine write_netCDF
    
END MODULE output_netCDF