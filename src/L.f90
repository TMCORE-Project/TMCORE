SUBROUTINE L_operator(wu,wh,wh_s,Lwu,Lwh)
    use module_para      , only : RKIND,config_enstrophy_conservation
    use mesh_info        , only : nCells      ,&
                                  nEdges      ,&
                                  nVertices   ,&
                                  areaCell    ,&
                                  areaEdge    ,&
                                  areaTriangle,&
                                  cellsOnEdge ,&
                                  edgesOnCell ,&
                                  nEdgesOnCell,&
                                  dvEdge      ,&
                                  dcEdge      ,&
                                  sum_areaEdge,&
                                  sum_areaCell,&
                                  n
    use diagnostics_tools, only : IIAP                            ,&
                                  interp_C2E                      ,&
                                  interp_C2V                      ,&
                                  inner_product                   ,&
                                  compute_vorticity               ,&
                                  compute_divergence              ,&
                                  compute_ke                      ,&
                                  compute_v                       ,&
                                  compute_pv_vertex               ,&
                                  compute_pv_cell                 ,&
                                  compute_pv_cell_v               ,&
                                  compute_pv_edge                 ,&
                                  compute_Q                       ,&
                                  compute_tend_vorticity_by_tend_u,&
                                  potentialEnstrophy_Conservation
    implicit none
    
    real(kind=RKIND), intent(in ) :: wu  (nEdges)
    real(kind=RKIND), intent(in ) :: wh  (nCells)
    real(kind=RKIND), intent(in ) :: wh_s(nCells)
    real(kind=RKIND), intent(out) :: Lwu (nEdges)
    real(kind=RKIND), intent(out) :: Lwh (nCells)
    
    real(kind=RKIND)              :: wh_edge   (nEdges   )  ! Working array
    real(kind=RKIND)              :: wh_vertex (nVertices)  ! Working array
    
    real(kind=RKIND)              :: u         (nEdges   )  ! Working array
    real(kind=RKIND)              :: v         (nEdges   )  ! Working array
    real(kind=RKIND)              :: vorticity (nVertices)  ! Working array
    real(kind=RKIND)              :: ke        (nCells   )  ! Working array
    real(kind=RKIND)              :: pv_vertex (nVertices)  ! Working array
    real(kind=RKIND)              :: pv_cell   (nCells   )  ! Working array
    real(kind=RKIND)              :: pv_edge   (nEdges   )  ! Working array
    real(kind=RKIND)              :: Q         (nEdges   )  ! Working array
    real(kind=RKIND)              :: h_edge    (nEdges   )  ! Working array
    
    real(kind=RKIND)              :: tend_u        (nEdges   )  ! Working array
    real(kind=RKIND)              :: tend_wu       (nEdges   )  ! Working array
    real(kind=RKIND)              :: tend_wh       (nCells   )  ! Working array
    real(kind=RKIND)              :: tend_wh_edge  (nEdges   )  ! Working array
    real(kind=RKIND)              :: tend_wh_vertex(nVertices)  ! Working array
    real(kind=RKIND)              :: tend_vorticity(nVertices)  ! Working array
    
    real(kind=RKIND)              :: dkedx       (nEdges   )  ! Working array
    real(kind=RKIND)              :: dwhdx       (nEdges   )  ! Working array
    
    real(kind=RKIND)              :: flux        (nEdges   )  ! Working array
    
    real(kind=RKIND)              :: epsilon
    
    integer iCell
    integer cell1(nEdges),cell2(nEdges)
    
    call IIAP              (wu        ,wh        ,u                            )
    call interp_C2E        (wh        ,wh_edge   ,2                            )
    call compute_ke        (u         ,ke                                      )
        
    h_edge = dsqrt(wh_edge)
    
    cell1 = cellsOnEdge(1,:)
    cell2 = cellsOnEdge(2,:)
    
    !
    ! Compute height tendency for each cell
    !
    tend_wh = 0.d0
    flux    = u * dvEdge * wh_edge
    do iCell = 1, nCells
       tend_wh(iCell) = -sum(n               (1:nEdgesOnCell(iCell),iCell) &
                            *flux(edgesOnCell(1:nEdgesOnCell(iCell),iCell)))
    end do 
    
    tend_wh = tend_wh / areaCell
    
    call interp_C2E(tend_wh,tend_wh_edge,2)
    call interp_C2V(tend_wh,tend_wh_vertex)
    call interp_C2V(wh     ,wh_vertex     )
    
    !
    ! Compute u (normal) velocity tendency for each edge (cell face)
    !
    
    tend_u = 0.d0   
    dkedx  = (ke(cell2)               - ke(cell1)              )/dcEdge
    dwhdx  = (wh(cell2) + wh_s(cell2) - wh(cell1) - wh_s(cell1))/dcEdge
    
    call compute_vorticity (u         ,vorticity                                                        )
    call compute_v         (u         ,v                                                                )
    call compute_pv_vertex (u         ,wh_vertex ,vorticity ,pv_vertex                                  )
    call compute_pv_cell   (pv_vertex ,pv_cell                                                          )
    !call compute_pv_cell_v (wh        ,v         ,pv_cell                                               )
    call compute_pv_edge   (u         ,v         ,wh_edge   ,tend_wh_vertex ,pv_vertex ,pv_cell ,pv_edge)
    call compute_Q         (u         ,wh_edge   ,pv_edge   ,Q                                          )
    
    if(config_enstrophy_conservation)then
        call potentialEnstrophy_Conservation(dkedx,dwhdx,Q,pv_vertex,tend_wh,epsilon)
    else
        epsilon = 1.d0
    endif
    
    tend_u = Q*epsilon - dkedx - dwhdx
    
    ! Compute U (normal) velocity tendency for each edge (cell face), used by square conservation        
    tend_wu = h_edge*tend_u + 0.5d0*u/h_edge*tend_wh_edge
    
    Lwu = -tend_wu
    Lwh = -tend_wh
    
    !! Mass conservation check
    !print*,sum(tend_wh*areaCell)/sum(areaCell)
    
    !! Coriolis force check
    !print*,sum(Q*u*wh_edge*areaEdge)/sum(areaCell)
    
    !! Antisymmtry Check
    !print*,sum(wu*h_edge*dkedx*areaEdge),sum(wu*0.5*u/h_edge*tend_wh_edge*areaEdge)
    !print*,sum(wu*h_edge*dwhdx*areaEdge),sum(wh*tend_wh*areaCell)
    !print*,inner_product(wu,wh,Lwu,Lwh)
    !print*,(sum(tend_u*u*wh_edge*areaEdge)+sum((ke+wh)*tend_wh*areaCell))/sum(areaCell)
    
    !call compute_tend_vorticity_by_tend_u(tend_u,tend_vorticity)
    !print*,sum((pv_vertex*tend_vorticity-0.5d0*pv_vertex**2*tend_wh_vertex)*areaTriangle)!/sum(areaTriangle)
    
END SUBROUTINE L_operator
