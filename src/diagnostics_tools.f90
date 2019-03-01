module diagnostics_tools
    use module_para,only: RKIND,IKIND           ,&
                          a                     ,&
                          pi                    ,&
                          omega                 ,&
                          config_PV_scheme      ,&
                          config_apvm_upwinding ,&
                          config_CLUST_upwinding,&
                          config_dt
    use mesh_info  ,only: nCells              ,&
                          nEdges              ,&
                          nVertices           ,&
                          vertexDegree        ,&
                          maxEdges            ,&
                          maxEdges2           ,&
                          latCell,lonCell     ,&
                          latEdge,lonEdge     ,&
                          latVertex,lonVertex ,&
                          dvEdge              ,&
                          dv1Edge             ,&
                          dv2Edge             ,&
                          dcEdge              ,&
                          angleEdge           ,&
                          areaCell            ,&
                          areaTriangle        ,&
                          areaEdge            ,&
                          nEdgesOnCell        ,&
                          nEdgesOnEdge        ,&
                          cellsOnEdge         ,&
                          edgesOnCell         ,&
                          edgesOnEdge         ,&
                          weightsOnEdge       ,&
                          cellsOnCell         ,&
                          verticesOnCell      ,&
                          verticesOnEdge      ,&
                          edgesOnVertex       ,&
                          cellsOnVertex       ,&
                          kiteAreasOnVertex   ,&
                          fCell               ,&
                          fEdge               ,&
                          fVertex             ,&
                          sum_areaEdge        ,&
                          sum_areaCell        ,&
                          cos_angleVertex12   ,&
                          cos_angleVertex13   ,&
                          cos_angleVertex23   ,&
                          t                   ,&
                          n
    implicit none
    
    contains
    
    !
    ! IAP transformation
    !
    subroutine IAP(u,wh,wu)
        implicit none
        real(kind=RKIND), intent(in )           :: u      (nEdges)
        real(kind=RKIND), intent(in )           :: wh     (nCells)
        real(kind=RKIND), intent(out)           :: wu     (nEdges)
        
        real(kind=RKIND)                        :: wh_edge(nEdges) ! Working array
        
        call interp_C2E(wh,wh_edge,2)
        wu = dsqrt(wh_edge)*u
    end subroutine IAP
    
    !
    ! Inverse IAP transformation
    !
    subroutine IIAP(wu,wh,u)
        implicit none
        real(kind=RKIND), intent(in ) :: wu(nEdges)
        real(kind=RKIND), intent(in ) :: wh(nCells)
        real(kind=RKIND), intent(out) :: u (nEdges)
        
        real(kind=RKIND)              :: wh_edge(nEdges) ! Working array
        
        call interp_C2E(wh,wh_edge,2)
        u = wu/dsqrt(wh_edge)
    end subroutine IIAP
    
    !
    ! Interpolate the fields on cells to edges
    !
    subroutine interp_C2E(f_cell,f_edge,order_in)
        implicit none
        real(kind=RKIND), intent(in )          :: f_cell(nCells)
        integer         , intent(in ),optional :: order_in
        real(kind=RKIND), intent(out)          :: f_edge(nEdges)
        integer                                :: order
        
        integer cell1,cell2
        
        integer iEdge
        
        order = 2
        if(present(order_in))order = order_in
        
        if(order==2)then
            do iEdge = 1, nEdges
                cell1 = cellsOnEdge(1,iEdge)
                cell2 = cellsOnEdge(2,iEdge)
                
                f_edge(iEdge) = 0.5d0 * (f_cell(cell1) + f_cell(cell2))
            end do
        else
            stop 'Now we support only interpolating fields from cell to edge with 2nd order accuacy'
        endif
            
    end subroutine interp_C2E

    !
    ! Interpolate the fields on cells to vertices
    !
    subroutine interp_C2V(f_cell,f_vertex)
        implicit none
        real(kind=RKIND), intent(in ) :: f_cell  (nCells   )
        real(kind=RKIND), intent(out) :: f_vertex(nVertices)
        
        integer i,iVertex
        
        do iVertex = 1,nVertices
            f_vertex(iVertex) = 0.d0
            do i = 1, vertexDegree
                f_vertex(iVertex) = f_vertex(iVertex) + f_cell(cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
            end do
        end do
        f_vertex = f_vertex / areaTriangle
    end subroutine interp_C2V

    !
    ! Interpolate the fields on vertices to cells
    !
    subroutine interp_V2C(f_vertex,f_cell)
        implicit none
        real(kind=RKIND),intent(in ) :: f_vertex(nVertices)
        real(kind=RKIND),intent(out) :: f_cell  (nCells)
        
        integer :: i,iCell,iVertex
        
        f_cell = 0.d0
        do iVertex = 1, nVertices
            do i = 1, vertexDegree
                iCell = cellsOnVertex(i,iVertex)
                
                f_cell(iCell) = f_cell(iCell) + kiteAreasOnVertex(i, iVertex) * f_vertex(iVertex) / areaCell(iCell)
            enddo
        enddo
    end subroutine interp_V2C
    
    !
    ! Compute vorticity
    !
    subroutine compute_vorticity(u,vorticity)
        implicit none
        real(kind=RKIND), intent(in ) :: u          (nEdges)
        real(kind=RKIND), intent(out) :: vorticity  (nVertices)
        real(kind=RKIND)              :: circulation(nVertices)
        
        integer vertex1,vertex2
        integer iEdge,iVertex
        
        circulation = 0.d0
        do iEdge = 1, nEdges
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)
            
            circulation(vertex1) = circulation(vertex1) - dcEdge(iEdge) * u(iEdge)
            circulation(vertex2) = circulation(vertex2) + dcEdge(iEdge) * u(iEdge)
        end do
        
        vorticity = circulation / areaTriangle
    end subroutine compute_vorticity
    
    !
    ! Compute divergence
    !
    subroutine compute_divergence(u,divergence)
        implicit none
        real(kind=RKIND),intent(in ) :: u         (nEdges)
        real(kind=RKIND),intent(out) :: divergence(nCells)
        
        integer cell1,cell2
        integer iEdge
        
        divergence = 0.d0
        do iEdge = 1, nEdges
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           
           divergence(cell1) = divergence(cell1) + u(iEdge)*dvEdge(iEdge)
           divergence(cell2) = divergence(cell2) - u(iEdge)*dvEdge(iEdge)
        end do
        
        divergence = divergence / areaCell
        
    end subroutine compute_divergence
    
    !
    ! Compute kinetic energy in each cell
    !
    subroutine compute_ke(u,ke)
        implicit none
        real(kind=RKIND),intent(in ) :: u (nEdges)
        real(kind=RKIND),intent(out) :: ke(nCells)
        
        integer i,iCell,iEdge
        
        ke = 0.d0
        do iCell = 1, nCells
            do i = 1, nEdgesOnCell(iCell)
                iEdge = edgesOnCell(i,iCell)
                ke(iCell) = ke(iCell) + 0.25d0 * areaEdge(iEdge) * u(iEdge)**2.d0
            end do
        end do
        
        ke= ke / areaCell
    end subroutine compute_ke
        
    !
    ! Compute v (tangential) velocities
    !
    subroutine compute_v(u,v)
        implicit none
        real(kind=RKIND),intent(in ) :: u(nEdges)
        real(kind=RKIND),intent(out) :: v(nEdges)
        
        integer eoe
        integer i,iEdge
        
        v = 0.d0
        do iEdge = 1,nEdges
            do i = 1, nEdgesOnEdge(iEdge)
                eoe = edgesOnEdge(i,iEdge)
                v(iEdge) = v(iEdge) + weightsOnEdge(i,iEdge) * u(eoe)
            end do
        end do
    end subroutine compute_v
    
    !
    ! Compute height at vertices, pv at vertices, and average pv to edge locations
    !  ( this computes pv_vertex at all vertices bounding real cells and distance-1 ghost cells )
    !
    subroutine compute_pv_vertex(u,wh_vertex,vorticity,pv_vertex)
        implicit none
        real(kind=RKIND),intent(in ) :: u        (nEdges   )
        real(kind=RKIND),intent(in ) :: wh_vertex(nVertices)
        real(kind=RKIND),intent(in ) :: vorticity(nVertices)
        real(kind=RKIND),intent(out) :: pv_vertex(nVertices)
        
        pv_vertex = (fVertex + vorticity) / wh_vertex
    end subroutine compute_pv_vertex

    !
    ! Compute pv at cell centers
    !    ( this computes pv_cell for all real cells and distance-1 ghost cells )
    !
    subroutine compute_pv_cell(pv_vertex,pv_cell)
        implicit none
        real(kind=RKIND),intent(in ) :: pv_vertex(nVertices)
        real(kind=RKIND),intent(out) :: pv_cell  (nCells)
        
        call interp_V2C(pv_vertex,pv_cell)
    end subroutine compute_pv_cell
    
    subroutine compute_pv_cell_v(wh,v,pv_cell_v)
        implicit none
        real(kind=RKIND),intent(in ) :: wh(nCells)
        real(kind=RKIND),intent(in ) :: v (nEdges)
        real(kind=RKIND),intent(out) :: pv_cell_v(nCells)
        
        integer(kind=IKIND)          :: iCell,iEdge
        
        pv_cell_v = 0.d0
        do iCell = 1,nCells
            do iEdge = 1,nEdgesOnCell(iCell)
                pv_cell_v(iCell) = pv_cell_v(iCell) + v(iEdge)*dvEdge(iEdge)
            enddo
        enddo
        pv_cell_v = (pv_cell_v/wh+fCell)/areaCell
    endsubroutine compute_pv_cell_v
    
    ! Compute q for each vertex   
    subroutine compute_Q(u,wh_edge,pv_edge,Q)
        implicit none
        real(kind=RKIND),intent(in ) :: u      (nEdges)
        real(kind=RKIND),intent(in ) :: wh_edge(nEdges)
        real(kind=RKIND),intent(in ) :: pv_edge(nEdges)
        real(kind=RKIND),intent(out) :: Q      (nEdges)
        
        !real(kind=RKIND)             :: v(nEdges)
        real(kind=RKIND)             :: workpv
        
        integer eoe
        integer i,j,iEdge
    
        !Energy Conservation Scheme
        do iEdge = 1, nEdges
            Q(iEdge) = 0.d0
            do j = 1, nEdgesOnEdge(iEdge)
                eoe      = edgesOnEdge(j,iEdge)
                workpv   = 0.5d0*(pv_edge(iEdge) + pv_edge(eoe))
                
                Q(iEdge) = Q(iEdge) + weightsOnEdge(j,iEdge)*u(eoe)*workpv*wh_edge(eoe)
            end do
        end do
        
        ! Potential Enstrophy Conservation Scheme
        ! call compute_v(u,v)
        ! Q = v*pv_edge*wh_edge
    end subroutine compute_Q
    
    subroutine compute_tend_vorticity_by_tend_u(tend_u,tend_vorticity)
        implicit none
        real(kind=RKIND),intent(in ) :: tend_u        (nEdges   )
        real(kind=RKIND),intent(out) :: tend_vorticity(nVertices)
        
        real(kind=RKIND)             :: temp
        
        integer(kind=IKIND)          :: vertex1,vertex2
        integer(kind=IKIND)          :: i,iEdge
        
        ! Compute vorticity tendency for each vertex
        tend_vorticity = 0.d0
        do iEdge = 1, nEdges
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)
            
            temp    = tend_u(iEdge)*dcEdge(iEdge)
            
            tend_vorticity(vertex1) = tend_vorticity(vertex1) - temp
            tend_vorticity(vertex2) = tend_vorticity(vertex2) + temp
        end do            
        tend_vorticity = tend_vorticity / areaTriangle
        
    end subroutine compute_tend_vorticity_by_tend_u
    
    !
    ! Compute pv at the edges
    !   ( this computes pv_edge at all edges bounding real cells )
    !
    subroutine compute_pv_edge(u,v,wh_edge,tend_wh_vertex,pv_vertex,pv_cell,pv_edge)
        implicit none
        real(kind=RKIND),intent(in ) :: u             (nEdges   )
        real(kind=RKIND),intent(in ) :: v             (nEdges   )
        real(kind=RKIND),intent(in ) :: wh_edge       (nEdges   )
        real(kind=RKIND),intent(in ) :: tend_wh_vertex(nVertices)
        real(kind=RKIND),intent(in ) :: pv_vertex     (nVertices)
        real(kind=RKIND),intent(in ) :: pv_cell       (nCells   )
        real(kind=RKIND),intent(out) :: pv_edge       (nEdges   )
        
        if(trim(adjustl(config_PV_scheme))=='APVM')then
            call compute_pv_edge_APVM(u,v,pv_vertex,pv_cell,pv_edge)
        elseif(trim(adjustl(config_PV_scheme))=='APVM_Conservation')then
            call compute_pv_edge_APVM_Conservation(u,v,wh_edge,tend_wh_vertex,pv_vertex,pv_cell,pv_edge)
        elseif(trim(adjustl(config_PV_scheme))=='order2')then
            call compute_pv_edge_order2(pv_vertex, pv_edge)
        elseif(trim(adjustl(config_PV_scheme))=='order2_smooth')then
            call compute_pv_edge_order2_smooth(pv_vertex, pv_cell, pv_edge)
        elseif(trim(adjustl(config_PV_scheme))=='LUST')then
            call compute_pv_edge_LUST(u, v, pv_vertex, pv_edge)
        elseif(trim(adjustl(config_PV_scheme))=='CLUST')then
            call compute_pv_edge_CLUST(u, v, pv_vertex, pv_edge)
        else
            stop 'Unknow PV scheme, please choose from APVM, order2, LUST and CLUST'
        endif
        
    end subroutine compute_pv_edge
    
    !
    ! Compute pv at the edges
    !   ( this computes pv_edge at all edges bounding real cells )
    !
    subroutine compute_pv_edge_APVM(u,v,pv_vertex,pv_cell,pv_edge)
        implicit none
        real(kind=RKIND),intent(in ) :: u        (nEdges   )
        real(kind=RKIND),intent(in ) :: v        (nEdges   )
        real(kind=RKIND),intent(in ) :: pv_vertex(nVertices)
        real(kind=RKIND),intent(in ) :: pv_cell  (nCells   )
        real(kind=RKIND),intent(out) :: pv_edge  (nEdges   )
        
        real(kind=RKIND)             :: gradPVt(nEdges)
        real(kind=RKIND)             :: gradPVn(nEdges)
        
        real(kind=RKIND)             :: lambda    (nEdges   )
        integer(kind=IKIND)          :: vertex1   (nEdges   ),&
                                        vertex2   (nEdges   ),&
                                        cell1     (nEdges   ),&
                                        cell2     (nEdges   )
        integer i,iEdge,iVertex
        
        vertex1 = verticesOnEdge(1,:)
        vertex2 = verticesOnEdge(2,:)
        cell1   = cellsOnEdge(1,:)
        cell2   = cellsOnEdge(2,:)
        
        !lambda = dv1Edge/dvEdge
        !pv_edge = dcEdge/(dcEdge+dvEdge)*(lambda*pv_vertex(vertex2) + (1.d0-lambda)*pv_vertex(vertex1)) + dvEdge/(dcEdge+dvEdge)*0.5d0*(pv_cell(cell1) + pv_cell(cell2))
        
        pv_edge = 0.5d0*(pv_vertex(vertex1)+pv_vertex(vertex2))
        
        if(config_apvm_upwinding/=0.d0)then
            !
            ! Compute gradient of PV in the tangent direction
            !   ( this computes gradPVt at all edges bounding real cells and distance-1 ghost cells )
            !
            do iEdge = 1, nEdges
                !gradPVt(iEdge) = (pv_vertex(verticesOnEdge(2,iEdge))*lambda(iEdge) - pv_vertex(verticesOnEdge(1,iEdge))*(1.d0-lambda(iEdge))) / (dvEdge(iEdge)*2.d0*lambda(iEdge)*(1.d0-lambda(iEdge)))
                gradPVt(iEdge) = (pv_vertex(verticesOnEdge(2,iEdge)) - pv_vertex(verticesOnEdge(1,iEdge))) / dvEdge(iEdge)
            enddo
            
            !
            ! Compute gradient of PV in normal direction
            !   ( this computes gradPVn for all edges bounding real cells )
            !
            gradPVn = 0.d0
            do iEdge = 1, nEdges
                gradPVn(iEdge) = (pv_cell(cellsOnEdge(2,iEdge)) - pv_cell(cellsOnEdge(1,iEdge))) / dcEdge(iEdge)
            enddo
            
            !
            ! Modify PV edge with upstream bias. 
            !
            do iEdge = 1, nEdges
                pv_edge(iEdge) = pv_edge(iEdge) - config_apvm_upwinding * v(iEdge) * config_dt * gradPVt(iEdge)
            enddo
            
            !
            ! Modify PV edge with upstream bias.
            !
            do iEdge = 1, nEdges
                pv_edge(iEdge) = pv_edge(iEdge) - config_apvm_upwinding * u(iEdge) * config_dt * gradPVn(iEdge)
            enddo
        endif
            
    end subroutine compute_pv_edge_APVM
    
    !
    ! Compute pv at the edges
    !   ( this computes pv_edge at all edges bounding real cells )
    !
    subroutine compute_pv_edge_APVM_Conservation(u,v,wh_edge,tend_wh_vertex,pv_vertex,pv_cell,pv_edge)
        implicit none
        real(kind=RKIND),intent(in ) :: u              (nEdges   )
        real(kind=RKIND),intent(in ) :: v              (nEdges   )
        real(kind=RKIND),intent(in ) :: wh_edge        (nEdges   )
        real(kind=RKIND),intent(in ) :: tend_wh_vertex (nVertices)
        real(kind=RKIND),intent(in ) :: pv_vertex      (nVertices)
        real(kind=RKIND),intent(in ) :: pv_cell        (nCells   )
        real(kind=RKIND),intent(out) :: pv_edge        (nEdges   )
        
        real(kind=RKIND)             :: pv_edge0       (nEdges   )
        real(kind=RKIND)             :: pv_edge_epsilon(nEdges   )
        real(kind=RKIND)             :: gradPVt        (nEdges   )
        real(kind=RKIND)             :: gradPVn        (nEdges   )
        
        real(kind=RKIND)             :: Q                    (nEdges   )
        real(kind=RKIND)             :: Q0                   (nEdges   )
        real(kind=RKIND)             :: Qepsilon             (nEdges   )
        real(kind=RKIND)             :: Qcombine             (nEdges   )
        real(kind=RKIND)             :: tend_u               (nEdges   )
        real(kind=RKIND)             :: tend_vorticity       (nVertices)
        real(kind=RKIND)             :: tend_vorticity0      (nVertices)
        real(kind=RKIND)             :: tend_vorticityEpsilon(nVertices)
        real(kind=RKIND)             :: tend_vorticityE      (nVertices)
        
        real(kind=RKIND)             :: lambda    (nEdges   )
        integer(kind=IKIND)          :: vertex1   (nEdges   ),&
                                        vertex2   (nEdges   ),&
                                        cell1     (nEdges   ),&
                                        cell2     (nEdges   )
        
        real(kind=RKIND)             :: epsilon
        integer i,iEdge,iVertex
        
        vertex1 = verticesOnEdge(1,:)
        vertex2 = verticesOnEdge(2,:)
        cell1   = cellsOnEdge(1,:)
        cell2   = cellsOnEdge(2,:)
        
        !lambda = dv1Edge/dvEdge
        !pv_edge = dcEdge/(dcEdge+dvEdge)*(lambda*pv_vertex(vertex2) + (1.d0-lambda)*pv_vertex(vertex1)) + dvEdge/(dcEdge+dvEdge)*0.5d0*(pv_cell(cell1) + pv_cell(cell2))
        
        pv_edge0 = 0.5d0*(pv_vertex(vertex1)+pv_vertex(vertex2))
        
        !
        ! Compute gradient of PV in the tangent direction
        !   ( this computes gradPVt at all edges bounding real cells and distance-1 ghost cells )
        !
        !gradPVt(iEdge) = (pv_vertex(verticesOnEdge(2,iEdge))*lambda(iEdge) - pv_vertex(verticesOnEdge(1,iEdge))*(1.d0-lambda(iEdge))) / (dvEdge(iEdge)*2.d0*lambda(iEdge)*(1.d0-lambda(iEdge)))
        gradPVt = (pv_vertex(verticesOnEdge(2,:)) - pv_vertex(verticesOnEdge(1,:))) / dvEdge
        
        !
        ! Compute gradient of PV in normal direction
        !   ( this computes gradPVn for all edges bounding real cells )
        !
        gradPVn = (pv_cell(cellsOnEdge(2,:)) - pv_cell(cellsOnEdge(1,:))) / dcEdge
        
        pv_edge_epsilon = v * gradPVt + u * gradPVn
        
        call compute_Q(u,wh_edge,pv_edge0       ,Q0      )
        call compute_Q(u,wh_edge,pv_edge_epsilon,Qepsilon)
        
        call compute_tend_vorticity_by_tend_u(Q0      ,tend_vorticity0      )
        call compute_tend_vorticity_by_tend_u(Qepsilon,tend_vorticityEpsilon)
        
        epsilon = -(sum(pv_vertex**2*tend_wh_vertex*areaTriangle)-2.d0*sum(pv_vertex*tend_vorticity0*areaTriangle))&
                  /(2.d0*config_dt*sum(pv_vertex*tend_vorticityEpsilon*areaTriangle))
        !epsilon = 0.5d0
        
        !print*,epsilon
        
        !
        ! Modify PV edge with upstream bias. 
        !
        !pv_edge = pv_edge0 - epsilon * v * config_dt * gradPVt
        
        !
        ! Modify PV edge with upstream bias.
        !
        !pv_edge = pv_edge0 - epsilon * u * config_dt * gradPVn
        
        pv_edge = pv_edge0 - epsilon * config_dt * pv_edge_epsilon
        
        !Qcombine = Q0-epsilon*config_dt*Qepsilon
        !call compute_Q(u,wh_edge,pv_edge       ,Q       )
        !tend_u = Q + flux
        !call compute_tend_vorticity_by_tend_u(tend_u ,tend_vorticity      )
        !print*,maxloc(abs((Q-Qcombine)/Q)),maxval(abs((Q-Qcombine)/Q))
        !print*,sum(2.d0*pv_vertex*tend_vorticity*areaTriangle)/sum(areaTriangle),sum(pv_vertex**2*tend_wh_vertex*areaTriangle)/sum(areaTriangle)
            
    end subroutine compute_pv_edge_APVM_Conservation
    
    ! Compute pv_edge with 2nd order accuacy, pv_cl in Weller,2012
    subroutine compute_pv_edge_order2_smooth(pv_vertex, pv_cell, pv_edge)
        implicit none
        real(kind=RKIND),intent(in ) :: pv_vertex(nVertices)
        real(kind=RKIND),intent(in ) :: pv_cell  (nCells   )
        real(kind=RKIND),intent(out) :: pv_edge  (nEdges   )
        
        real(kind=RKIND)             :: lambda    (nEdges   )
        integer(kind=IKIND)          :: vertex1   (nEdges   ),&
                                        vertex2   (nEdges   ),&
                                        cell1     (nEdges   ),&
                                        cell2     (nEdges   )
        
        vertex1 = verticesOnEdge(1,:)
        vertex2 = verticesOnEdge(2,:)
        cell1   = cellsOnEdge(1,:)
        cell2   = cellsOnEdge(2,:)
        
        !lambda = dv1Edge/dvEdge
        !pv_edge = dcEdge/(dcEdge+dvEdge)*(lambda*pv_vertex(vertex2) + (1.d0-lambda)*pv_vertex(vertex1)) + dvEdge/(dcEdge+dvEdge)*0.5d0*(pv_cell(cell1) + pv_cell(cell2))
        pv_edge = 0.5d0*(pv_cell(cell1) + pv_cell(cell2))
    end subroutine compute_pv_edge_order2_smooth
    
    ! Compute pv_edge with 2nd order accuacy, pv_cl in Weller,2012
    subroutine compute_pv_edge_order2(pv_vertex, pv_edge)
        implicit none
        real(kind=RKIND),intent(in ) :: pv_vertex(nVertices)
        real(kind=RKIND),intent(out) :: pv_edge  (nEdges   )
        
        real(kind=RKIND)             :: lambda    (nEdges   )
        integer(kind=IKIND)          :: vertex1   (nEdges   ),&
                                        vertex2   (nEdges   )
        integer iEdge,iVertex
        
        vertex1 = verticesOnEdge(1,:)
        vertex2 = verticesOnEdge(2,:)
        
        lambda = dv1Edge/dvEdge
        
        pv_edge = lambda*pv_vertex(vertex2) + (1.d0-lambda)*pv_vertex(vertex1)
    end subroutine compute_pv_edge_order2
    
    ! Compute pv_edge according to Weller,2012, pv_LUST
    subroutine compute_pv_edge_LUST(u, v, pv_vertex, pv_edge, pvcl_out)
        implicit none
        real(kind=RKIND),intent(in )          :: u        (nEdges   )
        real(kind=RKIND),intent(in )          :: v        (nEdges   )
        real(kind=RKIND),intent(in )          :: pv_vertex(nVertices)
        real(kind=RKIND),intent(out)          :: pv_edge  (nEdges   )
        real(kind=RKIND),intent(out),optional :: pvcl_out (nEdges   )
        
        real(kind=RKIND)             :: pvcl     (nEdges   )  ! Working array
        real(kind=RKIND)             :: pvlu     (nEdges   )  ! Working array
        real(kind=RKIND)             :: nabla_q  (nVertices)  ! Working array
        
        integer(kind=IKIND)          :: vertexUpwind
        real(kind=RKIND)             :: dvUpwind
        integer iEdge,iVertex,Edge1,Edge2,Edge3
        
        call compute_pv_edge_order2(pv_vertex, pvcl)
        if(present(pvcl_out))pvcl_out=pvcl
        
        do iEdge = 1,nEdges
            if(v(iEdge)>0)then
                vertexUpwind = verticesOnEdge(1,iEdge)
                dvUpwind     = dv1Edge(iEdge)
            elseif(v(iEdge)<0)then
                vertexUpwind = verticesOnEdge(2,iEdge)
                dvUpwind     = dv2Edge(iEdge)
            endif
            
            Edge1 = edgesOnVertex(1,vertexUpwind)
            Edge2 = edgesOnVertex(2,vertexUpwind)
            Edge3 = edgesOnVertex(3,vertexUpwind)
            if    (iEdge==edgesOnVertex(1,vertexUpwind))then
                pvlu(iEdge) = pv_vertex(vertexUpwind) + dvUpwind*(dcEdge(Edge1)*pvcl(Edge1)                                 &
                                                                 +dcEdge(Edge2)*pvcl(Edge2)*cos_angleVertex12(vertexUpwind) &
                                                                 +dcEdge(Edge3)*pvcl(Edge3)*cos_angleVertex13(vertexUpwind))&
                                                                 /areaTriangle(vertexUpwind)
            elseif(iEdge==edgesOnVertex(2,vertexUpwind))then
                pvlu(iEdge) = pv_vertex(vertexUpwind) + dvUpwind*(dcEdge(Edge1)*pvcl(Edge1)*cos_angleVertex12(vertexUpwind) &
                                                                 +dcEdge(Edge2)*pvcl(Edge2)                                 &
                                                                 +dcEdge(Edge3)*pvcl(Edge3)*cos_angleVertex23(vertexUpwind))&
                                                                 /areaTriangle(vertexUpwind)
            elseif(iEdge==edgesOnVertex(3,vertexUpwind))then
                pvlu(iEdge) = pv_vertex(vertexUpwind) + dvUpwind*(dcEdge(Edge1)*pvcl(Edge1)*cos_angleVertex13(vertexUpwind) &
                                                                 +dcEdge(Edge2)*pvcl(Edge2)*cos_angleVertex23(vertexUpwind) &
                                                                 +dcEdge(Edge3)*pvcl(Edge3)                                )&
                                                                 /areaTriangle(vertexUpwind)
            endif
            
            if(v(iEdge)==0)pvlu(iEdge) = pvcl(iEdge)
        enddo
        
        pv_edge    = config_CLUST_upwinding*pvlu+(1.d0-config_CLUST_upwinding)*pvcl
    end subroutine compute_pv_edge_LUST
    
    ! Compute pv_edge according to Weller,2012
    subroutine compute_pv_edge_CLUST(u, v, pv_vertex, pv_edge)
        implicit none
        real(kind=RKIND),intent(in ) :: u        (nEdges   )
        real(kind=RKIND),intent(in ) :: v        (nEdges   )
        real(kind=RKIND),intent(in ) :: pv_vertex(nVertices)
        real(kind=RKIND),intent(out) :: pv_edge  (nEdges   )
        
        real(kind=RKIND)             :: pvcl     (nEdges   )  ! Working array
        real(kind=RKIND)             :: pvlu     (nEdges   )  ! Working array
        real(kind=RKIND)             :: pvlu1    (nEdges   )  ! Working array
        real(kind=RKIND)             :: pvlu2    (nEdges   )  ! Working array
        real(kind=RKIND)             :: pvLUST   (nEdges   )  ! Working array
        real(kind=RKIND)             :: nabla_q  (nVertices)  ! Working array
        
        real(kind=RKIND)             :: lambda    (nEdges   )
        real(kind=RKIND)             :: CLUST_para(nEdges   )
        integer(kind=IKIND)          :: vertex1   (nEdges   ),&
                                        vertex2   (nEdges   )
        
        call compute_pv_edge_LUST(u, v, pv_vertex, pvLUST, pvcl)
        
        CLUST_para = abs(u)/dsqrt(u**2+v**2)
        
        pv_edge    = CLUST_para*pvLUST + (1.d0-CLUST_para)*pvcl
        
    end subroutine compute_pv_edge_CLUST
    
    subroutine potentialEnstrophy_Conservation(dkedx,dwhdx,Q,pv_vertex,tend_wh,epsilon)
        implicit none
        real(kind=RKIND),intent(in ) :: dkedx    (nEdges   )
        real(kind=RKIND),intent(in ) :: dwhdx    (nEdges   )
        real(kind=RKIND),intent(in ) :: Q        (nEdges   )
        real(kind=RKIND),intent(in ) :: pv_vertex(nVertices)
        real(kind=RKIND),intent(in ) :: tend_wh  (nCells   )
        real(kind=RKIND),intent(out) :: epsilon
        
        real(kind=RKIND)             :: tend_vorticity_E(nVertices)
        real(kind=RKIND)             :: tend_vorticity_Q(nVertices)
        real(kind=RKIND)             :: tend_wh_vertex  (nVertices)
        
        real(kind=RKIND)             :: flux               (nEdges   )
        real(kind=RKIND)             :: epsilon1,epsilon2,epsilon3
        
        integer(kind=IKIND)          :: vertex1,vertex2,cell1,cell2
        integer(kind=IKIND)          :: i,iVertex,iEdge
        
        ! Compute vorticity tendency for each vertex
        flux = -dkedx-dwhdx
        
        call compute_tend_vorticity_by_tend_u(flux,tend_vorticity_E)
        call compute_tend_vorticity_by_tend_u(Q   ,tend_vorticity_Q)
        
        ! Compute h tend on vertex
        call interp_C2V(tend_wh,tend_wh_vertex)
        
        ! Compute the inner product
        epsilon1 = sum(     tend_wh_vertex  * pv_vertex**2 *areaTriangle)
        epsilon2 = sum(2.d0*tend_vorticity_E* pv_vertex    *areaTriangle)
        epsilon3 = sum(2.d0*tend_vorticity_Q* pv_vertex    *areaTriangle)
        
        epsilon = (epsilon1-epsilon2)/epsilon3
        !if(epsilon<=0.7.or.epsilon>01.2)epsilon=1
        
        !print*,'epsilon   = ',epsilon,' epsilon1  = ',epsilon1,' epsilon2  = ',epsilon2,' epsilon3  = ',epsilon3
        !print*,'epsilon   = ',epsilon
    end subroutine potentialEnstrophy_Conservation
    
    function inner_product(u1,h1,u2,h2)
        implicit none
        real(kind=RKIND), intent(in) :: u1(nEdges)
        real(kind=RKIND), intent(in) :: u2(nEdges)
        real(kind=RKIND), intent(in) :: h1(nCells)
        real(kind=RKIND), intent(in) :: h2(nCells)
        real(kind=RKIND)             :: inner_product
        
        inner_product = (sum(u1*u2*areaEdge) + sum(h1*h2*areaCell))/sum_areaCell
        
    end function inner_product

    function totalMass(wh)
        implicit none
        real(kind=RKIND), intent(in) :: wh(nCells)
        real(kind=RKIND)             :: totalMass
        
        totalMass = sum(wh*areaCell)/sum(areaCell)
    end function totalMass
    
    function totalAbsoluteVorticity(wu,wh)
        implicit none
        real(kind=RKIND), intent(in) :: wu(nEdges)
        real(kind=RKIND), intent(in) :: wh(nCells)
        real(kind=RKIND)             :: totalAbsoluteVorticity
        
        real(kind=RKIND)             :: u(nEdges)
        real(kind=RKIND)             :: vorticity(nVertices)
        
        call IIAP(wu,wh,u)
        call compute_vorticity(u,vorticity)
        
        totalAbsoluteVorticity = sum((vorticity+fVertex)*areaTriangle)/sum(areaTriangle)
    end function totalAbsoluteVorticity
    
    function totalEnergy(wu,wh,wh_s)
        implicit none
        real(kind=RKIND), intent(in) :: wu  (nEdges)
        real(kind=RKIND), intent(in) :: wh  (nCells)
        real(kind=RKIND), intent(in) :: wh_s(nCells)
        real(kind=RKIND)             :: totalEnergy
        
        totalEnergy = (sum(wu**2*areaEdge) + sum((wh+wh_s)**2*areaCell))/sum(areaCell)
    end function totalEnergy
    
    function totalPotentialEnstrophy(wu,wh)
        implicit none
        real(kind=RKIND), intent(in) :: wu(nEdges)
        real(kind=RKIND), intent(in) :: wh(nCells)
        real(kind=RKIND)             :: totalPotentialEnstrophy
        
        real(kind=RKIND)             :: u(nEdges)
        real(kind=RKIND)             :: vorticity(nVertices)
        real(kind=RKIND)             :: pv_vertex(nVertices)
        real(kind=RKIND)             :: wh_vertex(nVertices)
        
        call IIAP(wu,wh,u)
        call compute_vorticity(u,vorticity)
        call interp_C2V(wh,wh_vertex)
        call compute_pv_vertex(u,wh_vertex,vorticity,pv_vertex)
        
        totalPotentialEnstrophy = sum(wh_vertex*pv_vertex**2*areaTriangle)/sum(areaTriangle)
    end function totalPotentialEnstrophy
    
    function totalAngularMomentum(wu,wh)
        implicit none
        real(kind=RKIND), intent(in) :: wu(nEdges)
        real(kind=RKIND), intent(in) :: wh(nCells)
        real(kind=RKIND)             :: totalAngularMomentum
        
        real(kind=RKIND)             :: u      (nEdges)
        real(kind=RKIND)             :: wh_edge(nEdges)
        
        call interp_C2E(wh,wh_edge)
        
        totalAngularMomentum = sum(wh_edge*(u*angleEdge+a*omega*dcos(latEdge)**2)*areaEdge)/sum(areaEdge)
    end function totalAngularMomentum
    
end module diagnostics_tools