module mesh_info
    use netcdf
    use module_para,only: RKIND,a,omega,pi,config_mesh_file
    implicit none
    include 'netcdf.inc'
    
    integer nCells       ,&
            nEdges       ,&
            nVertices    ,&
            vertexDegree ,&
            maxEdges     ,&
            maxEdges2
    
    integer latCell_id          ,&
            lonCell_id          ,&
            latEdge_id          ,&
            lonEdge_id          ,&
            latVertex_id        ,&
            lonVertex_id        ,&
            dvEdge_id           ,&
            dv1Edge_id          ,&
            dv2Edge_id           ,&
            dcEdge_id           ,&
            angleEdge_id        ,&
            areaCell_id         ,&
            areaTriangle_id     ,&
            nEdgesOnCell_id     ,&
            nEdgesOnEdge_id     ,&
            cellsOnEdge_id      ,&
            edgesOnCell_id      ,&
            edgesOnEdge_id      ,&
            weightsOnEdge_id    ,&
            cellsOnCell_id      ,&
            verticesOnCell_id   ,&
            verticesOnEdge_id   ,&
            edgesOnVertex_id    ,&
            cellsOnVertex_id    ,&
            kiteAreasOnVertex_id
            
    real(kind=RKIND),allocatable :: latCell          (:)
    real(kind=RKIND),allocatable :: lonCell          (:)
    real(kind=RKIND),allocatable :: latEdge          (:)
    real(kind=RKIND),allocatable :: lonEdge          (:)
    real(kind=RKIND),allocatable :: latVertex        (:)
    real(kind=RKIND),allocatable :: lonVertex        (:)
    real(kind=RKIND),allocatable :: dvEdge           (:)
    real(kind=RKIND),allocatable :: dv1Edge           (:)
    real(kind=RKIND),allocatable :: dv2Edge           (:)
    real(kind=RKIND),allocatable :: dcEdge           (:)
    real(kind=RKIND),allocatable :: angleEdge        (:)
    real(kind=RKIND),allocatable :: areaCell         (:)
    real(kind=RKIND),allocatable :: areaTriangle     (:)
    real(kind=RKIND),allocatable :: areaEdge         (:)
    real(kind=RKIND),allocatable :: kiteAreasOnVertex(:,:)
    real(kind=RKIND),allocatable :: fCell            (:)
    real(kind=RKIND),allocatable :: fEdge            (:)
    real(kind=RKIND),allocatable :: fVertex          (:)
    real(kind=RKIND),allocatable :: weightsOnEdge    (:,:)
    real(kind=RKIND),allocatable :: cos_angleVertex12(:)
    real(kind=RKIND),allocatable :: cos_angleVertex13(:)
    real(kind=RKIND),allocatable :: cos_angleVertex23(:)
    
    integer,allocatable :: nEdgesOnCell    (:)
    integer,allocatable :: nEdgesOnEdge    (:)
    
    integer,allocatable :: cellsOnEdge     (:,:)
    integer,allocatable :: edgesOnCell     (:,:)
    integer,allocatable :: edgesOnEdge     (:,:)
    integer,allocatable :: cellsOnCell     (:,:)
    integer,allocatable :: verticesOnCell  (:,:)
    integer,allocatable :: verticesOnEdge  (:,:)
    integer,allocatable :: edgesOnVertex   (:,:)
    integer,allocatable :: cellsOnVertex   (:,:)
    integer,allocatable :: t               (:,:)
    integer,allocatable :: n               (:,:)
    real(kind=RKIND)    :: sum_areaEdge
    real(kind=RKIND)    :: sum_areaCell
    
    integer                      :: iCell,iEdge,iVertex,vertex1,vertex2,Edge1,Edge2,Edge3
    real(kind=RKIND),allocatable :: angleEdgeOnVertex(:,:)
    
    contains
    
    subroutine read_mesh
    implicit none
    integer ncid,status
    integer nCells_dimid      ,&
            nEdges_dimid      ,&
            nVertices_dimid   ,&
            vertexDegree_dimid,&
            maxEdges_dimid    ,&
            maxEdges_dimid2
    
    Print*,'Reading mesh info'
    
    status = nf90_open(trim(adjustl(config_mesh_file)),NF90_NOWRITE,ncid)
    !if(status/=nf90_noerr) call handle_err(status)
    
    status = nf90_inq_dimid(ncid,'nCells'      ,nCells_dimid      )
    status = nf90_inq_dimid(ncid,'nEdges'      ,nEdges_dimid      )
    status = nf90_inq_dimid(ncid,'nVertices'   ,nVertices_dimid   )
    status = nf90_inq_dimid(ncid,'vertexDegree',vertexDegree_dimid)
    status = nf90_inq_dimid(ncid,'maxEdges'    ,maxEdges_dimid    )
    status = nf90_inq_dimid(ncid,'maxEdges2'   ,maxEdges_dimid2   )
    
    status = nf90_inquire_dimension(ncid,nCells_dimid      ,len=nCells      )
    status = nf90_inquire_dimension(ncid,nEdges_dimid      ,len=nEdges      )
    status = nf90_inquire_dimension(ncid,nVertices_dimid   ,len=nVertices   )
    status = nf90_inquire_dimension(ncid,vertexDegree_dimid,len=vertexDegree)
    status = nf90_inquire_dimension(ncid,maxEdges_dimid    ,len=maxEdges    )
    status = nf90_inquire_dimension(ncid,maxEdges_dimid2   ,len=maxEdges2   )
    
    ! Plot mesh info
    print*,'Using mesh file '//trim(adjustl(config_mesh_file)),' with'
    print*,'nCells       = ',nCells
    print*,'nEdges       = ',nEdges
    print*,'nVertices    = ',nVertices
    print*,'vertexDegree = ',vertexDegree
    print*,'maxEdges     = ',maxEdges
    print*,'maxEdges2    = ',maxEdges2
    
    allocate(                                 &
    latCell          (nCells                ),&
    lonCell          (nCells                ),&
    latEdge          (nEdges                ),&
    lonEdge          (nEdges                ),&
    latVertex        (nVertices             ),&
    lonVertex        (nVertices             ),&
    dvEdge           (nEdges                ),&
    dv1Edge          (nEdges                ),&
    dv2Edge          (nEdges                ),&
    dcEdge           (nEdges                ),&
    angleEdge        (nEdges                ),&
    areaCell         (nCells                ),&
    areaTriangle     (nVertices             ),&
    areaEdge         (nEdges                ),&
    nEdgesOnCell     (nCells                ),&
    nEdgesOnEdge     (nEdges                ),&
    cellsOnEdge      (2          ,nEdges    ),&
    edgesOnCell      (maxEdges   ,nCells    ),&
    edgesOnEdge      (maxEdges2  ,nEdges    ),&
    weightsOnEdge    (maxEdges2  ,nEdges    ),&
    cellsOnCell      (maxEdges   ,nCells    ),&
    verticesOnCell   (maxEdges   ,nCells    ),&
    verticesOnEdge   (2          ,nEdges    ),&
    edgesOnVertex    (vertexDegree,nVertices),&
    cellsOnVertex    (vertexDegree,nVertices),&
    kiteAreasOnVertex(vertexDegree,nVertices),&
    fCell            (nCells                ),&
    fEdge            (nEdges                ),&
    fVertex          (nVertices             ),&
    cos_angleVertex12(nVertices             ),&
    cos_angleVertex13(nVertices             ),&
    cos_angleVertex23(nVertices             ),&
    t                (vertexDegree,nVertices),&
    n                (maxEdges    ,nCells   ) &
    )
    
    status = nf90_inq_varid(ncid,'latCell'          ,latCell_id          )
    status = nf90_inq_varid(ncid,'lonCell'          ,lonCell_id          )
    status = nf90_inq_varid(ncid,'latEdge'          ,latEdge_id          )
    status = nf90_inq_varid(ncid,'lonEdge'          ,lonEdge_id          )
    status = nf90_inq_varid(ncid,'latVertex'        ,latVertex_id        )
    status = nf90_inq_varid(ncid,'lonVertex'        ,lonVertex_id        )
    status = nf90_inq_varid(ncid,'dvEdge'           ,dvEdge_id           )
    status = nf90_inq_varid(ncid,'dv1Edge'          ,dv1Edge_id          )
    status = nf90_inq_varid(ncid,'dv2Edge'          ,dv2Edge_id          )
    status = nf90_inq_varid(ncid,'dcEdge'           ,dcEdge_id           )
    status = nf90_inq_varid(ncid,'angleEdge'        ,angleEdge_id        )
    status = nf90_inq_varid(ncid,'areaCell'         ,areaCell_id         )
    status = nf90_inq_varid(ncid,'areaTriangle'     ,areaTriangle_id     )
    status = nf90_inq_varid(ncid,'nEdgesOnCell'     ,nEdgesOnCell_id     )
    status = nf90_inq_varid(ncid,'nEdgesOnEdge'     ,nEdgesOnEdge_id     )
    status = nf90_inq_varid(ncid,'cellsOnEdge'      ,cellsOnEdge_id      )
    status = nf90_inq_varid(ncid,'edgesOnCell'      ,edgesOnCell_id      )
    status = nf90_inq_varid(ncid,'edgesOnEdge'      ,edgesOnEdge_id      )
    status = nf90_inq_varid(ncid,'weightsOnEdge'    ,weightsOnEdge_id    )
    status = nf90_inq_varid(ncid,'cellsOnCell'      ,cellsOnCell_id      )
    status = nf90_inq_varid(ncid,'verticesOnCell'   ,verticesOnCell_id   )
    status = nf90_inq_varid(ncid,'verticesOnEdge'   ,verticesOnEdge_id   )
    status = nf90_inq_varid(ncid,'edgesOnVertex'    ,edgesOnVertex_id    )
    status = nf90_inq_varid(ncid,'cellsOnVertex'    ,cellsOnVertex_id    )
    status = nf90_inq_varid(ncid,'kiteAreasOnVertex',kiteAreasOnVertex_id)
    
    status = nf90_get_var(ncid,latCell_id          ,latCell          )
    status = nf90_get_var(ncid,lonCell_id          ,lonCell          )
    status = nf90_get_var(ncid,latEdge_id          ,latEdge          )
    status = nf90_get_var(ncid,lonEdge_id          ,lonEdge          )
    status = nf90_get_var(ncid,latVertex_id        ,latVertex        )
    status = nf90_get_var(ncid,lonVertex_id        ,lonVertex        )
    status = nf90_get_var(ncid,dvEdge_id           ,dvEdge           )
    status = nf90_get_var(ncid,dv1Edge_id          ,dv1Edge          )
    status = nf90_get_var(ncid,dv2Edge_id          ,dv2Edge          )
    status = nf90_get_var(ncid,dcEdge_id           ,dcEdge           )
    status = nf90_get_var(ncid,angleEdge_id        ,angleEdge        )
    status = nf90_get_var(ncid,areaCell_id         ,areaCell         )
    status = nf90_get_var(ncid,areaTriangle_id     ,areaTriangle     )
    status = nf90_get_var(ncid,nEdgesOnCell_id     ,nEdgesOnCell     )
    status = nf90_get_var(ncid,nEdgesOnEdge_id     ,nEdgesOnEdge     )
    status = nf90_get_var(ncid,cellsOnEdge_id      ,cellsOnEdge      )
    status = nf90_get_var(ncid,edgesOnCell_id      ,edgesOnCell      )
    status = nf90_get_var(ncid,edgesOnEdge_id      ,edgesOnEdge      )
    status = nf90_get_var(ncid,weightsOnEdge_id    ,weightsOnEdge    )
    status = nf90_get_var(ncid,cellsOnCell_id      ,cellsOnCell      )
    status = nf90_get_var(ncid,verticesOnCell_id   ,verticesOnCell   )
    status = nf90_get_var(ncid,verticesOnEdge_id   ,verticesOnEdge   )
    status = nf90_get_var(ncid,edgesOnVertex_id    ,edgesOnVertex    )
    status = nf90_get_var(ncid,cellsOnVertex_id    ,cellsOnVertex    )
    status = nf90_get_var(ncid,kiteAreasOnVertex_id,kiteAreasOnVertex)
    
    fCell    = 2.d0*omega*dsin(latCell)
    fEdge    = 2.d0*omega*dsin(latEdge)
    fVertex  = 2.d0*omega*dsin(latVertex)

    !
    ! Scale all distances and areas from a unit sphere to one with radius a
    !

    dvEdge            = dvEdge            * a
    dv1Edge           = dv1Edge           * a
    dv2Edge           = dv2Edge           * a
    dcEdge            = dcEdge            * a
    areaCell          = areaCell          * a**2.d0
    areaTriangle      = areaTriangle      * a**2.d0
    areaEdge          = dvEdge*dcEdge
    kiteAreasOnVertex = kiteAreasOnVertex * a**2.d0
    sum_areaEdge      = sum(areaEdge)
    sum_areaCell      = sum(areaCell)
    
    allocate(angleEdgeOnVertex(3,nVertices))
    
    do iEdge = 1,nEdges
        vertex1 = verticesOnEdge(1,iEdge)
        vertex2 = verticesOnEdge(2,iEdge)
        
        if(iEdge==edgesOnVertex(1,vertex1))angleEdgeOnVertex(1,Vertex1) = angleEdge(iEdge)
        if(iEdge==edgesOnVertex(2,vertex1))angleEdgeOnVertex(2,Vertex1) = angleEdge(iEdge)
        if(iEdge==edgesOnVertex(3,vertex1))angleEdgeOnVertex(3,Vertex1) = angleEdge(iEdge)
        
        if(iEdge==edgesOnVertex(1,vertex2))angleEdgeOnVertex(1,Vertex2) = angleEdge(iEdge) + pi
        if(iEdge==edgesOnVertex(2,vertex2))angleEdgeOnVertex(2,Vertex2) = angleEdge(iEdge) + pi
        if(iEdge==edgesOnVertex(3,vertex2))angleEdgeOnVertex(3,Vertex2) = angleEdge(iEdge) + pi
    enddo
    
    cos_angleVertex12 = dcos(angleEdgeOnVertex(1,:) - angleEdgeOnVertex(2,:))
    cos_angleVertex13 = dcos(angleEdgeOnVertex(1,:) - angleEdgeOnVertex(3,:))
    cos_angleVertex23 = dcos(angleEdgeOnVertex(2,:) - angleEdgeOnVertex(3,:))
    
    t=0
    do iVertex = 1,nVertices
        do iEdge = 1,vertexDegree
            if(iVertex==verticesOnEdge(1,edgesOnVertex(iEdge,iVertex))) t(iEdge,iVertex)=1
            if(iVertex==verticesOnEdge(2,edgesOnVertex(iEdge,iVertex))) t(iEdge,iVertex)=-1
        enddo
    enddo
    
    n=0
    do iCell = 1,nCells
        do iEdge = 1,nEdgesOnCell(iCell)
            if(iCell==cellsOnEdge(1,edgesOnCell(iEdge,iCell))) n(iEdge,iCell)=1
            if(iCell==cellsOnEdge(2,edgesOnCell(iEdge,iCell))) n(iEdge,iCell)=-1
        enddo
    enddo
    
    end subroutine read_mesh

    real(kind=RKIND) function sphere_distance(lat1, lon1, lat2, lon2, radius)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
    !   sphere with given radius.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        implicit none
        
        real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius
        
        real (kind=RKIND) :: arg1
        
        arg1            = sqrt( dsin(0.5d0*(lat2-lat1))**2.d0 + dcos(lat1)*dcos(lat2)*dsin(0.5d0*(lon2-lon1))**2.d0 )
        sphere_distance = 2.d0*radius*dasin(arg1)
    
    end function sphere_distance
    
    subroutine handle_err(status)
        implicit none
        integer,intent(in)::status
        
        if(status/=nf90_noerr)then
            print*, trim(nf90_strerror(status))
            stop "Stopped by netCDF"
        endif
                    
    endsubroutine handle_err
    
end module mesh_info
    