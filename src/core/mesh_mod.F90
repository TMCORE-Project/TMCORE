module mesh_mod

  use const_mod
  use params_mod

  implicit none

  ! Dimension sizes
  integer nCells
  integer nEdges
  integer nVertices
  integer vertexDegree
  integer maxEdges
  integer maxEdges2
  integer, allocatable :: nEdgesOnCell(:)                 ! Number of edges on a given cell
  integer, allocatable :: nEdgesOnEdge(:)                 ! Number of edges on a given edge to reconstruct tangential velocities
  integer, allocatable :: nCellsOnVertex(:)               ! Number of cells connected with a given vertex
  ! Coordinates
  real(real_kind), allocatable :: latCell(:)
  real(real_kind), allocatable :: lonCell(:)
  real(real_kind), allocatable :: xCell(:)
  real(real_kind), allocatable :: yCell(:)
  real(real_kind), allocatable :: zCell(:)
  real(real_kind), allocatable :: latEdge(:)
  real(real_kind), allocatable :: lonEdge(:)
  real(real_kind), allocatable :: xEdge(:)
  real(real_kind), allocatable :: yEdge(:)
  real(real_kind), allocatable :: zEdge(:)
  real(real_kind), allocatable :: latVertex(:)
  real(real_kind), allocatable :: lonVertex(:)
  real(real_kind), allocatable :: xVertex(:)
  real(real_kind), allocatable :: yVertex(:)
  real(real_kind), allocatable :: zVertex(:)
  ! Geometric measures
  real(real_kind), allocatable :: dvEdge(:)               ! Distance in meters between the vertices that saddle a given edge
  real(real_kind), allocatable :: dv1Edge(:)              ! Distance in meters between vertex 1 and edge point
  real(real_kind), allocatable :: dv2Edge(:)              ! Distance in meters between vertex 2 and edge point
  real(real_kind), allocatable :: dcEdge(:)               ! Distance in meters between the cells that saddle a given edge
  real(real_kind), allocatable :: areaCell(:)             ! Area in square meters for a given cell of the primary mesh
  real(real_kind), allocatable :: areaEdge(:)             ! Area in square meters for a given edge point
  real(real_kind), allocatable :: areaTriangle(:)         ! Area in square meters for a given triangle of the dual mesh
  real(real_kind), allocatable :: kiteAreasOnVertex(:,:)  ! The intersection area of areaTriangle with each cell that radiates from a given vertex
  real(real_kind), allocatable :: angleEdge(:)            ! Angle in radians an edge’s normal vector makes with the local eastward direction
  integer, allocatable :: nSignEdge(:,:)
  integer, allocatable :: tSignEdge(:,:)
  real(real_kind) totalArea
  ! Indices
  integer, allocatable :: indexToCellID(:)                ! Global cell ID for all cell centers
  integer, allocatable :: indexToEdgeID(:)                ! Global edge ID for all edge locations
  integer, allocatable :: indexToVertexID(:)              ! Global vertex ID for all cell vertices
  integer, allocatable :: cellsOnCell(:,:)                ! Cell indices that surround a given cell
  integer, allocatable :: cellsOnEdge(:,:)                ! Cell indices that saddle a given edge
  integer, allocatable :: cellsOnVertex(:,:)              ! Cell indices that radiate from a given vertex
  integer, allocatable :: edgesOnCell(:,:)                ! Edge indices that surround a given cell
  integer, allocatable :: edgesOnEdge(:,:)                ! Edge indices that are used to reconstruct tangential velocities
  integer, allocatable :: edgesOnVertex(:,:)              ! Edge indices that radiate from a given vertex
  integer, allocatable :: verticesOnCell(:,:)             ! Vertex indices that surround a given cell
  integer, allocatable :: verticesOnEdge(:,:)             ! Vertex indices that saddle a given edge
  ! Weights
  real(real_kind), allocatable :: weightsOnEdge(:,:)      ! Weights to reconstruct tangential velocities
  real(real_kind), allocatable :: meshDensity(:)          ! The value of the generating density function at each cell center

  real(real_kind), allocatable :: fCell(:)                ! Coriolis coefficients on a given cell
  real(real_kind), allocatable :: fVertex(:)              ! Coriolis coefficients on a given vertex

contains

  subroutine mesh_init()

    use netcdf

    integer ncid, ierr
    integer dimid, varid

    integer iCell, iEdge, iVertex, i

    ierr = nf90_open(mesh_file_path, nf90_nowrite, ncid)

    ierr = nf90_inq_dimid(ncid, 'nCells', dimid)

    ierr = nf90_inquire_dimension(ncid, dimid, len=nCells)

    ierr = nf90_inq_dimid(ncid, 'nEdges', dimid)

    ierr = nf90_inquire_dimension(ncid, dimid, len=nEdges)

    ierr = nf90_inq_dimid(ncid, 'nVertices', dimid)

    ierr = nf90_inquire_dimension(ncid, dimid, len=nVertices)

    ierr = nf90_inq_dimid(ncid, 'vertexDegree', dimid)

    ierr = nf90_inquire_dimension(ncid, dimid, len=vertexDegree)

    ierr = nf90_inq_dimid(ncid, 'maxEdges', dimid)

    ierr = nf90_inquire_dimension(ncid, dimid, len=maxEdges)

    ierr = nf90_inq_dimid(ncid, 'maxEdges2', dimid)

    ierr = nf90_inquire_dimension(ncid, dimid, len=maxEdges2)

    allocate(nEdgesOnCell(nCells))
    allocate(nEdgesOnEdge(nEdges))
    allocate(latCell(nCells))
    allocate(lonCell(nCells))
    allocate(xCell(nCells))
    allocate(yCell(nCells))
    allocate(zCell(nCells))
    allocate(latEdge(nEdges))
    allocate(lonEdge(nEdges))
    allocate(xEdge(nEdges))
    allocate(yEdge(nEdges))
    allocate(zEdge(nEdges))
    allocate(latVertex(nVertices))
    allocate(lonVertex(nVertices))
    allocate(xVertex(nVertices))
    allocate(yVertex(nVertices))
    allocate(zVertex(nVertices))
    allocate(dvEdge(nEdges))
    allocate(dv1Edge(nEdges))
    allocate(dv2Edge(nEdges))
    allocate(dcEdge(nEdges))
    allocate(areaCell(nCells))
    allocate(areaTriangle(nVertices))
    allocate(kiteAreasOnVertex(vertexDegree,nVertices))
    allocate(angleEdge(nEdges))
    allocate(indexToCellID(nCells))
    allocate(indexToEdgeID(nEdges))
    allocate(indexToVertexID(nVertices))
    allocate(cellsOnCell(maxEdges,nCells))
    allocate(cellsOnEdge(2,nEdges))
    allocate(cellsOnVertex(vertexDegree,nVertices))
    allocate(edgesOnCell(maxEdges,nCells))
    allocate(edgesOnEdge(maxEdges2,nEdges))
    allocate(edgesOnVertex(vertexDegree,nVertices))
    allocate(verticesOnCell(maxEdges,nCells))
    allocate(verticesOnEdge(2,nEdges))
    allocate(weightsOnEdge(maxEdges2,nEdges))
    allocate(meshDensity(nCells))

    ierr = nf90_inq_varid(ncid, 'nEdgesOnCell', varid)

    ierr = nf90_get_var(ncid, varid, nEdgesOnCell)

    ierr = nf90_inq_varid(ncid, 'nEdgesOnEdge', varid)

    ierr = nf90_get_var(ncid, varid, nEdgesOnEdge)

    ierr = nf90_inq_varid(ncid, 'latCell', varid)

    ierr = nf90_get_var(ncid, varid, latCell)

    ierr = nf90_inq_varid(ncid, 'lonCell', varid)

    ierr = nf90_get_var(ncid, varid, lonCell)

    ierr = nf90_inq_varid(ncid, 'xCell', varid)

    ierr = nf90_get_var(ncid, varid, xCell)

    ierr = nf90_inq_varid(ncid, 'yCell', varid)

    ierr = nf90_get_var(ncid, varid, yCell)

    ierr = nf90_inq_varid(ncid, 'zCell', varid)

    ierr = nf90_get_var(ncid, varid, zCell)

    ierr = nf90_inq_varid(ncid, 'latEdge', varid)

    ierr = nf90_get_var(ncid, varid, latEdge)

    ierr = nf90_inq_varid(ncid, 'lonEdge', varid)

    ierr = nf90_get_var(ncid, varid, lonEdge)

    ierr = nf90_inq_varid(ncid, 'xEdge', varid)

    ierr = nf90_get_var(ncid, varid, xEdge)

    ierr = nf90_inq_varid(ncid, 'yEdge', varid)

    ierr = nf90_get_var(ncid, varid, yEdge)

    ierr = nf90_inq_varid(ncid, 'zEdge', varid)

    ierr = nf90_get_var(ncid, varid, zEdge)

    ierr = nf90_inq_varid(ncid, 'latVertex', varid)

    ierr = nf90_get_var(ncid, varid, latVertex)

    ierr = nf90_inq_varid(ncid, 'lonVertex', varid)

    ierr = nf90_get_var(ncid, varid, lonVertex)

    ierr = nf90_inq_varid(ncid, 'xVertex', varid)

    ierr = nf90_get_var(ncid, varid, xVertex)

    ierr = nf90_inq_varid(ncid, 'yVertex', varid)

    ierr = nf90_get_var(ncid, varid, yVertex)

    ierr = nf90_inq_varid(ncid, 'zVertex', varid)

    ierr = nf90_get_var(ncid, varid, zVertex)

    ierr = nf90_inq_varid(ncid, 'dvEdge', varid)

    ierr = nf90_get_var(ncid, varid, dvEdge)

    ierr = nf90_inq_varid(ncid, 'dv1Edge', varid)

    ierr = nf90_get_var(ncid, varid, dv1Edge)

    ierr = nf90_inq_varid(ncid, 'dv2Edge', varid)

    ierr = nf90_get_var(ncid, varid, dv2Edge)

    ierr = nf90_inq_varid(ncid, 'dcEdge', varid)

    ierr = nf90_get_var(ncid, varid, dcEdge)

    ierr = nf90_inq_varid(ncid, 'areaCell', varid)

    ierr = nf90_get_var(ncid, varid, areaCell)

    ierr = nf90_inq_varid(ncid, 'areaTriangle', varid)

    ierr = nf90_get_var(ncid, varid, areaTriangle)

    ierr = nf90_inq_varid(ncid, 'kiteAreasOnVertex', varid)

    ierr = nf90_get_var(ncid, varid, kiteAreasOnVertex)

    ierr = nf90_inq_varid(ncid, 'angleEdge', varid)

    ierr = nf90_get_var(ncid, varid, angleEdge)

    ierr = nf90_inq_varid(ncid, 'indexToCellID', varid)

    ierr = nf90_get_var(ncid, varid, indexToCellID)

    ierr = nf90_inq_varid(ncid, 'indexToEdgeID', varid)

    ierr = nf90_get_var(ncid, varid, indexToEdgeID)

    ierr = nf90_inq_varid(ncid, 'indexToVertexID', varid)

    ierr = nf90_get_var(ncid, varid, indexToVertexID)

    ierr = nf90_inq_varid(ncid, 'cellsOnCell', varid)

    ierr = nf90_get_var(ncid, varid, cellsOnCell)

    ierr = nf90_inq_varid(ncid, 'cellsOnEdge', varid)

    ierr = nf90_get_var(ncid, varid, cellsOnEdge)

    ierr = nf90_inq_varid(ncid, 'cellsOnVertex', varid)

    ierr = nf90_get_var(ncid, varid, cellsOnVertex)

    ierr = nf90_inq_varid(ncid, 'edgesOnCell', varid)

    ierr = nf90_get_var(ncid, varid, edgesOnCell)

    ierr = nf90_inq_varid(ncid, 'edgesOnEdge', varid)

    ierr = nf90_get_var(ncid, varid, edgesOnEdge)

    ierr = nf90_inq_varid(ncid, 'edgesOnVertex', varid)

    ierr = nf90_get_var(ncid, varid, edgesOnVertex)

    ierr = nf90_inq_varid(ncid, 'verticesOnCell', varid)

    ierr = nf90_get_var(ncid, varid, verticesOnCell)

    ierr = nf90_inq_varid(ncid, 'verticesOnEdge', varid)

    ierr = nf90_get_var(ncid, varid, verticesOnEdge)

    ierr = nf90_inq_varid(ncid, 'weightsOnEdge', varid)

    ierr = nf90_get_var(ncid, varid, weightsOnEdge)

    ierr = nf90_inq_varid(ncid, 'meshDensity', varid)

    ierr = nf90_get_var(ncid, varid, meshDensity)

    ierr = nf90_close(ncid)

    ! Derived quantities
    allocate(nCellsOnVertex(nVertices))
    allocate(areaEdge(nEdges))
    allocate(fCell(nCells))
    allocate(fVertex(nVertices))
    allocate(nSignEdge(maxEdges,nCells))
    allocate(tSignEdge(vertexDegree,nVertices))

    nCellsOnVertex = vertexDegree
    areaEdge       = dvEdge(:) * dcEdge(:)
    fCell          = 2.0d0 * omega * sin(latCell(:))
    fVertex        = 2.0d0 * omega * sin(latVertex(:))

    nSignEdge = 0
    do iCell = 1, nCells
      do i = 1, nEdgesOnCell(iCell)
        if (iCell == cellsOnEdge(1,edgesOnCell(i,iCell))) nSignEdge(i,iCell) =  1
        if (iCell == cellsOnEdge(2,edgesOnCell(i,iCell))) nSignEdge(i,iCell) = -1
      end do
    end do

    tSignEdge = 0
    do iVertex = 1, nVertices
      do i = 1, vertexDegree
        if (iVertex == verticesOnEdge(1,edgesOnVertex(i,iVertex))) tSignEdge(i,iVertex) =  1
        if (iVertex == verticesOnEdge(2,edgesOnVertex(i,iVertex))) tSignEdge(i,iVertex) = -1
      end do
    end do

    ! Scale mesh parameters.
    xCell             = xCell             * radius
    yCell             = yCell             * radius
    zCell             = zCell             * radius
    xEdge             = xEdge             * radius
    yEdge             = yEdge             * radius
    zEdge             = zEdge             * radius
    xVertex           = xVertex           * radius
    yVertex           = yVertex           * radius
    zVertex           = zVertex           * radius
    dvEdge            = dvEdge            * radius
    dv1Edge           = dv1Edge           * radius
    dv2Edge           = dv2Edge           * radius
    dcEdge            = dcEdge            * radius
    areaCell          = areaCell          * radius**2
    areaTriangle      = areaTriangle      * radius**2
    areaEdge          = areaEdge          * radius**2
    kiteAreasOnVertex = kiteAreasOnVertex * radius**2

    totalArea         = sum(areaCell)

  end subroutine mesh_init

  subroutine mesh_final()

    if (allocated(nEdgesOnCell))      deallocate(nEdgesOnCell)
    if (allocated(nEdgesOnEdge))      deallocate(nEdgesOnEdge)
    if (allocated(latCell))           deallocate(latCell)
    if (allocated(lonCell))           deallocate(lonCell)
    if (allocated(xCell))             deallocate(xCell)
    if (allocated(yCell))             deallocate(yCell)
    if (allocated(zCell))             deallocate(zCell)
    if (allocated(latEdge))           deallocate(latEdge)
    if (allocated(lonEdge))           deallocate(lonEdge)
    if (allocated(xEdge))             deallocate(xEdge)
    if (allocated(yEdge))             deallocate(yEdge)
    if (allocated(zEdge))             deallocate(zEdge)
    if (allocated(latVertex))         deallocate(latVertex)
    if (allocated(lonVertex))         deallocate(lonVertex)
    if (allocated(xVertex))           deallocate(xVertex)
    if (allocated(yVertex))           deallocate(yVertex)
    if (allocated(zVertex))           deallocate(zVertex)
    if (allocated(dvEdge))            deallocate(dvEdge)
    if (allocated(dcEdge))            deallocate(dcEdge)
    if (allocated(areaCell))          deallocate(areaCell)
    if (allocated(areaTriangle))      deallocate(areaTriangle)
    if (allocated(kiteAreasOnVertex)) deallocate(kiteAreasOnVertex)
    if (allocated(angleEdge))         deallocate(angleEdge)
    if (allocated(indexToCellID))     deallocate(indexToCellID)
    if (allocated(indexToEdgeID))     deallocate(indexToEdgeID)
    if (allocated(indexToVertexID))   deallocate(indexToVertexID)
    if (allocated(cellsOnCell))       deallocate(cellsOnCell)
    if (allocated(cellsOnEdge))       deallocate(cellsOnEdge)
    if (allocated(cellsOnVertex))     deallocate(cellsOnVertex)
    if (allocated(edgesOnCell))       deallocate(edgesOnCell)
    if (allocated(edgesOnEdge))       deallocate(edgesOnEdge)
    if (allocated(edgesOnVertex))     deallocate(edgesOnVertex)
    if (allocated(verticesOnCell))    deallocate(verticesOnCell)
    if (allocated(verticesOnEdge))    deallocate(verticesOnEdge)
    if (allocated(weightsOnEdge))     deallocate(weightsOnEdge)
    if (allocated(meshDensity))       deallocate(meshDensity)
    if (allocated(nCellsOnVertex))    deallocate(nCellsOnVertex)
    if (allocated(areaEdge))          deallocate(areaEdge)
    if (allocated(fCell))             deallocate(fCell)
    if (allocated(fVertex))           deallocate(fVertex)
    if (allocated(nSignEdge))         deallocate(nSignEdge)
    if (allocated(tSignEdge))         deallocate(tSignEdge)

  end subroutine mesh_final

end module mesh_mod
