module poly_fit_mod

  use params_mod
  use mesh_mod
  use sphere_geometry_mod
  use math_mod

  implicit none

  private

  public poly_fit_init
  public poly_fit_final
  public nFit2Cells
  public nFit3Cells
  public nFit4Cells
  public fit2Cells
  public fit3Cells
  public fit4Cells
  public deriv2OnCell
  public deriv3OnCell
  public deriv4OnCell
  public nFit2Vertices
  public nFit3Vertices
  public nFit4Vertices
  public fit2Vertices
  public fit3Vertices
  public fit4Vertices
  public deriv2OnVertex
  public deriv3OnVertex
  public deriv4OnVertex

  integer, parameter :: maxFitCells = 30

  integer, allocatable :: nFit2Cells   (:)
  integer, allocatable :: nFit3Cells   (:)
  integer, allocatable :: nFit4Cells   (:)
  
  integer, allocatable :: nFit2Vertices(:)
  integer, allocatable :: nFit3Vertices(:)
  integer, allocatable :: nFit4Vertices(:)
  
  integer, allocatable :: fit2Cells    (:,:)
  integer, allocatable :: fit3Cells    (:,:)
  integer, allocatable :: fit4Cells    (:,:)
  
  integer, allocatable :: fit2Vertices (:,:)
  integer, allocatable :: fit3Vertices (:,:)
  integer, allocatable :: fit4Vertices (:,:)

  real(real_kind), allocatable :: deriv2OnCell  (:,:,:)
  real(real_kind), allocatable :: deriv3OnCell  (:,:,:)
  real(real_kind), allocatable :: deriv4OnCell  (:,:,:)
  
  real(real_kind), allocatable :: deriv2OnVertex(:,:,:)
  real(real_kind), allocatable :: deriv3OnVertex(:,:,:)
  real(real_kind), allocatable :: deriv4OnVertex(:,:,:)

contains

  subroutine poly_fit_init()

    if (.not. allocated(nFit2Cells     )) allocate(nFit2Cells    (nCells   ))
    if (.not. allocated(nFit3Cells     )) allocate(nFit3Cells    (nCells   ))
    if (.not. allocated(nFit4Cells     )) allocate(nFit4Cells    (nCells   ))
    
    if (.not. allocated(nFit2Vertices  )) allocate(nFit2Vertices (nVertices))
    if (.not. allocated(nFit3Vertices  )) allocate(nFit3Vertices (nVertices))
    if (.not. allocated(nFit4Vertices  )) allocate(nFit4Vertices (nVertices))
    
    if (.not. allocated(fit2Cells      )) allocate(fit2Cells     (0:maxFitCells,nCells   ))
    if (.not. allocated(fit3Cells      )) allocate(fit3Cells     (0:maxFitCells,nCells   ))
    if (.not. allocated(fit4Cells      )) allocate(fit4Cells     (0:maxFitCells,nCells   ))
    
    if (.not. allocated(fit2Vertices   )) allocate(fit2Vertices  (0:maxFitCells,nVertices))
    if (.not. allocated(fit3Vertices   )) allocate(fit3Vertices  (0:maxFitCells,nVertices))
    if (.not. allocated(fit4Vertices   )) allocate(fit4Vertices  (0:maxFitCells,nVertices))
    
    if (.not. allocated(deriv2OnCell   )) allocate(deriv2OnCell  (0:maxFitCells,2,nEdges ))
    if (.not. allocated(deriv3OnCell   )) allocate(deriv3OnCell  (0:maxFitCells,2,nEdges ))
    if (.not. allocated(deriv4OnCell   )) allocate(deriv4OnCell  (0:maxFitCells,2,nEdges ))
    
    if (.not. allocated(deriv2OnVertex )) allocate(deriv2OnVertex(0:maxFitCells,2,nEdges ))
    if (.not. allocated(deriv3OnVertex )) allocate(deriv3OnVertex(0:maxFitCells,2,nEdges ))
    if (.not. allocated(deriv4OnVertex )) allocate(deriv4OnVertex(0:maxFitCells,2,nEdges ))
    
    call poly_fit_run()

  end subroutine poly_fit_init

  subroutine poly_fit_final()

    if (allocated(nFit2Cells     )) deallocate(nFit2Cells    )
    if (allocated(nFit3Cells     )) deallocate(nFit3Cells    )
    if (allocated(nFit4Cells     )) deallocate(nFit4Cells    )
    
    if (allocated(nFit2Vertices  )) deallocate(nFit2Vertices )
    if (allocated(nFit3Vertices  )) deallocate(nFit3Vertices )
    if (allocated(nFit4Vertices  )) deallocate(nFit4Vertices )
    
    if (allocated(fit2Cells      )) deallocate(fit2Cells     )
    if (allocated(fit3Cells      )) deallocate(fit3Cells     )
    if (allocated(fit4Cells      )) deallocate(fit4Cells     )
    
    if (allocated(fit2Vertices   )) deallocate(fit2Vertices  )
    if (allocated(fit3Vertices   )) deallocate(fit3Vertices  )
    if (allocated(fit4Vertices   )) deallocate(fit4Vertices  )
    
    if (allocated(deriv2OnCell   )) deallocate(deriv2OnCell  )
    if (allocated(deriv3OnCell   )) deallocate(deriv3OnCell  )
    if (allocated(deriv4OnCell   )) deallocate(deriv4OnCell  )
    
    if (allocated(deriv2OnVertex )) deallocate(deriv2OnVertex)
    if (allocated(deriv3OnVertex )) deallocate(deriv3OnVertex)
    if (allocated(deriv4OnVertex )) deallocate(deriv4OnVertex)

  end subroutine poly_fit_final

  subroutine poly_fit_run()
    
    real(real_kind) xc   (0:maxFitCells)
    real(real_kind) yc   (0:maxFitCells)
    real(real_kind) zc   (0:maxFitCells)
    real(real_kind) xp   (0:maxFitCells)  ! Projected coordinate x
    real(real_kind) yp   (0:maxFitCells)  ! Projected coordinate y
    real(real_kind) theta(  maxFitCells) ! Anticlockwise angle from the first line to other lines connecting cell center and neighbor cells
    real(real_kind) d
    real(real_kind) P     (maxFitCells,maxFitCells)
    real(real_kind) B     (maxFitCells,maxFitCells)
    real(real_kind) W     (maxFitCells,maxFitCells)
    real(real_kind) PT    (maxFitCells,maxFitCells)
    real(real_kind) WTW   (maxFitCells,maxFitCells)
    real(real_kind) PTWTW (maxFitCells,maxFitCells)
    real(real_kind) PTWTWP(maxFitCells,maxFitCells)
    real(real_kind) cos_theta, sin_theta
    real(real_kind) d2fdx2, d2fdxdy , d2fdy2
    real(real_kind) d3fdx3, d3fdx2dy, d3fdxdy2 , d3fdy3
    real(real_kind) d4fdx4, d4fdx3dy, d4fdx2dy2, d4fdxdy3, d4fdy4
    integer m ! Number of involving cells except for center cell
    integer n ! Number of parameters for polynomial fit
    integer iCell, iEdge, iVertex
    integer iNgbCell  , iRecordedNgbCell
    integer iNgbVertex, iRecordedNgbVertex
    integer i, j, k
    integer iVertexOnVertex, haloVertex, haloCell
    logical recorded
    
    !!!!!!!!!!!!!!!!!!!!
    ! Fitting on cells !
    !!!!!!!!!!!!!!!!!!!!
    do iCell = 1, nCells
      ! Record cell indices for fitting.
      fit2Cells(0,iCell) = iCell
      fit4Cells(0,iCell) = iCell
      ! First halo
      do i = 1, nEdgesOnCell(iCell)
        fit2Cells(i,iCell) = cellsOnCell(i,iCell)
        fit4Cells(i,iCell) = cellsOnCell(i,iCell)
      end do
      nFit2Cells(iCell) = nEdgesOnCell(iCell) + 1
      nFit4Cells(iCell) = nEdgesOnCell(iCell) + 1
      ! Second halo
      do i = 1, nEdgesOnCell(iCell)
        do j = 1, nEdgesOnCell(fit4Cells(i,iCell))
          iNgbCell = cellsOnCell(j,fit4Cells(i,iCell))
          ! Check if cell has been recorded.
          recorded = .false.
          do k = 0, nFit4Cells(iCell) - 1
            iRecordedNgbCell = fit4Cells(k,iCell)
            if (iNgbCell == iRecordedNgbCell) then
              recorded = .true.
              exit
            end if
          end do
          if (.not. recorded) then
            fit4Cells(nFit4Cells(iCell),iCell) = iNgbCell
            nFit4Cells(iCell) = nFit4Cells(iCell) + 1
          end if
        end do
      end do
      fit3Cells (:,iCell) = fit4Cells (:,iCell)
      nFit3Cells(iCell  ) = nFit4Cells(iCell  )

      !!!!!!!!!!!!!!!!!!!!
      ! Second-order fit !
      !!!!!!!!!!!!!!!!!!!!
      ! Get Cartesian coordinates of fit cells (remove Earth radius).
      do i = 0, nFit2Cells(iCell) - 1
        xc(i) = xCell(fit2Cells(i,iCell)) / radius
        yc(i) = yCell(fit2Cells(i,iCell)) / radius
        zc(i) = zCell(fit2Cells(i,iCell)) / radius
      end do

      ! Calculate angles and projected coordinates.
      xp(0) = 0.0d0
      yp(0) = 0.0d0
      do i = 1, nFit2Cells(iCell) - 1
        theta(i) = calc_sphere_angle([xc(0),yc(0),zc(0)], &
                                     [xc(1),yc(1),zc(1)], &
                                     [xc(i),yc(i),zc(i)])
        
        d        = radius * calc_arc_length([xc(0),yc(0),zc(0)], &
                                            [xc(i),yc(i),zc(i)])
        
        xp(i) = cos(theta(i)) * d
        yp(i) = sin(theta(i)) * d
      end do

      ! Set matrices for least square fit.
      m = nFit2Cells(iCell) - 1
      n = 5
      P = 0.0d0 ! m x n
      W = 0.0d0 ! m x m
      B = 0.0d0 ! n x m
      do i = 1, m
        P(i,1) = xp(i)
        P(i,2) = yp(i)
        P(i,3) = xp(i)**2
        P(i,4) = xp(i) * yp(i)
        P(i,5) = yp(i)**2
        W(i,i) = 1.0d0
      end do

      WTW   (1:m,1:m) = matmul(transpose(W(1:m,1:m)), W(1:m,1:m))
      PT    (1:n,1:m) = transpose(P(1:m,1:n))
      PTWTW (1:n,1:m) = matmul(PT(1:n,1:m), WTW(1:m,1:m))
      PTWTWP(1:n,1:n) = matmul(PTWTW(1:n,1:m), P(1:m,1:n))
      call math_inv_matrix(n, PTWTWP(1:n,1:n), B(1:n,1:n))
      !call lapack_inv_matrix(n, PTWTWP(1:n,1:n), B(1:n,1:n))
      B     (1:n,1:m) = matmul(B(1:n,1:n), PTWTW(1:n,1:m))

      ! Calculate second-order derivative weights.
      do i = 1, nEdgesOnCell(iCell)
        iEdge = edgesOnCell(i,iCell)
        k = merge(1, 2, iCell == cellsOnEdge(1,iEdge))
        cos_theta     = cos(theta(i))
        sin_theta     = sin(theta(i))
        do j = 1, nFit2Cells(iCell) - 1
          d2fdx2  = 2.0d0 * B(3,j) * cos_theta**2
          d2fdxdy =         B(4,j) * cos_theta * sin_theta
          d2fdy2  = 2.0d0 * B(5,j) * sin_theta**2
          
          deriv2OnCell(j,k,iEdge) = d2fdx2 + 2.0d0 * d2fdxdy + d2fdy2
        end do
      end do
      
      !!!!!!!!!!!!!!!!!!!
      ! Third-order fit !
      !!!!!!!!!!!!!!!!!!!
      ! Get Cartesian coordinates of fit cells (remove Earth radius).
      do i = 0, nFit3Cells(iCell) - 1
        xc(i) = xCell(fit3Cells(i,iCell)) / radius
        yc(i) = yCell(fit3Cells(i,iCell)) / radius
        zc(i) = zCell(fit3Cells(i,iCell)) / radius
      end do

      ! Calculate angles and projected coordinates.
      xp(0) = 0.0d0
      yp(0) = 0.0d0
      do i = 1, nFit3Cells(iCell) - 1
        theta(i) = calc_sphere_angle([xc(0),yc(0),zc(0)], &
                                     [xc(1),yc(1),zc(1)], &
                                     [xc(i),yc(i),zc(i)])
        
        d        = radius * calc_arc_length([xc(0),yc(0),zc(0)], &
                                            [xc(i),yc(i),zc(i)])
        
        xp(i) = cos(theta(i)) * d
        yp(i) = sin(theta(i)) * d
      end do

      ! Set matrices for least square fit.
      m = nFit3Cells(iCell) - 1
      n = 9
      P = 0.0d0 ! m x n
      W = 0.0d0 ! m x m
      B = 0.0d0 ! n x m
      do i = 1, m
          P(i,1) = xp(i)
          P(i,2) = yp(i)
   
          P(i,3) = xp(i)**2
          P(i,4) = xp(i) * yp(i)
          P(i,5) = yp(i)**2
   
          P(i,6) = xp(i)**3
          P(i,7) = yp(i) * (xp(i)**2)
          P(i,8) = xp(i) * (yp(i)**2)
          P(i,9) = yp(i)**3
          
          W(i,i) = 1.0d0
      end do

      WTW   (1:m,1:m) = matmul(transpose(W(1:m,1:m)), W(1:m,1:m))
      PT    (1:n,1:m) = transpose(P(1:m,1:n))
      PTWTW (1:n,1:m) = matmul(PT(1:n,1:m), WTW(1:m,1:m))
      PTWTWP(1:n,1:n) = matmul(PTWTW(1:n,1:m), P(1:m,1:n))
      call math_inv_matrix(n, PTWTWP(1:n,1:n), B(1:n,1:n))
      !call lapack_inv_matrix(n, PTWTWP(1:n,1:n), B(1:n,1:n))
      B     (1:n,1:m) = matmul(B(1:n,1:n), PTWTW(1:n,1:m))

      ! Calculate second-order derivative weights.
      do i = 1, nEdgesOnCell(iCell)
        iEdge = edgesOnCell(i,iCell)
        k = merge(1, 2, iCell == cellsOnEdge(1,iEdge))
        cos_theta     = cos(theta(i))
        sin_theta     = sin(theta(i))
        do j = 1, nFit3Cells(iCell) - 1
          d3fdx3   = 6.0d0 * B(6,j) *  cos_theta**3
          d3fdx2dy = 2.0d0 * B(7,j) * (cos_theta**2) *  sin_theta
          d3fdx2dy = 2.0d0 * B(8,j) *  cos_theta     * (sin_theta)**2
          d3fdy3   = 6.0d0 * B(9,j) *  sin_theta**3
          
          deriv3OnCell(j,k,iEdge) = d3fdx3 + 3.0d0 * d3fdx2dy + 3.0d0 * d3fdxdy2 + d3fdy3
        end do
      end do

      !!!!!!!!!!!!!!!!!!!!
      ! Fourth-order fit !
      !!!!!!!!!!!!!!!!!!!!
      ! Get Cartesian coordinates of fit cells (remove Earth radius).
      do i = 0, nFit4Cells(iCell) - 1
        xc(i) = xCell(fit4Cells(i,iCell)) / radius
        yc(i) = yCell(fit4Cells(i,iCell)) / radius
        zc(i) = zCell(fit4Cells(i,iCell)) / radius
      end do

      ! Calculate angles and projected coordinates.
      xp(0) = 0.0d0
      yp(0) = 0.0d0
      do i = 1, nFit4Cells(iCell) - 1
        theta(i) = calc_sphere_angle([xc(0),yc(0),zc(0)], &
                                     [xc(1),yc(1),zc(1)], &
                                     [xc(i),yc(i),zc(i)])
        
        d        = radius * calc_arc_length([xc(0),yc(0),zc(0)], &
                                            [xc(i),yc(i),zc(i)])
        
        xp(i) = cos(theta(i)) * d
        yp(i) = sin(theta(i)) * d
      end do

      ! Set matrices for least square fit.
      m = nFit4Cells(iCell) - 1
      n = 14
      P = 0.0d0 ! m x n
      W = 0.0d0 ! m x m
      B = 0.0d0 ! n x m
      do i = 1, m
          P(i,1 ) = xp(i)
          P(i,2 ) = yp(i)
   
          P(i,3 ) = xp(i)**2
          P(i,4 ) = xp(i) * yp(i)
          P(i,5 ) = yp(i)**2
   
          P(i,6 )  = xp(i)**3
          P(i,7 )  = yp(i) * (xp(i)**2)
          P(i,8 )  = xp(i) * (yp(i)**2)
          P(i,9 ) = yp(i)**3
   
          P(i,10) = xp(i)**4
          P(i,11) = yp(i) * (xp(i)**3)
          P(i,12) = (xp(i)**2)*(yp(i)**2)
          P(i,13) = xp(i) * (yp(i)**3)
          P(i,14) = yp(i)**4
          W(i,i ) = 1.0d0
      end do

      WTW   (1:m,1:m) = matmul(transpose(W(1:m,1:m)), W(1:m,1:m))
      PT    (1:n,1:m) = transpose(P(1:m,1:n))
      PTWTW (1:n,1:m) = matmul(PT(1:n,1:m), WTW(1:m,1:m))
      PTWTWP(1:n,1:n) = matmul(PTWTW(1:n,1:m), P(1:m,1:n))
      call math_inv_matrix(n, PTWTWP(1:n,1:n), B(1:n,1:n))
      !call lapack_inv_matrix(n, PTWTWP(1:n,1:n), B(1:n,1:n))
      B     (1:n,1:m) = matmul(B(1:n,1:n), PTWTW(1:n,1:m))

      ! Calculate second-order derivative weights.
      do i = 1, nEdgesOnCell(iCell)
        iEdge = edgesOnCell(i,iCell)
        k = merge(1, 2, iCell == cellsOnEdge(1,iEdge))
        cos_theta     = cos(theta(i))
        sin_theta     = sin(theta(i))
        do j = 1, nFit4Cells(iCell) - 1
          d4fdx4    = 24.0d0 * B(10,j) *  cos_theta**4
          d4fdx3dy  =  6.0d0 * B(11,j) * (cos_theta**3) *  sin_theta
          d4fdx2dy2 =  4.0d0 * B(12,j) * (cos_theta**2) * (sin_theta**2)
          d4fdxdy3  =  6.0d0 * B(13,j) *  cos_theta     * (sin_theta**3)
          d4fdy4    = 24.0d0 * B(14,j) *  sin_theta**4
          
          deriv4OnCell(j,k,iEdge) = d4fdx4 + 4.0d0 * d4fdx3dy + 6.0d0 * d4fdx2dy2 + 4.0d0 * d4fdxdy3 + d4fdy4
        end do
      end do
    end do
    
    !!!!!!!!!!!!!!!!!!!!!!!
    ! Fitting on vertices !
    !!!!!!!!!!!!!!!!!!!!!!!
    do iVertex = 1, nVertices
      ! Record cell indices for fitting.
      fit2Vertices(0,iVertex) = iVertex
      fit3Vertices(0,iVertex) = iVertex
      fit4Vertices(0,iVertex) = iVertex
      
      ! Choose the neighbor vertices for calculate the local coordicate
      do i = 1, vertexDegree
        fit2Vertices(i,iVertex) = verticesOnVertex(i,iVertex)
      end do
      nFit2Vertices(iVertex)  = 4
      
      ! First halo
      haloVertex = iVertex
      do i = 1, vertexDegree
        haloCell = cellsOnVertex(i,haloVertex)
        do j = 1, nEdgesOnCell(haloCell)
          iNgbVertex = verticesOnCell(j,haloCell)
          
          ! Check if vertex has been recorded.
          recorded = .false.
          do k = 0, nFit2Vertices(iVertex) - 1
            iRecordedNgbVertex = fit2Vertices(k,iVertex)
            if(iNgbVertex == iRecordedNgbVertex)then
              recorded = .true.
              exit
            end if
          end do
          
          if (.not. recorded) then
            fit2Vertices(nFit2Vertices(iVertex),iVertex) = iNgbVertex
            nFit2Vertices(iVertex) = nFit2Vertices(iVertex) + 1
          end if
        end do
      end do
      nFit3Vertices(iVertex)  = nFit2Vertices(iVertex)
      nFit4Vertices(iVertex)  = nFit2Vertices(iVertex)
      fit3Vertices(:,iVertex) = fit2Vertices(:,iVertex)
      fit4Vertices(:,iVertex) = fit2Vertices(:,iVertex)
      
      ! Second halo
      do iVertexOnVertex = 1, vertexDegree
        haloVertex = verticesOnVertex(iVertexOnVertex,iVertex)
        do i = 1, vertexDegree
          haloCell = cellsOnVertex(i,haloVertex)
          do j = 1, nEdgesOnCell(haloCell)
            iNgbVertex = verticesOnCell(j,haloCell)
            
            ! Check if vertex has been recorded.
            recorded = .false.
            do k = 0, nFit4Vertices(iVertex) - 1
              iRecordedNgbVertex = fit4Vertices(k,iVertex)
              if(iNgbVertex == iRecordedNgbVertex)then
                recorded = .true.
                exit
              end if
            end do

            if (.not. recorded) then
              fit4Vertices(nFit4Vertices(iVertex),iVertex) = iNgbVertex
              nFit4Vertices(iVertex) = nFit4Vertices(iVertex) + 1
            end if
          end do
        end do
      end do

      !!!!!!!!!!!!!!!!!!!!
      ! Second-order fit !
      !!!!!!!!!!!!!!!!!!!!
      ! Get Cartesian coordinates of fit cells (remove Earth radius).
      do i = 0, nFit2Vertices(iVertex) - 1
        xc(i) = xVertex(fit2Vertices(i,iVertex)) / radius
        yc(i) = yVertex(fit2Vertices(i,iVertex)) / radius
        zc(i) = zVertex(fit2Vertices(i,iVertex)) / radius
      end do

      ! Calculate angles and projected coordinates.
      xp(0) = 0.0d0
      yp(0) = 0.0d0
      do i = 1, nFit2Vertices(iVertex) - 1
        theta(i) = calc_sphere_angle([xc(0),yc(0),zc(0)], &
                                     [xc(1),yc(1),zc(1)], &
                                     [xc(i),yc(i),zc(i)])
        
        d        = radius * calc_arc_length([xc(0),yc(0),zc(0)], &
                                            [xc(i),yc(i),zc(i)])
        
        xp(i) = cos(theta(i)) * d
        yp(i) = sin(theta(i)) * d
      end do

      ! Set matrices for least square fit.
      m = nFit2Vertices(iVertex) - 1
      n = 5
      P = 0.0d0 ! m x n
      W = 0.0d0 ! m x m
      B = 0.0d0 ! n x m
      do i = 1, m
        P(i,1) = xp(i)
        P(i,2) = yp(i)
        P(i,3) = xp(i)**2
        P(i,4) = xp(i) * yp(i)
        P(i,5) = yp(i)**2
        W(i,i) = 1.0d0
      end do

      WTW   (1:m,1:m) = matmul(transpose(W(1:m,1:m)), W(1:m,1:m))
      PT    (1:n,1:m) = transpose(P(1:m,1:n))
      PTWTW (1:n,1:m) = matmul(PT(1:n,1:m), WTW(1:m,1:m))
      PTWTWP(1:n,1:n) = matmul(PTWTW(1:n,1:m), P(1:m,1:n))
      call math_inv_matrix(n, PTWTWP(1:n,1:n), B(1:n,1:n))
      !call lapack_inv_matrix(n, PTWTWP(1:n,1:n), B(1:n,1:n))
      B     (1:n,1:m) = matmul(B(1:n,1:n), PTWTW(1:n,1:m))

      ! Calculate second-order derivative weights.
      do i = 1, vertexDegree
        iEdge = edgesOnVertex(i,iVertex)
        k = merge(1, 2, iVertex == verticesOnEdge(1,iEdge))
        
        cos_theta = cos(theta(i))
        sin_theta = sin(theta(i))
        do j = 1, nFit2Vertices(iVertex) - 1
          d2fdx2  = 2.0d0 * B(3,j) * cos_theta**2
          d2fdxdy =         B(4,j) * cos_theta * sin_theta
          d2fdy2  = 2.0d0 * B(5,j) * sin_theta**2
          
          deriv2OnVertex(j,k,iEdge) = d2fdx2 + 2.0d0 * d2fdxdy + d2fdy2
        end do
      end do
      
      !!!!!!!!!!!!!!!!!!!!
      ! Fourth-order fit !
      !!!!!!!!!!!!!!!!!!!!
      ! Get Cartesian coordinates of fit cells (remove Earth radius).
      do i = 0, nFit4Vertices(iVertex) - 1
        xc(i) = xVertex(fit4Vertices(i,iVertex)) / radius
        yc(i) = yVertex(fit4Vertices(i,iVertex)) / radius
        zc(i) = zVertex(fit4Vertices(i,iVertex)) / radius
      end do

      ! Calculate angles and projected coordinates.
      xp(0) = 0.0d0
      yp(0) = 0.0d0
      do i = 1, nFit4Vertices(iVertex) - 1
        theta(i) = calc_sphere_angle([xc(0),yc(0),zc(0)], &
                                     [xc(1),yc(1),zc(1)], &
                                     [xc(i),yc(i),zc(i)])
        
        d        = radius * calc_arc_length([xc(0),yc(0),zc(0)], &
                                            [xc(i),yc(i),zc(i)])
        
        xp(i) = cos(theta(i)) * d
        yp(i) = sin(theta(i)) * d
      end do

      ! Set matrices for least square fit.
      m = nFit4Vertices(iVertex) - 1
      n = 14
      P = 0.0d0 ! m x n
      W = 0.0d0 ! m x m
      B = 0.0d0 ! n x m
      do i = 1, m
          P(i,1 ) = xp(i)
          P(i,2 ) = yp(i)
   
          P(i,3 ) = xp(i)**2
          P(i,4 ) = xp(i) * yp(i)
          P(i,5 ) = yp(i)**2
   
          P(i,6 )  = xp(i)**3
          P(i,7 )  = yp(i) * (xp(i)**2)
          P(i,8 )  = xp(i) * (yp(i)**2)
          P(i,9 ) = yp(i)**3
   
          P(i,10) = xp(i)**4
          P(i,11) = yp(i) * (xp(i)**3)
          P(i,12) = (xp(i)**2)*(yp(i)**2)
          P(i,13) = xp(i) * (yp(i)**3)
          P(i,14) = yp(i)**4
          W(i,i ) = 1.0d0
      end do

      WTW   (1:m,1:m) = matmul(transpose(W(1:m,1:m)), W(1:m,1:m))
      PT    (1:n,1:m) = transpose(P(1:m,1:n))
      PTWTW (1:n,1:m) = matmul(PT(1:n,1:m), WTW(1:m,1:m))
      PTWTWP(1:n,1:n) = matmul(PTWTW(1:n,1:m), P(1:m,1:n))
      call math_inv_matrix(n, PTWTWP(1:n,1:n), B(1:n,1:n))
      !call lapack_inv_matrix(n, PTWTWP(1:n,1:n), B(1:n,1:n))
      B     (1:n,1:m) = matmul(B(1:n,1:n), PTWTW(1:n,1:m))

      ! Calculate second-order derivative weights.
      do i = 1, vertexDegree
        iEdge = edgesOnVertex(i,iVertex)
        k = merge(1, 2, iVertex == verticesOnEdge(1,iEdge))
        
        cos_theta = cos(theta(i))
        sin_theta = sin(theta(i))
        do j = 1, nFit4Vertices(iVertex) - 1
          d4fdx4    = 24.0d0 * B(10,j) *  cos_theta**4
          d4fdx3dy  =  6.0d0 * B(11,j) * (cos_theta**3) *  sin_theta
          d4fdx2dy2 =  4.0d0 * B(12,j) * (cos_theta**2) * (sin_theta**2)
          d4fdxdy3  =  6.0d0 * B(13,j) *  cos_theta     * (sin_theta**3)
          d4fdy4    = 24.0d0 * B(14,j) *  sin_theta**4
          
          deriv4OnVertex(j,k,iEdge) = d4fdx4 + 4.0d0 * d4fdx3dy + 6.0d0 * d4fdx2dy2 + 4.0d0 * d4fdxdy3 + d4fdy4
        end do
      end do
      
    end do ! end do loop iVertex = 1, nVertices
    
  end subroutine poly_fit_run

end module poly_fit_mod
