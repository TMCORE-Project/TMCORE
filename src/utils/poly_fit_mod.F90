module poly_fit_mod

  use params_mod
  use mesh_mod
  use sphere_geometry_mod
  use math_mod

  implicit none

  private

  public poly_fit_init
  public poly_fit_final
  public nFitCellsOnCell
  public fitCellsOnCell
  public derivOnCell

  integer, parameter :: maxFitCells = 25
  integer, parameter :: maxFitOrder = 4

  integer, allocatable :: nFitCellsOnCell     (:,:)
  integer, allocatable :: nFitVerticesOnVertex(:,:)
  integer, allocatable :: fitCellsOnCell      (:,:,:)
  integer, allocatable :: fitVerticesOnVertex (:,:,:)

  real(real_kind), allocatable :: derivOnCell  (:,:,:,:)
  real(real_kind), allocatable :: derivOnVertex(:,:,:,:)

contains

  subroutine poly_fit_init()

    if (.not. allocated(nFitCellsOnCell     )) allocate(nFitCellsOnCell     (                2:maxFitOrder,nCells   ))
    if (.not. allocated(nFitVerticesOnVertex)) allocate(nFitVerticesOnVertex(                2:maxFitOrder,nVertices))
    if (.not. allocated(fitCellsOnCell      )) allocate(fitCellsOnCell      (0:maxFitCells,  2:maxFitOrder,nCells   ))
    if (.not. allocated(fitVerticesOnVertex )) allocate(fitVerticesOnVertex (0:maxFitCells,  2:maxFitOrder,nVertices))
    if (.not. allocated(derivOnCell         )) allocate(derivOnCell         (0:maxFitCells,2,2:maxFitOrder,nEdges   ))
    if (.not. allocated(derivOnVertex       )) allocate(derivOnVertex       (0:maxFitCells,2,2:maxFitOrder,nEdges   ))
    call poly_fit_run()

  end subroutine poly_fit_init

  subroutine poly_fit_final()

    if (allocated(nFitCellsOnCell     )) deallocate(nFitCellsOnCell     )
    if (allocated(nFitVerticesOnVertex)) deallocate(nFitVerticesOnVertex)
    if (allocated(fitCellsOnCell      )) deallocate(fitCellsOnCell      )
    if (allocated(fitVerticesOnVertex )) deallocate(fitVerticesOnVertex )
    if (allocated(derivOnCell         )) deallocate(derivOnCell         )
    if (allocated(derivOnVertex       )) deallocate(derivOnVertex       )

  end subroutine poly_fit_final

  subroutine poly_fit_run()
    implicit none
    
    real(real_kind) xc(0:maxFitCells)
    real(real_kind) yc(0:maxFitCells)
    real(real_kind) zc(0:maxFitCells)
    real(real_kind) xp(0:maxFitCells)  ! Projected coordinate x
    real(real_kind) yp(0:maxFitCells)  ! Projected coordinate y
    real(real_kind) theta(maxFitCells) ! Anticlockwise angle from the first line to other lines connecting cell center and neighbor cells
    real(real_kind) d
    real(real_kind) P     (maxFitCells,maxFitCells)
    real(real_kind) B     (maxFitCells,maxFitCells)
    real(real_kind) W     (maxFitCells,maxFitCells)
    real(real_kind) PT    (maxFitCells,maxFitCells)
    real(real_kind) WTW   (maxFitCells,maxFitCells)
    real(real_kind) PTWTW (maxFitCells,maxFitCells)
    real(real_kind) PTWTWP(maxFitCells,maxFitCells)
    real(real_kind) cos_theta, sin_theta
    real(real_kind) d2fdx2, d2fdxdy, d2fdy2
    real(real_kind) d4fdx4, d4fdx3dy, d4fdx2dy2, d4fdxdy3, d4fdy4
    integer m ! Number of involving cells except for center cell
    integer n ! Number of parameters for polynomial fit
    integer iCell, iNgbCell, iRecordedNgbCell, iEdge, iOrder, i, j, k
    logical recorded

    ! Fitting on cells
    do iCell = 1, nCells
      ! Record cell indices for fitting.
      fitCellsOnCell(0,2:4,iCell) = iCell
      ! First halo
      do i = 1, nEdgesOnCell(iCell)
        fitCellsOnCell(i,2:4,iCell) = cellsOnCell(i,iCell)
      end do
      nFitCellsOnCell(2:4,iCell) = nEdgesOnCell(iCell) + 1
      ! Second halo
      do i = 1, nEdgesOnCell(iCell)
        do j = 1, nEdgesOnCell(fitCellsOnCell(i,4,iCell))
          iNgbCell = cellsOnCell(j,fitCellsOnCell(i,4,iCell))
          ! Check if cell has been recorded.
          recorded = .false.
          do k = 0, nFitCellsOnCell(4,iCell) - 1
            iRecordedNgbCell = fitCellsOnCell(k,4,iCell)
            if (iNgbCell == iRecordedNgbCell) then
              recorded = .true.
              exit
            end if
          end do
          if (.not. recorded) then
            fitCellsOnCell(nFitCellsOnCell(4,iCell),3:4,iCell) = iNgbCell
            nFitCellsOnCell(3:4,iCell) = nFitCellsOnCell(3:4,iCell) + 1
          end if
        end do
      end do

      ! Second-order fit
      ! Get Cartesian coordinates of fit cells (remove Earth radius).
      do i = 0, nFitCellsOnCell(2,iCell) - 1
        xc(i) = xCell(fitCellsOnCell(i,2,iCell)) / radius
        yc(i) = yCell(fitCellsOnCell(i,2,iCell)) / radius
        zc(i) = zCell(fitCellsOnCell(i,2,iCell)) / radius
      end do

      ! Calculate angles and projected coordinates.
      xp(0) = 0.0d0
      yp(0) = 0.0d0
      do i = 1, nFitCellsOnCell(2,iCell) - 1
        theta(i) = calc_sphere_angle([xc(0),yc(0),zc(0)], &
                                     [xc(1),yc(1),zc(1)], &
                                     [xc(i),yc(i),zc(i)])
        
        d        = radius * calc_arc_length([xc(0),yc(0),zc(0)], &
                                            [xc(i),yc(i),zc(i)])
        
        xp(i) = cos(theta(i)) * d
        yp(i) = sin(theta(i)) * d
      end do

      ! Set matrices for least square fit.
      m = nFitCellsOnCell(2,iCell) - 1
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
        do j = 1, nFitCellsOnCell(2,iCell) - 1
          d2fdx2  = 2.0d0 * B(3,j) * cos_theta**2
          d2fdxdy =         B(4,j) * cos_theta * sin_theta
          d2fdy2  = 2.0d0 * B(5,j) * sin_theta**2
          derivOnCell(j,k,2,iEdge) = d2fdx2 + 2.0d0 * d2fdxdy + d2fdy2
        end do
      end do

      ! Fourth-order fit
      ! Get Cartesian coordinates of fit cells (remove Earth radius).
      do i = 0, nFitCellsOnCell(4,iCell) - 1
        xc(i) = xCell(fitCellsOnCell(i,4,iCell)) / radius
        yc(i) = yCell(fitCellsOnCell(i,4,iCell)) / radius
        zc(i) = zCell(fitCellsOnCell(i,4,iCell)) / radius
      end do

      ! Calculate angles and projected coordinates.
      xp(0) = 0.0d0
      yp(0) = 0.0d0
      do i = 1, nFitCellsOnCell(4,iCell) - 1
        theta(i) = calc_sphere_angle([xc(0),yc(0),zc(0)], &
                                     [xc(1),yc(1),zc(1)], &
                                     [xc(i),yc(i),zc(i)])
        
        d        = radius * calc_arc_length([xc(0),yc(0),zc(0)], &
                                            [xc(i),yc(i),zc(i)])
        
        xp(i) = cos(theta(i)) * d
        yp(i) = sin(theta(i)) * d
      end do

      ! Set matrices for least square fit.
      m = nFitCellsOnCell(4,iCell) - 1
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
        do j = 1, nFitCellsOnCell(4,iCell) - 1
          d4fdx4    = 24.0d0 * B(10,j) *  cos_theta**4
          d4fdx3dy  =  6.0d0 * B(11,j) * (cos_theta**3) *  sin_theta
          d4fdx2dy2 =  4.0d0 * B(12,j) * (cos_theta**2) * (sin_theta**2)
          d4fdxdy3  =  6.0d0 * B(13,j) *  cos_theta     * (sin_theta**3)
          d4fdy4    = 24.0d0 * B(14,j) *  sin_theta**4
          derivOnCell(j,k,4,iEdge) = d4fdx4 + 4.0d0 * d4fdx3dy + 6.0d0 * d4fdx2dy2 + 4.0d0 * d4fdxdy3 + d4fdy4
        end do
      end do
    end do
    
  end subroutine poly_fit_run

end module poly_fit_mod