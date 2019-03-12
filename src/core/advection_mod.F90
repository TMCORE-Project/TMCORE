! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module advection_mod
   
   use params_mod
   use mesh_mod
   use const_mod, sphere_radius => radius
   
   implicit none
   
   real(real_kind), dimension(:,:,:), allocatable :: deriv_two   ! 2nd order derivative for high order advection
   real(real_kind), dimension(:,:,:), allocatable :: deriv_four  ! 4nd order derivative for high order advection
   integer,         dimension(:    ), allocatable :: nFitCells3
   integer,         dimension(:    ), allocatable :: nFitCells5
   integer,         dimension(:,:  ), allocatable :: adv3Cells
   integer,         dimension(:,:  ), allocatable :: adv5Cells
   
   real(real_kind), dimension(:,:  ), allocatable :: defc_a, defc_b
   
   contains

   subroutine advecion_init
     implicit none
     
     if(adv_order>=3)then
         
         call calc_deriv_two
         
         if(adv_order>=5)then
             call calc_deriv_four
         endif
         
     endif
     
   end subroutine advecion_init

   subroutine calc_deriv_two
                                      
!
! compute the cell coefficients for the polynomial fit.
! this is performed during setup for model integration.
! WCS, 31 August 2009
!
      implicit none


!  local variables
      integer                                     , parameter   :: maxCellNum       = 25
      integer                                     , parameter   :: polynomial_order = 2  ! For now polynomial_order = 2 is the best choice, 3 or 4 will lead to bad results
      real (kind=real_kind), dimension(:  )       , allocatable :: theta_abs

      real (kind=real_kind), dimension(  maxCellNum)            :: xc, yc, zc ! cell center coordinates
      real (kind=real_kind), dimension(  maxCellNum)            :: thetav, thetaCell, dl_sphere
      real (kind=real_kind), dimension(0:maxCellNum)            :: xp, yp
      real (kind=real_kind)                                     :: xec, yec, zec
      real (kind=real_kind)                                     :: thetae_tmp
      integer                                                   :: i, j, k, ip1, ip2
      integer                                                   :: nFitCells
      integer                                                   :: iCell, iEdge
      integer                                                   :: indexCOE
      real (kind=real_kind)                                     :: d2fdx2, d2fdxdy, d2fdy2
      
      real (kind=real_kind) :: amatrix(maxCellNum,maxCellNum), &
                               bmatrix(maxCellNum,maxCellNum), &
                               wmatrix(maxCellNum,maxCellNum)
      
      integer :: polyParamNum, cell_add
      integer, dimension(maxCellNum) :: cell_list

      integer :: cell1, cell2
      logical :: add_the_cell, do_the_cell
      
      real (kind=real_kind) :: cost , sint ,                                 &
                               cost2, sint2, costsint ,                      &
                               cost3, sint3, cost2sint, costsint2 ,          &
                               cost4, sint4, cost3sint, cost2sint2, costsint3

!----------------------------------------------------------------------------
      
      polyParamNum = (polynomial_order+1)*(polynomial_order+2)/2
      
      allocate(deriv_two (maxCellNum,2,nEdges))
      allocate(adv3Cells (maxCellNum,nCells))
      allocate(nFitCells3(nCells))
      allocate(theta_abs(nCells))
      
      deriv_two(:,:,:) = 0.d0

      do iCell = 1, nCells !  is this correct? - we need first halo cell also...

        cell_list(1) = iCell
        do i = 2, nEdgesOnCell(iCell)+1
          cell_list(i) = cellsOnCell(i-1,iCell)
        end do
        nFitCells3(iCell) = nEdgesOnCell(iCell) + 1

        if ( polynomial_order > 2 ) then
          do i = 2, nEdgesOnCell(iCell) + 1
            do j = 1, nEdgesOnCell( cell_list(i) )
              cell_add = CellsOnCell(j,cell_list(i))
              add_the_cell = .true.
              
              do k=1,nFitCells3(iCell)
                if ( cell_add == cell_list(k) ) add_the_cell = .false.
              end do
              
              if (add_the_cell) then
                 nFitCells3(iCell) = nFitCells3(iCell)+1
                 cell_list(nFitCells3(iCell)) = cell_add
              end if
              
            end do
          end do
        end if

!  check to see if we are reaching outside the halo

        do_the_cell = .true.
        do i = 1, nFitCells3(iCell)
          if (cell_list(i) > nCells) do_the_cell = .false.
        end do

        if ( .not. do_the_cell ) cycle

!  compute polynomial fit for this cell if all needed neighbors exist

        do i = 1, nFitCells3(iCell)
          adv3Cells(i,iCell) = cell_list(i)
          
          xc(i) = xCell(adv3Cells(i,iCell)) / sphere_radius
          yc(i) = yCell(adv3Cells(i,iCell)) / sphere_radius
          zc(i) = zCell(adv3Cells(i,iCell)) / sphere_radius
        end do

        theta_abs(iCell) = angleEdge(edgesOnCell(1,iCell))

! angles from cell center to neighbor centers (thetav)

        do i = 1, nFitCells3(iCell)-1
   
          ip2 = i+2
          if (ip2 > nFitCells3(iCell)) ip2 = 2
    
          thetav(i) = sphere_angle( xc(1  ), yc(1  ), zc(1  ),&
                                    xc(i+1), yc(i+1), zc(i+1),&
                                    xc(ip2), yc(ip2), zc(ip2) )

          dl_sphere(i) = sphere_radius*arc_length( xc(1  ), yc(1  ), zc(1  ),&
                                                   xc(i+1), yc(i+1), zc(i+1) )
        end do

        thetaCell(1) = theta_abs(iCell)  !  this defines the x direction, longitude line
        do i=2,nFitCells3(iCell)-1
          thetaCell(i) = thetaCell(i-1) + thetav(i-1)
        end do

        xp(0) = 0.d0
        yp(0) = 0.d0
        do i=1,nFitCells3(iCell)-1
          xp(i) = cos(thetaCell(i)) * dl_sphere(i)
          yp(i) = sin(thetaCell(i)) * dl_sphere(i)
        end do

        bmatrix = 0.d0
        amatrix = 0.d0
        wmatrix = 0.d0
  
        do i=1,nFitCells3(iCell)
          amatrix(i,1) = 1.d0
          amatrix(i,2) = xp(i-1)
          amatrix(i,3) = yp(i-1)
          amatrix(i,4) = xp(i-1)**2
          amatrix(i,5) = xp(i-1) * yp(i-1)
          amatrix(i,6) = yp(i-1)**2
   
          wmatrix(i,i) = 1.d0
        end do
        
        call sw_poly_fit_2( amatrix, bmatrix, wmatrix, nFitCells3(iCell), polyParamNum, maxCellNum )

!  fill second derivative stencil for rk advection 

        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i,iCell)
  
          if     (iCell == cellsOnEdge(1,iEdge)) then
            indexCOE = 1
          elseif (iCell == cellsOnEdge(2,iEdge)) then
            indexCOE = 2
          endif
          
          cost     = cos(angleEdge(iEdge))
          sint     = sin(angleEdge(iEdge))
          costsint = cost*sint
          cost2    = cost**2
          sint2    = sint**2
   
          do j = 1, nFitCells3(iCell)
            d2fdx2  = 2.d0 * bmatrix(4 ,j) * cost2
            
            d2fdxdy =        bmatrix(5 ,j) * costsint
            
            d2fdy2  = 2.d0 * bmatrix(6 ,j) * sint2
            
            deriv_two(j,indexCOE,iEdge) = d2fdx2 + 2.d0*d2fdxdy + d2fdy2
          end do
        end do
 
      end do ! end of loop over cells
      
   end subroutine calc_deriv_two
   

   subroutine calc_deriv_four
                                      
!
! compute the cell coefficients for the polynomial fit.
! this is performed during setup for model integration.
! WCS, 31 August 2009
!
      implicit none

!  local variables
      integer                                     , parameter   :: maxCellNum       = 25
      integer                                     , parameter   :: polynomial_order = 4
      real (kind=real_kind), dimension(:  )       , allocatable :: theta_abs

      real (kind=real_kind), dimension(  maxCellNum)            :: xc, yc, zc ! cell center coordinates
      real (kind=real_kind), dimension(  maxCellNum)            :: thetav, thetaCell, dl_sphere
      real (kind=real_kind), dimension(0:maxCellNum)            :: xp, yp
      real (kind=real_kind)                                     :: xec, yec, zec
      real (kind=real_kind)                                     :: thetae_tmp
      integer                                                   :: i, j, k, ip1, ip2
      integer                                                   :: nFitCells
      integer                                                   :: iCell, iEdge
      integer                                                   :: indexCOE
      real (kind=real_kind)                                     :: d4fdx4, d4fdx3dy, d4fdx2dy2, d4fdxdy3, d4fdy4
      
      real (kind=real_kind) :: amatrix(maxCellNum,maxCellNum), &
                               bmatrix(maxCellNum,maxCellNum), &
                               wmatrix(maxCellNum,maxCellNum)
      
      integer :: polyParamNum, cell_add
      integer, dimension(maxCellNum) :: cell_list

      integer :: cell1, cell2
      logical :: add_the_cell, do_the_cell
      
      real (kind=real_kind) :: cost , sint ,                                 &
                               cost2, sint2, costsint ,                      &
                               cost3, sint3, cost2sint, costsint2 ,          &
                               cost4, sint4, cost3sint, cost2sint2, costsint3

!----------------------------------------------------------------------------
      polyParamNum = (polynomial_order+1)*(polynomial_order+2)/2
      
      allocate(deriv_four(maxCellNum,2,nEdges))
      allocate(adv5Cells (maxCellNum,nCells))
      allocate(nFitCells5(nCells))
      allocate(theta_abs(nCells))
      
      deriv_four(:,:,:) = 0.d0

      do iCell = 1, nCells !  is this correct? - we need first halo cell also...

        cell_list(1) = iCell
        do i = 2, nEdgesOnCell(iCell)+1
          cell_list(i) = cellsOnCell(i-1,iCell)
        end do
        nFitCells5(iCell) = nEdgesOnCell(iCell) + 1

        if ( polynomial_order > 2 ) then
          do i = 2, nEdgesOnCell(iCell) + 1
            do j = 1, nEdgesOnCell( cell_list(i) )
              cell_add = CellsOnCell(j,cell_list(i))
              add_the_cell = .true.
              
              do k=1,nFitCells5(iCell)
                if ( cell_add == cell_list(k) ) add_the_cell = .false.
              end do
              
              if (add_the_cell) then
                 nFitCells5(iCell) = nFitCells5(iCell)+1
                 cell_list(nFitCells5(iCell)) = cell_add
              end if
              
            end do
          end do
        end if

!  check to see if we are reaching outside the halo

        do_the_cell = .true.
        do i = 1, nFitCells5(iCell)
           if (cell_list(i) > nCells) do_the_cell = .false.
        end do

        if ( .not. do_the_cell ) cycle

!  compute polynomial fit for this cell if all needed neighbors exist

        do i = 1, nFitCells5(iCell)
          adv5Cells(i,iCell) = cell_list(i)
          
          xc(i) = xCell(adv5Cells(i,iCell)) / sphere_radius
          yc(i) = yCell(adv5Cells(i,iCell)) / sphere_radius
          zc(i) = zCell(adv5Cells(i,iCell)) / sphere_radius
        end do

        theta_abs(iCell) = angleEdge(edgesOnCell(1,iCell))

! angles from cell center to neighbor centers (thetav)

        do i = 1, nFitCells5(iCell)-1
   
          ip2 = i+2
          if (ip2 > nFitCells5(iCell)) ip2 = 2
    
          thetav(i) = sphere_angle( xc(1  ), yc(1  ), zc(1  ),&
                                    xc(i+1), yc(i+1), zc(i+1),&
                                    xc(ip2), yc(ip2), zc(ip2) )

          dl_sphere(i) = sphere_radius * arc_length( xc(1  ), yc(1  ), zc(1  ),&
                                                     xc(i+1), yc(i+1), zc(i+1) )
        end do

        thetaCell(1) = theta_abs(iCell)  !  this defines the x direction, longitude line
        do i=2,nFitCells5(iCell)-1
          thetaCell(i) = thetaCell(i-1) + thetav(i-1)
        end do
   
        xp(0) = 0.d0
        yp(0) = 0.d0
        do i=1,nFitCells5(iCell)-1
          xp(i) = cos(thetaCell(i)) * dl_sphere(i)
          yp(i) = sin(thetaCell(i)) * dl_sphere(i)
        end do

        bmatrix = 0.d0
        amatrix = 0.d0
        wmatrix = 0.d0
  
        do i=1,nFitCells5(iCell)
          amatrix(i,1) = 1.d0
          amatrix(i,2) = xp(i-1)
          amatrix(i,3) = yp(i-1)
   
          amatrix(i,4) = xp(i-1)**2
          amatrix(i,5) = xp(i-1) * yp(i-1)
          amatrix(i,6) = yp(i-1)**2
   
          amatrix(i,7)  = xp(i-1)**3
          amatrix(i,8)  = yp(i-1) * (xp(i-1)**2)
          amatrix(i,9)  = xp(i-1) * (yp(i-1)**2)
          amatrix(i,10) = yp(i-1)**3
   
          amatrix(i,11) = xp(i-1)**4
          amatrix(i,12) = yp(i-1) * (xp(i-1)**3)
          amatrix(i,13) = (xp(i-1)**2)*(yp(i-1)**2)
          amatrix(i,14) = xp(i-1) * (yp(i-1)**3)
          amatrix(i,15) = yp(i-1)**4
   
          wmatrix(i,i) = 1.d0
        end do
          
        call sw_poly_fit_2( amatrix, bmatrix, wmatrix, nFitCells5(iCell), polyParamNum, maxCellNum )

!  fill fourth derivative stencil for rk advection 

        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i,iCell)
  
          if     (iCell == cellsOnEdge(1,iEdge)) then
            indexCOE = 1
          elseif (iCell == cellsOnEdge(2,iEdge)) then
            indexCOE = 2
          endif
          
          cost       = cos(angleEdge(iEdge))
          sint       = sin(angleEdge(iEdge))
          cost4      = cost**4
          sint4      = sint**4
          cost3sint  = cost**3 * sint
          cost2sint2 = cost**2 * sint**2
          costsint3  = cost    * sint**3
   
          do j = 1, nFitCells5(iCell)
              
            d4fdx4    = 24.d0 * bmatrix(11 ,j) * cost4
            
            d4fdx3dy  =  6.d0 * bmatrix(12 ,j) * cost3sint
            
            d4fdx2dy2 =  4.d0 * bmatrix(13 ,j) * cost2sint2
            
            d4fdxdy3  =  6.d0 * bmatrix(14 ,j) * costsint3
            
            d4fdy4    = 24.d0 * bmatrix(15 ,j) * sint4
            
            deriv_four(j,indexCOE,iEdge) = d4fdx4 + 4.d0*d4fdx3dy + 6.d0*d4fdx2dy2 + 4.d0*d4fdxdy3 + d4fdy4
            
          end do
        end do ! nEdgesOnCell(iCell)
 
      end do ! end of loop over cells
      
   end subroutine calc_deriv_four

   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION SPHERE_ANGLE
   !
   ! Computes the angle between arcs AB and AC, given points A, B, and C
   ! Equation numbers w.r.t. http://mathworld.wolfram.com/SphericalTrigonometry.html
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (kind=real_kind) function sphere_angle(ax, ay, az, bx, by, bz, cx, cy, cz)
   
      implicit none
   
      real (kind=real_kind), intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz
   
      real (kind=real_kind) :: a, b, c          ! Side lengths of spherical triangle ABC
   
      real (kind=real_kind) :: ABx, ABy, ABz    ! The components of the vector AB
      real (kind=real_kind) :: mAB              ! The magnitude of AB
      real (kind=real_kind) :: ACx, ACy, ACz    ! The components of the vector AC
      real (kind=real_kind) :: mAC              ! The magnitude of AC
   
      real (kind=real_kind) :: Dx               ! The i-components of the cross product AB x AC
      real (kind=real_kind) :: Dy               ! The j-components of the cross product AB x AC
      real (kind=real_kind) :: Dz               ! The k-components of the cross product AB x AC
   
      real (kind=real_kind) :: s                ! Semiperimeter of the triangle
      real (kind=real_kind) :: sin_angle
   
      a = acos(max(min(bx*cx + by*cy + bz*cz,1.0_real_kind),-1.0_real_kind))      ! Eqn. (3)
      b = acos(max(min(ax*cx + ay*cy + az*cz,1.0_real_kind),-1.0_real_kind))      ! Eqn. (2)
      c = acos(max(min(ax*bx + ay*by + az*bz,1.0_real_kind),-1.0_real_kind))      ! Eqn. (1)
   
      ABx = bx - ax
      ABy = by - ay
      ABz = bz - az
   
      ACx = cx - ax
      ACy = cy - ay
      ACz = cz - az
   
      Dx =   (ABy * ACz) - (ABz * ACy)
      Dy = -((ABx * ACz) - (ABz * ACx))
      Dz =   (ABx * ACy) - (ABy * ACx)
   
      s = 0.5*(a + b + c)
!      sin_angle = sqrt((sin(s-b)*sin(s-c))/(sin(b)*sin(c)))   ! Eqn. (28)
      sin_angle = sqrt(min(1.0_real_kind,max(0.0_real_kind,(sin(s-b)*sin(s-c))/(sin(b)*sin(c)))))   ! Eqn. (28)
   
      if ((Dx*ax + Dy*ay + Dz*az) >= 0.0) then
         sphere_angle =  2.0 * asin(max(min(sin_angle,1.0_real_kind),-1.0_real_kind))
      else
         sphere_angle = -2.0 * asin(max(min(sin_angle,1.0_real_kind),-1.0_real_kind))
      end if
   
   end function sphere_angle
   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION PLANE_ANGLE
   !
   ! Computes the angle between vectors AB and AC, given points A, B, and C, and
   !   a vector (u,v,w) normal to the plane.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (kind=real_kind) function plane_angle(ax, ay, az, bx, by, bz, cx, cy, cz, u, v, w)
   
      implicit none
   
      real (kind=real_kind), intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz, u, v, w
   
      real (kind=real_kind) :: ABx, ABy, ABz    ! The components of the vector AB
      real (kind=real_kind) :: mAB              ! The magnitude of AB
      real (kind=real_kind) :: ACx, ACy, ACz    ! The components of the vector AC
      real (kind=real_kind) :: mAC              ! The magnitude of AC
   
      real (kind=real_kind) :: Dx               ! The i-components of the cross product AB x AC
      real (kind=real_kind) :: Dy               ! The j-components of the cross product AB x AC
      real (kind=real_kind) :: Dz               ! The k-components of the cross product AB x AC
   
      real (kind=real_kind) :: cos_angle
   
      ABx = bx - ax
      ABy = by - ay
      ABz = bz - az
      mAB = sqrt(ABx**2.0 + ABy**2.0 + ABz**2.0)
   
      ACx = cx - ax
      ACy = cy - ay
      ACz = cz - az
      mAC = sqrt(ACx**2.0 + ACy**2.0 + ACz**2.0)
   
   
      Dx =   (ABy * ACz) - (ABz * ACy)
      Dy = -((ABx * ACz) - (ABz * ACx))
      Dz =   (ABx * ACy) - (ABy * ACx)
   
      cos_angle = (ABx*ACx + ABy*ACy + ABz*ACz) / (mAB * mAC)
   
      if ((Dx*u + Dy*v + Dz*w) >= 0.0) then
         plane_angle =  acos(max(min(cos_angle,1.0_real_kind),-1.0_real_kind))
      else
         plane_angle = -acos(max(min(cos_angle,1.0_real_kind),-1.0_real_kind))
      end if
   
   end function plane_angle


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION ARC_LENGTH
   !
   ! Returns the length of the great circle arc from A=(ax, ay, az) to 
   !    B=(bx, by, bz). It is assumed that both A and B lie on the surface of the
   !    same sphere centered at the origin.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (kind=real_kind) function arc_length(ax, ay, az, bx, by, bz)
   
      implicit none
   
      real (kind=real_kind), intent(in) :: ax, ay, az, bx, by, bz
   
      real (kind=real_kind) :: r, c
      real (kind=real_kind) :: cx, cy, cz
   
      cx = bx - ax
      cy = by - ay
      cz = bz - az

!      r = ax*ax + ay*ay + az*az
!      c = cx*cx + cy*cy + cz*cz
!
!      arc_length = sqrt(r) * acos(1.0 - c/(2.0*r))

      r = sqrt(ax*ax + ay*ay + az*az)
      c = sqrt(cx*cx + cy*cy + cz*cz)
!      arc_length = sqrt(r) * 2.0 * asin(c/(2.0*r))
      arc_length = r * 2.0 * asin(c/(2.0*r))

   end function arc_length
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! subroutine sw_arc_bisect
   !
   ! Returns the point C=(cx, cy, cz) that bisects the great circle arc from
   !   A=(ax, ay, az) to B=(bx, by, bz). It is assumed that A and B lie on the
   !   surface of a sphere centered at the origin.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine sw_arc_bisect(ax, ay, az, bx, by, bz, cx, cy, cz)
   
      implicit none
   
      real (kind=real_kind), intent(in) :: ax, ay, az, bx, by, bz
      real (kind=real_kind), intent(out) :: cx, cy, cz
   
      real (kind=real_kind) :: r           ! Radius of the sphere
      real (kind=real_kind) :: d           
   
      r = sqrt(ax*ax + ay*ay + az*az)
   
      cx = 0.5d0*(ax + bx)
      cy = 0.5d0*(ay + by)
      cz = 0.5d0*(az + bz)
   
      if (cx == 0.d0 .and. cy == 0.d0 .and. cz == 0.d0) then
         stop 'arc_bisect: A and B are diametrically opposite'
      else
         d = sqrt(cx*cx + cy*cy + cz*cz)
         cx = r * cx / d
         cy = r * cy / d
         cz = r * cz / d
      end if
   
   end subroutine sw_arc_bisect


   subroutine sw_poly_fit_2(aMatrix_input, bMatrix_output, weightsMatrix_input, cellNum, polyParamNum, maxCellNum)

     implicit none

     integer                                                , intent(in ) :: cellNum       , polyParamNum       , maxCellNum
     real(kind=real_kind), dimension(maxCellNum, maxCellNum), intent(in ) :: aMatrix_input , weightsMatrix_input
     real(kind=real_kind), dimension(maxCellNum, maxCellNum), intent(out) :: bMatrix_output
   
     ! local storage
   
     real (kind=real_kind), dimension(cellNum     , polyParamNum) :: aMatrix
     real (kind=real_kind), dimension(polyParamNum, cellNum     ) :: bMatrix
     
     real (kind=real_kind), dimension(cellNum     , cellNum     ) :: weightsMatrix, weightsTransposeMatrix, wTw  ! Working arrays
     real (kind=real_kind), dimension(polyParamNum, cellNum     ) :: aTransposeMatrix, aTh                       ! Working arrays
     real (kind=real_kind), dimension(polyParamNum, polyParamNum) :: aTa, ata_inv, aTha, atha_inv                ! Working arrays
     integer              , dimension(polyParamNum              ) :: indx
   
     if ( (maxCellNum<polyParamNum) .or. (maxCellNum<cellNum) ) then
       stop 'error in poly_fit_2 inversion, set a larger maxCellNum'
     end if
     
     aMatrix       (1:cellNum, 1:polyParamNum) = aMatrix_input      (1:cellNum,1:polyParamNum)
     weightsMatrix (1:cellNum, 1:cellNum     ) = weightsMatrix_input(1:cellNum,1:cellNum     ) 
     bMatrix_output                            = 0.d0

     weightsTransposeMatrix = transpose(weightsMatrix                        )
     wTw                    = matmul   (weightsTransposeMatrix, weightsMatrix)
     aTransposeMatrix       = transpose(aMatrix                              )
     aTh                    = matmul   (aTransposeMatrix      , wTw          )
     aTha                   = matmul   (aTh                   , aMatrix      )
     aTa                    = matmul   (aTransposeMatrix      , aMatrix      )
     
     call sw_migs(aTha, polyParamNum, atha_inv, indx)

     bMatrix = matmul(atha_inv, aTh)

     bMatrix_output(1:polyParamNum,1:cellNum) = bMatrix(1:polyParamNum,1:cellNum)

   end subroutine sw_poly_fit_2


! Updated 10/24/2001.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!   Program 4.4   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! Please Note:                                                          !
!                                                                       !
! (1) This computer program is written by Tao Pang in conjunction with  !
!     his book, "An Introduction to Computational Physics," published   !
!     by Cambridge University Press in 1997.                            !
!                                                                       !
! (2) No warranties, express or implied, are made for this program.     !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine sw_migs (A,N,X,INDX)
!
! subroutine to invert matrix A(N,N) with the inverse stored
! in X(N,N) in the output.  Copyright (c) Tao Pang 2001.
!
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I,J,K
  INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
  REAL (kind=real_kind), INTENT (INOUT), DIMENSION (N,N):: A
  REAL (kind=real_kind), INTENT (OUT), DIMENSION (N,N):: X
  REAL (kind=real_kind), DIMENSION (N,N) :: B
!
  DO I = 1, N
    DO J = 1, N
      B(I,J) = 0.0
    END DO
  END DO
  DO I = 1, N
    B(I,I) = 1.0
  END DO
!
  call sw_elgs (A,N,INDX)
!
  DO I = 1, N-1
    DO J = I+1, N
      DO K = 1, N
        B(INDX(J),K) = B(INDX(J),K)-A(INDX(J),I)*B(INDX(I),K)
      END DO
    END DO
  END DO
!
  DO I = 1, N
    X(N,I) = B(INDX(N),I)/A(INDX(N),N)
    DO J = N-1, 1, -1
      X(J,I) = B(INDX(J),I)
      DO K = J+1, N
        X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
      END DO
      X(J,I) =  X(J,I)/A(INDX(J),J)
    END DO
  END DO
end subroutine sw_migs


subroutine sw_elgs (A,N,INDX)
!
! subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
!
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I,J,K,ITMP
  INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
  REAL (kind=real_kind) :: C1,PI,PI1,PJ
  REAL (kind=real_kind), INTENT (INOUT), DIMENSION (N,N) :: A
  REAL (kind=real_kind), DIMENSION (N) :: C
!
! Initialize the index
!
  DO I = 1, N
    INDX(I) = I
  END DO
!
! Find the rescaling factors, one from each row
!
  DO I = 1, N
    C1= 0.0
    DO J = 1, N
      C1 = MAX(C1,ABS(A(I,J)))
    END DO
    C(I) = C1
  END DO
!
! Search the pivoting (largest) element from each column
!
  DO J = 1, N-1
    PI1 = 0.0
    DO I = J, N
      PI = ABS(A(INDX(I),J))/C(INDX(I))
      IF (PI.GT.PI1) THEN
        PI1 = PI
        K   = I
      ENDIF
    END DO
!
! Interchange the rows via INDX(N) to record pivoting order
!
    ITMP    = INDX(J)
    INDX(J) = INDX(K)
    INDX(K) = ITMP
    DO I = J+1, N
      PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
      A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
      DO K = J+1, N
        A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
      END DO
    END DO
  END DO
!
end subroutine sw_elgs

   subroutine sw_initialize_deformation_weights
                                      
!
! compute the cell coefficients for the deformation calculations
! WCS, 13 July 2010
!
      implicit none

!  local variables

      real (kind=real_kind), dimension(:,:), allocatable :: thetae
      real (kind=real_kind), dimension(:), allocatable :: theta_abs

      real (kind=real_kind), dimension(25) :: xc, yc, zc ! cell center coordinates
      real (kind=real_kind), dimension(25) :: thetav, thetat, dl_sphere
      real (kind=real_kind) :: dl, xec, yec, zec
      real (kind=real_kind) :: thetae_tmp, xe_tmp, ye_tmp
      integer :: i, j, k, ip1, ip2, m, n
      integer :: iCell, iEdge
      real (kind=real_kind), dimension(25) :: xp, yp
      
      integer :: cell_add
      integer, dimension(25) :: cell_list

      integer :: cell1, cell2, iv
      logical :: do_the_cell
      real (kind=real_kind) :: area_cell, sint2, cost2, sint_cost, area_cellt

      allocate(defc_a(maxEdges,nCells))
      allocate(defc_b(maxEdges,nCells))
      allocate(thetae(2, nEdges))
      allocate(theta_abs(nCells))

      defc_a(:,:) = 0.d0
      defc_b(:,:) = 0.d0

      do iCell = 1, nCells

         cell_list(1) = iCell
         do i = 2, nEdgesOnCell(iCell)+1
            cell_list(i) = CellsOnCell(i-1,iCell)
         end do
         n = nEdgesOnCell(iCell) + 1

!  check to see if we are reaching outside the halo

         do_the_cell = .true.
         do i = 1, n
            if (cell_list(i) > nCells) do_the_cell = .false.
         end do

         if (.not. do_the_cell) cycle


!  compute poynomial fit for this cell if all needed neighbors exist

         xc(1) = xCell(iCell) / sphere_radius
         yc(1) = yCell(iCell) / sphere_radius
         zc(1) = zCell(iCell) / sphere_radius


         do i = 2, n
            iv = verticesOnCell(i-1,iCell)
            xc(i) = xVertex(iv) / sphere_radius
            yc(i) = yVertex(iv) / sphere_radius
            zc(i) = zVertex(iv) / sphere_radius
         end do

         theta_abs(iCell) =  pi / 2. - sphere_angle( xc(1), yc(1), zc(1),  &
                                                    xc(2), yc(2), zc(2),  &
                                                    0.0_real_kind, 0.0_real_kind, 1.0_real_kind ) 

! angles from cell center to neighbor centers (thetav)

         do i = 1, n-1
   
            ip2 = i+2
            if (ip2 > n) ip2 = 2
    
            thetav(i) = sphere_angle( xc(1),   yc(1),   zc(1),    &
                                      xc(i+1), yc(i+1), zc(i+1),  &
                                      xc(ip2), yc(ip2), zc(ip2)   )

            dl_sphere(i) = sphere_radius * arc_length( xc(1  ), yc(1  ), zc(1  ),  &
                                                       xc(i+1), yc(i+1), zc(i+1) )
         end do

         thetat(1) = 0.  !  this defines the x direction, cell center 1 -> 
!         thetat(1) = theta_abs(iCell)  !  this defines the x direction, longitude line
         do i = 2, n-1
            thetat(i) = thetat(i-1) + thetav(i-1)
         end do
   
         do i = 1, n-1
            xp(i) = cos(thetat(i)) * dl_sphere(i)
            yp(i) = sin(thetat(i)) * dl_sphere(i)
         end do

!         thetat(1) = 0.
         thetat(1) = theta_abs(iCell)
         do i = 2, n-1
            ip1 = i+1
            if (ip1 == n) ip1 = 1
            thetat(i) = plane_angle( 0.0_real_kind, 0.0_real_kind, 0.0_real_kind,  &
                                     xp(i)-xp(i-1), yp(i)-yp(i-1), 0.0_real_kind,  &
                                     xp(ip1)-xp(i), yp(ip1)-yp(i), 0.0_real_kind,  &
                                     0.0_real_kind, 0.0_real_kind, 1.0_real_kind)
            thetat(i) = thetat(i) + thetat(i-1)
         end do

         area_cell = 0.
         area_cellt = 0.
         do i = 1, n-1
            ip1 = i+1
            if (ip1 == n) ip1 = 1
            dl = sqrt((xp(ip1)-xp(i))**2 + (yp(ip1)-yp(i))**2)
            area_cell  = area_cell + 0.25*(xp(i)+xp(ip1))*(yp(ip1)-yp(i)) - 0.25*(yp(i)+yp(ip1))*(xp(ip1)-xp(i))
            area_cellt = area_cellt + (0.25*(xp(i)+xp(ip1))*cos(thetat(i)) + 0.25*(yp(i)+yp(ip1))*sin(thetat(i)))*dl
         end do

         do i = 1, n-1
            ip1 = i+1
            if (ip1 == n) ip1 = 1
            dl = sqrt((xp(ip1)-xp(i))**2 + (yp(ip1)-yp(i))**2)
            sint2 = (sin(thetat(i)))**2
            cost2 = (cos(thetat(i)))**2
            sint_cost = sin(thetat(i))*cos(thetat(i))
            defc_a(i,iCell) = dl*(cost2 - sint2)/area_cell
            defc_b(i,iCell) = dl*2.*sint_cost/area_cell
            if (cellsOnEdge(1,edgesOnCell(i,iCell)) /= iCell) then
               defc_a(i,iCell) = - defc_a(i,iCell)
               defc_b(i,iCell) = - defc_b(i,iCell)
            end if
 
         end do

      end do

   end subroutine sw_initialize_deformation_weights

end module advection_mod
