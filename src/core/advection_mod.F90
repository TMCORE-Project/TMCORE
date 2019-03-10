! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module advection_mod

   use mesh_mod
   use const_mod, sphere_radius => radius

   contains

   subroutine advecion_init
                                      
!
! compute the cell coefficients for the polynomial fit.
! this is performed during setup for model integration.
! WCS, 31 August 2009
!
      implicit none

      integer              , dimension(:,:), allocatable :: advCells

!  local variables

      real (kind=real_kind), dimension(:,:), allocatable :: thetae
      real (kind=real_kind), dimension(:  ), allocatable :: xe, ye
      real (kind=real_kind), dimension(:  ), allocatable :: theta_abs

      real (kind=real_kind), dimension(25) :: xc, yc, zc ! cell center coordinates
      real (kind=real_kind), dimension(25) :: thetav, thetat, dl_sphere
      real (kind=real_kind) :: xm, ym, zm, dl, xec, yec, zec
      real (kind=real_kind) :: thetae_tmp, xe_tmp, ye_tmp
      real (kind=real_kind) :: xv1, xv2, yv1, yv2, zv1, zv2
      integer :: i, j, k, ip1, ip2, m, nReconstructCells, ip1a, ii
      integer :: iCell, iEdge
      real (kind=real_kind) :: x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5
      real (kind=real_kind) :: pdx1, pdx2, pdx3, pdy1, pdy2, pdy3, dx1, dx2, dy1, dy2
      real (kind=real_kind) :: angv1, angv2, dl1, dl2
      real (kind=real_kind), dimension(25) :: dxe, dye, x2v, y2v, xp, yp
      
      real (kind=real_kind) :: amatrix(25,25), bmatrix(25,25), wmatrix(25,25)
      real (kind=real_kind) :: length_scale
      integer :: ma,na, cell_add, mw, nn
      integer, dimension(25) :: cell_list

      integer :: cell1, cell2
      integer, parameter :: polynomial_order = 4
      logical, parameter :: least_squares = .true.
      logical :: add_the_cell, do_the_cell

      logical, parameter :: reset_poly = .true.

      real (kind=real_kind) :: rcell, cos2t, costsint, sin2t
      real (kind=real_kind), dimension(:), allocatable :: angle_2d

!----------------------------------------------------------------------------
      
      allocate(advCells (21,nCells))
      allocate(angle_2d (maxEdges))
      allocate(thetae   (2, nEdges))
      allocate(xe       (nEdges))
      allocate(ye       (nEdges))
      allocate(theta_abs(nCells))
      
      deriv_two(:,:,:) = 0.d0

      do iCell = 1, nCells !  is this correct? - we need first halo cell also...

        cell_list(1) = iCell
        do i = 2, nEdgesOnCell(iCell)+1
           cell_list(i) = cellsOnCell(i-1,iCell)
        end do
        nReconstructCells = nEdgesOnCell(iCell) + 1

        if ( polynomial_order > 2 .and. polynomial_order <= 4) then
           do i = 2, nEdgesOnCell(iCell) + 1
              do j = 1, nEdgesOnCell( cell_list(i) )
                 cell_add = CellsOnCell(j,cell_list(i))
                 add_the_cell = .true.
                 do k=1,nReconstructCells
                    if ( cell_add == cell_list(k) ) add_the_cell = .false.
                 end do
                 if (add_the_cell) then
                    nReconstructCells = nReconstructCells+1
                    cell_list(nReconstructCells) = cell_add
                 end if
              end do
           end do
        else
          stop 'Unknown polynomial order'
        end if

!  check to see if we are reaching outside the halo

        do_the_cell = .true.
        do i = 1, nReconstructCells
          if (cell_list(i) > nCells) do_the_cell = .false.
        end do

        if ( .not. do_the_cell ) cycle

!  compute polynomial fit for this cell if all needed neighbors exist

        do i = 1, nReconstructCells
          advCells(i,iCell) = cell_list(i)
          
          xc(i) = xCell(advCells(i,iCell)) / sphere_radius
          yc(i) = yCell(advCells(i,iCell)) / sphere_radius
          zc(i) = zCell(advCells(i,iCell)) / sphere_radius
        end do

        theta_abs(iCell) = 0.5d0*pi - sphere_angle( xc(1), yc(1), zc(1),&
                                                    xc(2), yc(2), zc(2),&
                                                    0.d0 , 0.d0 , 1.d0 ) 

! angles from cell center to neighbor centers (thetav)

        do i = 1, nReconstructCells-1
   
          ip2 = i+2
          if (ip2 > nReconstructCells) ip2 = 2
    
          thetav(i) = sphere_angle( xc(1  ), yc(1  ), zc(1  ), &
                                    xc(i+1), yc(i+1), zc(i+1), &
                                    xc(ip2), yc(ip2), zc(ip2) )

          dl_sphere(i) = sphere_radius*arc_length( xc(1  ), yc(1  ), zc(1  ), &
                                                   xc(i+1), yc(i+1), zc(i+1) )
        end do

        length_scale = 1.0_real_kind
        do i = 1, nReconstructCells-1
          dl_sphere(i) = dl_sphere(i) / length_scale
        end do

!        thetat(1) = 0.  !  this defines the x direction, cell center 1 -> 
        thetat(1) = theta_abs(iCell)  !  this defines the x direction, longitude line
        do i=2,nReconstructCells-1
          thetat(i) = thetat(i-1) + thetav(i-1)
        end do
   
        do i=1,nReconstructCells-1
          xp(i) = cos(thetat(i)) * dl_sphere(i)
          yp(i) = sin(thetat(i)) * dl_sphere(i)
        end do

        ma = nReconstructCells-1
        mw = nEdgesOnCell(iCell)

        bmatrix = 0.d0
        amatrix = 0.d0
        wmatrix = 0.d0

        if (polynomial_order == 2) then
          na = 6
          ma = ma+1
  
          amatrix(1,1) = 1.d0
          wmatrix(1,1) = 1.d0
          do i=2,ma
             amatrix(i,1) = 1.d0
             amatrix(i,2) = xp(i-1)
             amatrix(i,3) = yp(i-1)
             amatrix(i,4) = xp(i-1)**2
             amatrix(i,5) = xp(i-1) * yp(i-1)
             amatrix(i,6) = yp(i-1)**2
   
             wmatrix(i,i) = 1.d0
          end do
 
        else if (polynomial_order == 3) then
          na = 10
          ma = ma+1
  
          amatrix(1,1) = 1.d0
          wmatrix(1,1) = 1.d0
          do i=2,ma
             amatrix(i,1) = 1.d0
             amatrix(i,2) = xp(i-1)
             amatrix(i,3) = yp(i-1)
   
             amatrix(i,4) = xp(i-1)**2
             amatrix(i,5) = xp(i-1) * yp(i-1)
             amatrix(i,6) = yp(i-1)**2
   
             amatrix(i,7) = xp(i-1)**3
             amatrix(i,8) = yp(i-1) * (xp(i-1)**2)
             amatrix(i,9) = xp(i-1) * (yp(i-1)**2)
             amatrix(i,10) = yp(i-1)**3
   
             wmatrix(i,i) = 1.d0
 
          end do

        else if (polynomial_order == 4) then
          na = 15
          ma = ma+1
          
          amatrix(1,1) = 1.d0
          wmatrix(1,1) = 1.d0
          do i=2,ma
            amatrix(i,1) = 1.
            amatrix(i,2) = xp(i-1)
            amatrix(i,3) = yp(i-1)
          
            amatrix(i,4) = xp(i-1)**2
            amatrix(i,5) = xp(i-1) * yp(i-1)
            amatrix(i,6) = yp(i-1)**2
          
            amatrix(i,7) = xp(i-1)**3
            amatrix(i,8) = yp(i-1) * (xp(i-1)**2)
            amatrix(i,9) = xp(i-1) * (yp(i-1)**2)
            amatrix(i,10) = yp(i-1)**3
          
            amatrix(i,11) = xp(i-1)**4
            amatrix(i,12) = yp(i-1) * (xp(i-1)**3)
            amatrix(i,13) = (xp(i-1)**2)*(yp(i-1)**2)
            amatrix(i,14) = xp(i-1) * (yp(i-1)**3)
            amatrix(i,15) = yp(i-1)**4
          
            wmatrix(i,i) = 1.d0
          
          end do
          
          do i=1,mw
            wmatrix(i,i) = 1.d0
          end do
          
        else
          stop 'Unknown polynomial order'
        end if
 
        call swm_poly_fit_2( amatrix, bmatrix, wmatrix, ma, na, 25 )

        do i = 1, nEdgesOnCell(iCell)
          ip1 = i+1
          if (ip1 > nReconstructCells-1) ip1 = 1
  
          iEdge = edgesOnCell(i,iCell)
          xv1 = xVertex(verticesOnEdge(1, iEdge)) / sphere_radius
          yv1 = yVertex(verticesOnEdge(1, iEdge)) / sphere_radius
          zv1 = zVertex(verticesOnEdge(1, iEdge)) / sphere_radius
          xv2 = xVertex(verticesOnEdge(2, iEdge)) / sphere_radius
          yv2 = yVertex(verticesOnEdge(2, iEdge)) / sphere_radius
          zv2 = zVertex(verticesOnEdge(2, iEdge)) / sphere_radius
  
          call swm_arc_bisect( xv1, yv1, zv1,  &
                               xv2, yv2, zv2,  &
                               xec, yec, zec   )
  
          thetae_tmp = sphere_angle( xc(1),   yc(1),   zc(1),    &
                                     xc(i+1), yc(i+1), zc(i+1),  &
                                     xec,     yec,     zec       )
          thetae_tmp = thetae_tmp + thetat(i)
          if (iCell == cellsOnEdge(1,iEdge)) then
             thetae(1, iEdge) = thetae_tmp
          else
             thetae(2, iEdge) = thetae_tmp
          end if
  
        end do

!  fill second derivative stencil for rk advection 

        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i,iCell)
  
          if (iCell == cellsOnEdge(1,iEdge)) then
            cos2t    = cos(thetae(1, iEdge))
            sin2t    = sin(thetae(1, iEdge))
            costsint = cos2t*sin2t
            cos2t    = cos2t**2
            sin2t    = sin2t**2
   
            do j = 1, nReconstructCells
              deriv_two(j,1,iEdge) =   2.*cos2t    * bmatrix(4,j)  &
                                     + 2.*costsint * bmatrix(5,j)  &
                                     + 2.*sin2t    * bmatrix(6,j)
            end do
          else
            cos2t    = cos(thetae(2, iEdge))
            sin2t    = sin(thetae(2, iEdge))
            costsint = cos2t*sin2t
            cos2t    = cos2t**2
            sin2t    = sin2t**2
      
            do j=1,nReconstructCells
              deriv_two(j,2,iEdge) =   2.*cos2t    * bmatrix(4,j)  &
                                     + 2.*costsint * bmatrix(5,j)  &
                                     + 2.*sin2t    * bmatrix(6,j)
            end do
          end if
        end do
 
      end do ! end of loop over cells

   end subroutine advecion_init

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
   ! subroutine swm_arc_bisect
   !
   ! Returns the point C=(cx, cy, cz) that bisects the great circle arc from
   !   A=(ax, ay, az) to B=(bx, by, bz). It is assumed that A and B lie on the
   !   surface of a sphere centered at the origin.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine swm_arc_bisect(ax, ay, az, bx, by, bz, cx, cy, cz)
   
      implicit none
   
      real (kind=real_kind), intent(in) :: ax, ay, az, bx, by, bz
      real (kind=real_kind), intent(out) :: cx, cy, cz
   
      real (kind=real_kind) :: r           ! Radius of the sphere
      real (kind=real_kind) :: d           
   
      r = sqrt(ax*ax + ay*ay + az*az)
   
      cx = 0.5*(ax + bx)
      cy = 0.5*(ay + by)
      cz = 0.5*(az + bz)
   
      if (cx == 0. .and. cy == 0. .and. cz == 0.) then
         stop 'arc_bisect: A and B are diametrically opposite'
      else
         d = sqrt(cx*cx + cy*cy + cz*cz)
         cx = r * cx / d
         cy = r * cy / d
         cz = r * cz / d
      end if
   
   end subroutine swm_arc_bisect


   subroutine swm_poly_fit_2(a_in,b_out,weights_in,m,n,ne)

      implicit none

      integer, intent(in) :: m,n,ne
      real (kind=real_kind), dimension(ne,ne), intent(in) :: a_in, weights_in
      real (kind=real_kind), dimension(ne,ne), intent(out) :: b_out
   
      ! local storage
   
      real (kind=real_kind), dimension(m,n)  :: a
      real (kind=real_kind), dimension(n,m)  :: b
      real (kind=real_kind), dimension(m,m)  :: w,wt,h
      real (kind=real_kind), dimension(n,m)  :: at, ath
      real (kind=real_kind), dimension(n,n)  :: ata, ata_inv, atha, atha_inv
      integer, dimension(n) :: indx
      integer :: i,j
   
      if ( (ne<n) .or. (ne<m) ) then
         stop 'error in poly_fit_2 inversion'
      end if
   
!      a(1:m,1:n) = a_in(1:n,1:m) 
      a(1:m,1:n) = a_in(1:m,1:n)
      w(1:m,1:m) = weights_in(1:m,1:m) 
      b_out(:,:) = 0.   

      wt = transpose(w)
      h = matmul(wt,w)
      at = transpose(a)
      ath = matmul(at,h)
      atha = matmul(ath,a)
      
      ata = matmul(at,a)

!      if (m == n) then
!         call swm_migs(a,n,b,indx)
!      else

         call swm_migs(atha,n,atha_inv,indx)

         b = matmul(atha_inv,ath)

!         call swm_migs(ata,n,ata_inv,indx)
!         b = matmul(ata_inv,at)
!      end if
      b_out(1:n,1:m) = b(1:n,1:m)

!     do i=1,n
!        write(6,*) ' i, indx ',i,indx(i)
!     end do
!
!     write(6,*) ' '

   end subroutine swm_poly_fit_2


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
subroutine swm_migs (A,N,X,INDX)
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
  call swm_elgs (A,N,INDX)
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
end subroutine swm_migs


subroutine swm_elgs (A,N,INDX)
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
end subroutine swm_elgs

end module advection_mod
