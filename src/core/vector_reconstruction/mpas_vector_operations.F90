! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_vector_operations
!
!> \brief MPAS vector operations
!> \author Mark Petersen
!> \date   April 2013
!> \details
!>  This module contains the routines involving vector operations
!
!-----------------------------------------------------------------------
module mpas_vector_operations

   use const_mod
   use mesh_mod
   use log_mod

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: mpas_vec_mag_in_r3, &
             mpas_initialize_vectors, &
             mpas_initialize_tangent_vectors, &
             mpas_unit_vec_in_r3, &
             mpas_cross_product_in_r3, &
             mpas_tangential_velocity, &
             mpas_vector_R3Cell_to_2DEdge, &
             mpas_vector_R3Cell_to_normalVectorEdge, &
             mpas_vector_R3_to_LonLatR, &
             mpas_vector_LonLatR_to_R3, &
             mpas_zonal_meridional_vectors, &
             mpas_fix_periodicity, &
             mpas_unit_test_fix_periodicity


   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains


!***********************************************************************
!
!  routine mpas_vec_mag_in_r3
!
!> \brief   MPAS 3D vector magnitude routine
!> \author  Matt Hoffman
!> \date    13 Jan. 2015
!> \details
!> This routine calculates the magnitude of a 3d vector.
!-----------------------------------------------------------------------
  real(real_kind) function mpas_vec_mag_in_r3(xin)!{{{
    implicit none
    real(real_kind), dimension(3), intent(in) :: xin !< Input: Vector
    mpas_vec_mag_in_r3 = sqrt(xin(1)**2 + xin(2)**2 + xin(3)**2)
  end function mpas_vec_mag_in_r3!}}}

!***********************************************************************
!
!  routine mpas_unit_vec_in_r3
!
!> \brief   MPAS 3D unit vector routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details
!> This routine creates a unit vector out of an input point.
!-----------------------------------------------------------------------
  subroutine mpas_unit_vec_in_r3(xin)!{{{
    implicit none
    real(real_kind), dimension(3), intent(inout) :: xin !< Input/Output: Vector and unit vector
    real(real_kind) :: mag
    mag = sqrt(xin(1)**2+xin(2)**2+xin(3)**2)
    xin(:) = xin(:) / mag
  end subroutine mpas_unit_vec_in_r3!}}}

!***********************************************************************
!
!  routine mpas_cross_product_in_r3
!
!> \brief   MPAS 3D cross product routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details
!> This routine computes the cross product of two input vectors.
!-----------------------------------------------------------------------
  subroutine mpas_cross_product_in_r3(p_1,p_2,p_out)!{{{
    real(real_kind), intent(in ) :: p_1 (3) !< Input: Vector 1
    real(real_kind), intent(in ) :: p_2 (3) !< Input: Vector 2
    real(real_kind), intent(out) :: p_out (3) !< Output: Cross product of vector 1 and vector 2

    p_out(1) = p_1(2)*p_2(3)-p_1(3)*p_2(2)
    p_out(2) = p_1(3)*p_2(1)-p_1(1)*p_2(3)
    p_out(3) = p_1(1)*p_2(2)-p_1(2)*p_2(1)
  end subroutine mpas_cross_product_in_r3!}}}

!***********************************************************************
!
!  routine mpas_vector_R3Cell_to_2DEdge
!
!> \brief   Convert an R3 cell-centered vector field to 2D vectors at edges
!> \author  Mark Petersen
!> \date    April 2013
!> \details
!>  Convert an R3 cell-centered vector field to 2D vectors at edges, where
!>  the local coordinate system at the edge is normal and tangent to that edge.
!
!-----------------------------------------------------------------------

   subroutine mpas_vector_R3Cell_to_2DEdge(vectorR3Cell, edgeTangentVectors, normalVectorEdge, tangentialVectorEdge)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real(real_kind), dimension(:,:), intent(in) :: vectorR3Cell       !< Input: 3-vector located at cell centers, in x,y,z coordinates
      real(real_kind), dimension(:,:), intent(in) :: edgeTangentVectors !< Input: unit vector tangent to an edge

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real(real_kind), dimension(:), intent(out) :: normalVectorEdge     !< Output: normal component of vector at edge
      real(real_kind), dimension(:), intent(out) :: tangentialVectorEdge !< Output: tangential component of vector at edge

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, cell1, cell2

      real(real_kind), dimension(3) :: vectorR3Edge

      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         ! average neighboring cell-centered vectors to the edge
         vectorR3Edge(:) = 0.5d0*(vectorR3Cell(:,cell1) + vectorR3Cell(:,cell2))

         ! normal component at edge: take dot products with unit vectors at edge
         normalVectorEdge    (iEdge) = sum(vectorR3Edge(:)*edgeNormalVectors (:,iEdge))
         tangentialVectorEdge(iEdge) = sum(vectorR3Edge(:)*edgeTangentVectors(:,iEdge))
      enddo

   end subroutine mpas_vector_R3Cell_to_2DEdge!}}}

!***********************************************************************
!
!  routine mpas_vector_R3Cell_to_normalVectorEdge
!
!> \brief   Interpolate from R3 cell-centered vector field vector normal to edge
!> \author  Mark Petersen
!> \date    April 2013
!> \details
!>  Convert an R3 cell-centered vector field to normal vectors at edges.
!
!-----------------------------------------------------------------------

   subroutine mpas_vector_R3Cell_to_normalVectorEdge(vectorR3Cell, normalVectorEdge)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real(real_kind), dimension(:,:), intent(in) :: vectorR3Cell  !< Input: 3-vector located at cell centers, in x,y,z coordinates

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real(real_kind), dimension(:,:), intent(out) :: normalVectorEdge   !< Output: normal component of vector at edge

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, cell1, cell2

      real(real_kind), dimension(3) :: vectorR3Edge
      real(real_kind), dimension(:,:), pointer :: edgeNormalVectors

      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         
         ! average neighboring cell-centered vectors to the edge
         vectorR3Edge(:) = 0.5d0*(vectorR3Cell(:,cell1) + vectorR3Cell(:,cell2))

         ! normal component at edge: take dot products with unit vectors at edge
         normalVectorEdge(:,iEdge) = sum(vectorR3Edge(:)*edgeNormalVectors(:,iEdge))
      enddo

   end subroutine mpas_vector_R3Cell_to_normalVectorEdge!}}}

!***********************************************************************
!
!  routine mpas_tangential_velocity
!
!> \brief   compute tangential velocity at an edge from the normal velocities
!> \author  Mark Petersen
!> \date    April 2013
!> \details
!>  Compute tangential velocity at an edge from the normal velocities on the
!>  edges of the two neighboring cells.
!
!-----------------------------------------------------------------------

   subroutine mpas_tangential_velocity(normalVelocity, tangentialVelocity)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real(real_kind), dimension(:), intent(in ) :: normalVelocity   !< Input: Horizontal velocity normal to edge

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real(real_kind), dimension(:), intent(out) :: tangentialVelocity   !< Output: Horizontal velocity tangent to edge

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, i, k, eoe

      do iEdge = 1,nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            
            tangentialVelocity(iEdge) = tangentialVelocity(iEdge) + weightsOnEdge(i,iEdge) * normalVelocity(eoe)
         end do
      end do

   end subroutine mpas_tangential_velocity!}}}


!***********************************************************************
!
!  routine mpas_zonal_meridional_vectors
!
!> \brief   Computes zonal, meridional, and vertical unit vectors
!> \author  Mark Petersen
!> \date    1 May 2013
!> \details
!>  Given a latitude and longitude location, compute unit vectors pointing
!>  in the zonal, meridional, and vertical directions.
!
!-----------------------------------------------------------------------

   subroutine mpas_zonal_meridional_vectors(lon, lat, zonalUnitVector, meridionalUnitVector, verticalUnitVector)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real(real_kind), intent(in) :: &
         lon, &!< Input: longitude, in radians, ranging [0,2*pi]
         lat   !< Input: latitude,  in radians, ranging [-pi,pi]

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real(real_kind), dimension(3), intent(out) :: &
         zonalUnitVector,     &!< Output: zonal unit vector
         meridionalUnitVector,&!< Output: meridional unit vector
         verticalUnitVector    !< Output: vertical unit vector

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real(real_kind) :: sin_lat, sin_lon, cos_lat, cos_lon


      sin_lat = sin(lat)
      cos_lat = cos(lat)
      sin_lon = sin(lon)
      cos_lon = cos(lon)

      zonalUnitVector(1) = - sin_lon
      zonalUnitVector(2) =   cos_lon
      zonalUnitVector(3) =   0

      meridionalUnitVector(1) = - sin_lat * cos_lon
      meridionalUnitVector(2) = - sin_lat * sin_lon
      meridionalUnitVector(3) =   cos_lat

      verticalUnitVector(1) = cos_lat * cos_lon
      verticalUnitVector(2) = cos_lat * sin_lon
      verticalUnitVector(3) = sin_lat

   end subroutine mpas_zonal_meridional_vectors!}}}


!***********************************************************************
!
!  routine mpas_vector_R3_to_LonLatR
!
!> \brief   Convert an R3 vector to a vector in spherical coordinates
!> \author  Mark Petersen
!> \date    April 2013
!> \details
!>  Given an R3 vector, rotate it to spherical coordinates
!
!-----------------------------------------------------------------------

   subroutine mpas_vector_R3_to_LonLatR(vectorR3, lon, lat, vectorLonLatR)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real(real_kind), dimension(3), intent(in) :: &
         vectorR3  !< Input: vector in R3

      real(real_kind), intent(in) :: &
         lon, &!< Input: longitude, in radians, ranging [0,2*pi]
         lat   !< Input: latitude,  in radians, ranging [-pi,pi]

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real(real_kind), dimension(3), intent(out) :: &
         vectorLonLatR   !< Output: vector in spherical coordinates

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: i,j

      real(real_kind), dimension(3  ) :: zonalUnitVector, meridionalUnitVector, verticalUnitVector
      real(real_kind), dimension(3,3) :: rotationMatrix

      call mpas_zonal_meridional_vectors(lon, lat, zonalUnitVector, meridionalUnitVector, verticalUnitVector)

      rotationMatrix(:,1) = zonalUnitVector
      rotationMatrix(:,2) = meridionalUnitVector
      rotationMatrix(:,3) = verticalUnitVector

      vectorLonLatR = 0.0d0
      do i=1,3
         do j=1,3
            ! xi = R^T x
            ! where ^T is the transpose.  Note indices on R are (j,i)
            vectorLonLatR(i) = vectorLonLatR(i) + rotationMatrix(j,i)*vectorR3(j)
         enddo
      enddo

   end subroutine mpas_vector_R3_to_LonLatR!}}}


!***********************************************************************
!
!  routine mpas_vector_LonLatR_to_R3
!
!> \brief   Convert a vector in spherical coordinates to an R3 vector
!> \author  Mark Petersen
!> \date    April 2013
!> \details
!>  Given a vector in spherical coordinates, rotate it R3
!
!-----------------------------------------------------------------------

   subroutine mpas_vector_LonLatR_to_R3(vectorLonLatR, lon, lat, vectorR3)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real(real_kind), dimension(3), intent(in) :: &
         vectorLonLatR   !< Input: vector in spherical coordinates

      real(real_kind), intent(in) :: &
         lon, &!< Input: longitude, in radians, ranging [0,2*pi]
         lat   !< Input: latitude,  in radians, ranging [-pi,pi]

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real(real_kind), dimension(3), intent(out) :: &
         vectorR3  !< Output: vector in R3

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: i,j

      real(real_kind), dimension(3) :: zonalUnitVector, meridionalUnitVector, verticalUnitVector
      real(real_kind), dimension(3,3) :: rotationMatrix

      call mpas_zonal_meridional_vectors(lon, lat, zonalUnitVector, meridionalUnitVector, verticalUnitVector)

      rotationMatrix(:,1) = zonalUnitVector
      rotationMatrix(:,2) = meridionalUnitVector
      rotationMatrix(:,3) = verticalUnitVector

      vectorR3 = 0.0d0
      do i=1,3
         do j=1,3
            vectorR3(i) = vectorR3(i) + rotationMatrix(i,j)*vectorLonLatR(j)
         enddo
      enddo

   end subroutine mpas_vector_LonLatR_to_R3!}}}

!***********************************************************************
!
!  routine mpas_initialize_vectors
!
!> \brief   MPAS RBF interpolation initialization routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details
!> This routine computes geometric fields that will be potentially useful for calling
!> the interpolation routines.
!> Input: the mesh
!> Output:
!>      edgeNormalVectors - the unit vector normal to the edge and tangent to the sphere
!>      cellTangentPlane - 2 orthogonal unit vectors in the tangent plane of each cell
!>                         The first unit vector is chosen to point toward the center of the first
!>                         edge on the cell.
!>      localVerticalUnitVectors - the unit normal vector of the tangent plane at the center
!>                         of each cell
!-----------------------------------------------------------------------
  subroutine mpas_initialize_vectors!{{{

    implicit none

    integer :: iEdge, iCell, cell1, cell2
    real(real_kind), dimension(3) :: xHatPlane, yHatPlane, rHat
    real(real_kind)               :: normalDotRHat

    ! init arrays
    edgeNormalVectors = 0
    localVerticalUnitVectors = 0

    ! loop over all cells to be solved on this block
    do iCell = 1, nCells
      localVerticalUnitVectors(1,iCell) = xCell(iCell)
      localVerticalUnitVectors(2,iCell) = yCell(iCell)
      localVerticalUnitVectors(3,iCell) = zCell(iCell)
      
      call mpas_unit_vec_in_r3(localVerticalUnitVectors(:,iCell))
    end do

    ! Initialize normal unit vectors at each edge
    ! These vectors point from cell to cell.
    ! At boundaries, one cell does not exist, so it points from cell to edge or from edge to cell.
    do iEdge = 1,nEdges
      cell1 = cellsOnEdge(1,iEdge)
      cell2 = cellsOnEdge(2,iEdge)

      if (cell1 == nCells+1) then ! this is a boundary edge
        ! the normal points from the edge location to the cell location
        edgeNormalVectors(1,iEdge) = xCell(cell2) - xEdge(iEdge)
        edgeNormalVectors(2,iEdge) = yCell(cell2) - yEdge(iEdge)
        edgeNormalVectors(3,iEdge) = zCell(cell2) - zEdge(iEdge)
      else if (cell2 == nCells+1) then ! this is a boundary edge
        ! the normal points from the cell location to the edge location
        edgeNormalVectors(1,iEdge) = xEdge(iEdge) - xCell(cell1)
        edgeNormalVectors(2,iEdge) = yEdge(iEdge) - yCell(cell1)
        edgeNormalVectors(3,iEdge) = zEdge(iEdge) - zCell(cell1)
      else ! this is not a boundary cell
        ! the normal points from the cell 1 to cell2
        ! mrp problem: on periodic domains, vectors on edges of domain point the wrong way.
        edgeNormalVectors(1,iEdge) = xCell(cell2) - xCell(cell1)
        edgeNormalVectors(2,iEdge) = yCell(cell2) - yCell(cell1)
        edgeNormalVectors(3,iEdge) = zCell(cell2) - zCell(cell1)
      endif
      
      call mpas_unit_vec_in_r3(edgeNormalVectors(:,iEdge))
    end do

    do iCell=1,nCells
      iEdge = edgesOnCell(1,iCell)
      ! xHat and yHat are a local basis in the plane of the horizontal cell
      ! we arbitrarily choose xHat to point toward the first edge
      rHat          = localVerticalUnitVectors(:,iCell)
      normalDotRHat = sum(edgeNormalVectors(:,iEdge) * rHat)
      xHatPlane     = edgeNormalVectors(:,iEdge) - normalDotRHat * rHat
      call mpas_unit_vec_in_r3(xHatPlane)

      call mpas_cross_product_in_r3(rHat, xHatPlane, yHatPlane)
      call mpas_unit_vec_in_r3(yHatPlane) ! just to be sure...
      cellTangentPlane(:,1,iCell) = xHatPlane
      cellTangentPlane(:,2,iCell) = yHatPlane
    end do

  end subroutine mpas_initialize_vectors!}}}

!***********************************************************************
!
!  routine mpas_initialize_tangent_vectors
!
!> \brief   Initialize edgeTangentVectors
!> \author  Mark Petersen
!> \date    June 2013
!> \details
!> This routine initializes the array edgeTangentVectors.  It is not
!> included in the mpas_initialize_vectors subroutine because the
!> array edgeTangentVectors is not included in all cores.
!> Input: the mesh
!> Output:
!>      edgeTangentVectors - the unit vector tangent to the edge and tangent to the sphere
!-----------------------------------------------------------------------
   subroutine mpas_initialize_tangent_vectors(edgeTangentVectors)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real(real_kind), dimension(:,:), intent(out) :: edgeTangentVectors   !< Output: unit vector tangent to an edge

      integer :: iEdge, vertex1, vertex2

      do iEdge = 1,nEdges
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)
         
         edgeTangentVectors(1,iEdge) = xVertex(vertex2) - xVertex(vertex1)
         edgeTangentVectors(2,iEdge) = yVertex(vertex2) - yVertex(vertex1)
         edgeTangentVectors(3,iEdge) = zVertex(vertex2) - zVertex(vertex1)
           
         call mpas_unit_vec_in_r3(edgeTangentVectors(:,iEdge))
      end do

  end subroutine mpas_initialize_tangent_vectors!}}}

!***********************************************************************
!
!  subroutine mpas_fix_periodicity
!
!> \brief   Fixes periodicity of point pxi relative to a point xci and xiRef periodicity.
!> \author  Phillip Wolfram & Doug Jacobsen
!> \date    06/29/2015
!> \details
!>  This routine recomputes the location of a point pxi relative to a
!>  point xci with xiRef periodicity. The calculation ensures that the recomputed
!>  point and xci are spatially nearby, permitting reasonable calculation of
!>  relative geometry.  This function operates only on a single dimension "i" for
!>  a point pxi relative to a specific location xci for dimension xiRef
!>  in the periodic direction.  Note, pxi can only be adjusted by at most
!>  a single period, e.g., pxi = 260 + 360 = 620, xci = 0, and xiRef = 360
!>  returns 260, not -100.
!-----------------------------------------------------------------------
    real(real_kind) function mpas_fix_periodicity(pxi, xci, xiRef) !{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(real_kind), intent(in) :: pxi, xci, xiRef
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      ! real(real_kind), intent(out) :: mpas_fix_periodicity

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real(real_kind) :: dist

      dist = pxi - xci

      if (abs(dist) > xiRef * 0.5d0) then
        mpas_fix_periodicity = pxi - (dist/abs(dist)) * xiRef
      else
        mpas_fix_periodicity = pxi
      end if

    end function mpas_fix_periodicity !}}}

!***********************************************************************
!
!  subroutine mpas_unit_test_fix_periodicity
!
!> \brief   Simple unit test for testing the mpas_fix_periodicity routine.
!> \author  Phillip Wolfram
!> \date    06/29/2015
!> \details
!>  This routine tests the mpas_fix_periodicity routine.
!-----------------------------------------------------------------------
   subroutine mpas_unit_test_fix_periodicity(ierr)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: ierr
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real(real_kind) :: x1, x2, xc, xLen, xnew
      real(real_kind) :: eps = 1.0e-12

      ierr = 0

      xLen = 2d0*pi
      x1 = pi
      x2 = 3.d0*pi
      xc = pi

      xnew = mpas_fix_periodicity(x1, xc, xLen)
      if (abs(xnew - pi ) > eps) then
         call log_error("Error in mpas_unit_test_fix_periodicity: " // "x1's periodicity fix has error greater than tolerance -- FAILED.", __FILE__, __LINE__)
         ierr = 1
      else
         call log_notice('mpas_unit_test_fix_periodicity: x1 test - SUCCESS')
      endif

      xnew = mpas_fix_periodicity(x2, xc, xLen)
      if (abs(xnew - pi ) > eps) then
         call log_error("Error in mpas_unit_test_fix_periodicity: " // "x2's periodicity fix has error greater than tolerance -- FAILED.", __FILE__, __LINE__)
         ierr = 1
      else
         call log_notice('mpas_unit_test_fix_periodicity: x2 test - SUCCESS')
      endif

   end subroutine mpas_unit_test_fix_periodicity!}}}

end module mpas_vector_operations


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker

