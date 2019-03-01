module test_cases
    use module_para      , only: RKIND,a,g,pi,omega,config_test_case
    use module_array     , only: u,wu,wh,wh_s
    use diagnostics_tools, only: IAP
    use mesh_info        , only: nCells              ,&
                                 nEdges              ,&
                                 nVertices           ,&
                                 vertexDegree        ,&
                                 maxEdges            ,&
                                 maxEdges2           ,&
                                 latCell,lonCell     ,&
                                 latEdge,lonEdge     ,&
                                 latVertex,lonVertex ,&
                                 dvEdge,dcEdge       ,&
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
                                 wh_s
    
    use shallow_water_waves_test, only: getFields
    
    contains
    ! Choose test case
    subroutine choose_case
        if(config_test_case==2)then
            call sw_test_case_2
        elseif(config_test_case==5)then
            call sw_test_case_5
        elseif(config_test_case==6)then
            call sw_test_case_6
        elseif(config_test_case==7)then
            call sw_Shamir2016
        endif
        
        call IAP(u,wh,wu)
    end subroutine choose_case
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Setup shallow water test case 2: Global Steady State Nonlinear Zonal 
    !                                  Geostrophic Flow
    !
    ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
    !            Approximations to the Shallow Water Equations in Spherical 
    !            Geometry" J. of Comp. Phys., 102, pp. 211--224
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine sw_test_case_2
        implicit none

        real (kind=RKIND), parameter   :: u0    = 2.d0 * pi * a / (12.d0 * 86400.d0)
        !real (kind=RKIND), parameter   :: u0    = omega*a/12.d0
        real (kind=RKIND), parameter   :: gh0   = 29400.d0
        real (kind=RKIND), parameter   :: alpha = 0.d0
        
        real (kind=RKIND), allocatable :: psiVertex(:)

        integer                        :: iCell, iEdge, iVtx
        
        !
        ! Initialize wind field
        !
        allocate(psiVertex(nVertices))
        do iVtx = 1, nVertices
           psiVertex(iVtx) = -a * u0 * ( dsin(latVertex(iVtx)) * dcos(alpha)                         &
                                        -dcos(lonVertex(iVtx)) * dcos(latVertex(iVtx)) * dsin(alpha))
        end do
        
        do iEdge = 1,nEdges
           u(iEdge) = -( psiVertex(verticesOnEdge(2,iEdge)) &
                        -psiVertex(verticesOnEdge(1,iEdge)) ) / dvEdge(iEdge)
        end do
        deallocate(psiVertex)

        !
        ! Initialize height field (actually, fluid thickness field)
        !
        do iCell = 1, nCells
           wh(iCell) = gh0 - (a * omega * u0 + 0.5d0 * u0**2.d0) * &
                       (-dcos(lonCell(iCell)) * dcos(latCell(iCell)) * dsin(alpha) &
                        +dsin(latCell(iCell)) * dcos(alpha))**2.d0
        end do
        
        wh_s = 0.d0
        
    end subroutine sw_test_case_2


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Setup shallow water test case 5: Zonal Flow over an Isolated Mountain
    !
    ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
    !            Approximations to the Shallow Water Equations in Spherical 
    !            Geometry" J. of Comp. Phys., 102, pp. 211--224
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine sw_test_case_5
        implicit none
        
        real (kind=RKIND), parameter :: u0 = 20.d0
        real (kind=RKIND), parameter :: gh0 = 5960.d0*g
        real (kind=RKIND), parameter :: hs0 = 2000.d0
        real (kind=RKIND), parameter :: theta_c = pi/6.d0
        real (kind=RKIND), parameter :: lambda_c = 3.d0*pi/2.d0
        real (kind=RKIND), parameter :: rr = pi/9.d0
        real (kind=RKIND), parameter :: alpha = 0.d0
        
        integer :: iCell, iEdge, iVtx
        real (kind=RKIND) :: r
        real (kind=RKIND), allocatable, dimension(:) :: psiVertex
        !
        ! Initialize wind field
        !
        allocate(psiVertex(nVertices))
        do iVtx = 1, nVertices
           psiVertex(iVtx) = -a * u0 * ( sin(latVertex(iVtx)) * cos(alpha) - &
                                         cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) &
                                       )
        end do
        do iEdge = 1, nEdges
           u(iEdge) = -(  psiVertex(verticesOnEdge(2,iEdge)) &
                        - psiVertex(verticesOnEdge(1,iEdge)) &
                       ) / dvEdge(iEdge)
        end do
        deallocate(psiVertex)
        
        !
        ! Initialize mountain
        !
        do iCell = 1, nCells
           if (lonCell(iCell) < 0.0) lonCell(iCell) = lonCell(iCell) + 2.0 * pi
           r = sqrt(min(rr**2.0, (lonCell(iCell) - lambda_c)**2.0 + (latCell(iCell) - theta_c)**2.0))
           wh_s(iCell) = hs0 * (1.0 - r/rr)*g
        end do
        
        !
        ! Initialize height field (actually, fluid thickness field)
        !
        do iCell = 1, nCells
           wh(iCell) = gh0 - (a * omega * u0 + 0.5 * u0**2.0) * &
                             (-cos(lonCell(iCell)) * cos(latCell(iCell)) * sin(alpha) + &
                               sin(latCell(iCell)) * cos(alpha         )                &
                             )**2.0
           wh(iCell) = wh(iCell) - wh_s(iCell)
        end do

    end subroutine sw_test_case_5
    
    subroutine sw_test_case_6
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Setup shallow water test case 6: Rossby-Haurwitz Wave
    !
    ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
    !            Approximations to the Shallow Water Equations in Spherical 
    !            Geometry" J. of Comp. Phys., 102, pp. 211--224
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        implicit none
        
        real (kind=RKIND), parameter :: h0 = 8000.d0
        real (kind=RKIND), parameter :: w = 7.848e-6
        real (kind=RKIND), parameter :: K = 7.848e-6
        real (kind=RKIND), parameter :: R = 4.d0
        
        integer :: iCell, iEdge, iVtx
        real (kind=RKIND), dimension(nVertices) :: psiVertex
        real (kind=RKIND), dimension(nCells   ) :: aa,bb,cc
        !
        ! Initialize wind field
        !
        !allocate(psiVertex(nVertices))
        do iVtx = 1, nVertices
           psiVertex(iVtx) = -a * a * w * dsin(latVertex(iVtx)) + &
                              a *a * K * (dcos(latVertex(iVtx))**R) * &
                              dsin(latVertex(iVtx)) * dcos(R * lonVertex(iVtx))
        end do
        
        do iEdge = 1, nEdges
           u(iEdge) = -( psiVertex(verticesOnEdge(2,iEdge)) - &
                         psiVertex(verticesOnEdge(1,iEdge))) / dvEdge(iEdge)
        end do
        
        aa = 0.5d0 * w * (2.d0 * omega + w) * cos(latCell)**2.d0 + &
             0.25d0 * K**2.d0 * dcos(latCell)**(2.d0*R) * ((R+1.d0)*dcos(latCell)**2.d0 + 2.d0*R**2.d0 - R - 2.d0 - 2.d0*R**2.d0 * dcos(latCell)**(-2.d0))
        bb = (2.0*(omega + w)*K / ((R+1.0)*(R+2.0))) * dcos(latCell)**R * ((R**2.0 + 2.0*R + 2.0) - ((R+1.0)*dcos(latCell))**2.0)
        cc = 0.25 * K**2.0 * dcos(latCell)**(2.0*R) * ((R+1.0)*dcos(latCell)**2.0 - R - 2.0)
        
        
        !
        ! Initialize height field (actually, fluid thickness field)
        !
        do iCell = 1, nCells
           wh(iCell) = ( g * h0 + a*a*aa(iCell)                      &
                                + a*a*bb(iCell) * dcos(    R*lonCell(iCell))  &
                                + a*a*cc(iCell) * dcos(2.0*R*lonCell(iCell)))
        end do
        
        wh_s  = 0.d0
        
    end subroutine sw_test_case_6

    !========================================================================
    ! This test is according to :
    ! Shamir O, Paldor N. A quantitative test case for global-scale dynamical
    ! cores based on analytic wave solutions of the Shallow Water Equations.
    ! Submitted to: Q J ROY METEOR SOC, Feb 2016.
    !========================================================================
    subroutine sw_Shamir2016
        implicit none
        
        real (kind=RKIND) :: time(1)
        real (kind=RKIND) :: uu(nEdges,1),vv(nEdges,1),hh(nCells,1)
        real (kind=RKIND) :: theta(nEdges),beta(nEdges)
        
        time(1) = 0
        
        call getFields(latCell,lonCell,latEdge,lonEdge,time,uu,vv,hh,0)
        
        theta = datan(vv(:,1)/uu(:,1))
        where(vv(:,1)>0.and.uu(:,1)<0)theta = theta+pi
        where(vv(:,1)<0.and.uu(:,1)<0)theta = theta-pi
        beta  = theta-angleEdge
        u     = dcos(beta)*dsqrt(uu(:,1)*uu(:,1)+vv(:,1)*vv(:,1))
            
        wh    = hh(:,1)*g
        
        wh_s  = 0.d0
        
    end subroutine sw_Shamir2016
    
    
end module test_cases