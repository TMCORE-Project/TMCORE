MODULE module_array
    use mesh_info , only: nCells,nEdges
    
    implicit none

    real*8,allocatable,dimension(:) :: u               ! zonal wind
    real*8,allocatable,dimension(:) :: v               ! meridional wind
    real*8,allocatable,dimension(:) :: wh              ! geopotential height
    real*8,allocatable,dimension(:) :: wh_s            ! surface height
                                  
    real*8,allocatable,dimension(:) :: wu              ! u*sqrt(wh)
    real*8,allocatable,dimension(:) :: wv              ! v*sqrt(wh)
    real*8,allocatable,dimension(:) :: h               ! sqrt(wh)
    
    contains
    
    subroutine init_arrays
    implicit none
    allocate(             &
             u   (nEdges),&
             v   (nEdges),&
             wh  (nCells),&
             wh_s(nCells),&
             wu  (nEdges),&
             wv  (nEdges),&
             h   (nCells) &
            )
    
    wh_s = 0.d0
    
    end subroutine init_arrays
    
END MODULE module_array