SUBROUTINE B_operator(wu,wv,wh,du,dv,dh,bu,bv,bh,dt)
!
! This subroutine is to calculate the tendency of geopotential height
!
    use module_para,only : nx,nx1,ny,ny1
!
 	implicit none
!
	real*8,dimension(0:nx1,0:ny1),intent(in )  :: wu               ! wu = U = h*u
	real*8,dimension(0:nx1,0:ny1),intent(in )  :: wv               ! wv = V = h*v 
	real*8,dimension(0:nx1,0:ny1),intent(in )  :: wh               ! Geopotential height
	real*8,dimension(0:nx1,0:ny1),intent(in )  :: du               ! Tendency of wu
	real*8,dimension(0:nx1,0:ny1),intent(in )  :: dv               ! Tendency of wv
	real*8,dimension(0:nx1,0:ny1),intent(in )  :: dh               ! Tendency of wh
	real*8,dimension(0:nx1,0:ny1),intent(out)  :: bu               ! Dissipation of wu
	real*8,dimension(0:nx1,0:ny1),intent(out)  :: bv               ! Dissipation of wv
	real*8,dimension(0:nx1,0:ny1),intent(out)  :: bh               ! Dissipation od wh
	real*8,dimension(0:nx1,0:ny1)              :: su               ! working array
	real*8,dimension(0:nx1,0:ny1)              :: sv               ! working array
	real*8,dimension(0:nx1,0:ny1)              :: sh               ! working array
	real*8,dimension(0:nx1,0:ny1)              :: tu               ! working array
	real*8,dimension(0:nx1,0:ny1)              :: tv               ! working array
	real*8,dimension(0:nx1,0:ny1)              :: th               ! working array
	real*8                                     :: dt               ! time stepsize
	real*8                                     :: dt2              ! working variable
	integer                                    :: i,j              ! working variables
    real*8,external                            :: inner
    
    ! Antisymmetry Check
	!print*,inner(wu,wv,wh,du,dv,dh)
!
!   To calculate the consisitent dissipation operator using the 4th-order Runge-Kutta scheme
!
    dt2 = dt*0.5d0
    
    do j=1,ny
        do i=1,nx
            su(i,j) = wu(i,j)-dt2*du(i,j)
            sv(i,j) = wv(i,j)-dt2*dv(i,j)
            sh(i,j) = wh(i,j)-dt2*dh(i,j)
        enddo
    enddo
    
	call L_operator(su,sv,sh,tu,tv,th)
    !print*,inner(su,sv,sh,tu,tv,th)
    
    do j=1,ny
        do i=1,nx
             bu(i,j) = du(i,j)+2.0d0*tu(i,j)
             bv(i,j) = dv(i,j)+2.0d0*tv(i,j)
             bh(i,j) = dh(i,j)+2.0d0*th(i,j)
        
             su(i,j) = wu(i,j)-dt2*tu(i,j)
             sv(i,j) = wv(i,j)-dt2*tv(i,j)
             sh(i,j) = wh(i,j)-dt2*th(i,j)
        enddo
    enddo
    
	call L_operator(su,sv,sh,tu,tv,th)
    !print*,inner(su,sv,sh,tu,tv,th)
    
    do j=1,ny
        do i=1,nx
            bu(i,j) = bu(i,j)+2.0d0*tu(i,j)
            bv(i,j) = bv(i,j)+2.0d0*tv(i,j)
            bh(i,j) = bh(i,j)+2.0d0*th(i,j)
            
            su(i,j) = wu(i,j)-dt*tu(i,j)
            sv(i,j) = wv(i,j)-dt*tv(i,j)
            sh(i,j) = wh(i,j)-dt*th(i,j)
        enddo
    enddo
    
	call L_operator(su,sv,sh,tu,tv,th)
    !print*,inner(su,sv,sh,tu,tv,th)
    
    do j=1,ny
        do i=1,nx
            bu(i,j) = (du(i,j) - (bu(i,j)+tu(i,j))/6.0d0)/dt2 
            bv(i,j) = (dv(i,j) - (bv(i,j)+tv(i,j))/6.0d0)/dt2 
            bh(i,j) = (dh(i,j) - (bh(i,j)+th(i,j))/6.0d0)/dt2 
        enddo
    enddo
    
END SUBROUTINE B_operator
    
