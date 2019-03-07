!========================================================================
! This code is part of the online supporting material for:
! Shamir O, Paldor N. A quantitative test case for global-scale dynamical
! cores based on analytic wave solutions of the Shallow Water Equations.
! Submitted to: Q J ROY METEOR SOC, Feb 2016.
!========================================================================

!========================================================================
! The following subroutines calculate the analytic fields for the
! proposed test case on arbitrary lat x lon grids (See common arguments
! list below).
!
! getPhaseSpeed(C,waveFlag):
!   Returns the analytic phase speed for the desired wave.
!
! getFields(lat,lon,lev,time,u,v,h,waveFlag):
!   Returns the analytic velocity and free-surface height anomaly fields.
!
! getSurfacePressure(h,Ps):
!   Transforms the free-surface height field into surface pressure.
!
! getInitialTemperature(T):
!   Returns the initial temperature field used in the text.
!
! Arguments:
!     lat       - 1D array of desired latitudes (radians) for output
!     lon       - 1D array of desired longitudes (radians) for output
!     lev       - 1D array of desired vertical levels for output
!     time      - 1D array of desired times (sec) for output
!     waveFlag  - integer flag for desired wave: WIG=-1, Rossby=0, EIG=1
!     C         - phase speed (rad/sec)
!     u         - 4D array (lat,lon,lev,time) for output analytic zonal velocity field (m/sec)
!     v         - 4D array (lat,lon,lev,time) for output analytic meridional velocity field (m/sec)
!     h         - 3D array (lat,lon,time) for output analytic free-surface height anomaly field (m)
!     Ps        - 3D array (lat,lon,time) for output surface pressure (Pa)
!     T         - 4D array (lat,lon,lev,time) for output initial temperature field field (Kelvin)
!
! Notes:
!   1)  Since the analytic velocity fields of the proposed test case are
!       vertically homogeneous, the choice of vertical levels in
!       getFields is arbitrary. Only the length of the 1D lev array is
!       actually used.
!
!   2)  In accordance with the proposed test's parameters, this code is
!       only valid for gH=5e4 m^2/s^2 and (n,k)=(5,10). Attempting to
!       use this code as is with different values would yield erroneous
!       analytic solutions.
!========================================================================

module shallow_water_waves_test

  use const_mod, a => radius

  implicit none

  integer, parameter :: dp = real_kind

!========================================================================
! parameters
!========================================================================
  real(kind=dp), parameter :: H0    = 5.0e3_dp                 ! Layer's mean depth (m)
  real(kind=dp), parameter :: P0    = 1.0e5_dp                 ! Reference pressure (Pa)
  real(kind=dp), parameter :: Rd    = 287.0_dp                 ! Gas const for dry air (J/kg/K)

  integer,       parameter :: n     = 5                        ! chosen mode number
  integer,       parameter :: k     = 10                       ! chosen wave number
  real(kind=dp), parameter :: &
       sigma = 0.5_dp + ( 0.25_dp + k**2 ) ** 0.5_dp           ! see Eq.10 in text

  private :: H0, P0, Rd, n, k, sigma

  private :: getPsi,                &
             getAmplitudes

  public  :: getPhaseSpeed,         &
             getFields,             &
             getSurfacePressure,    &
             getInitialTemperature

contains

!========================================================================
! phase speed
!========================================================================
  subroutine getPhaseSpeed(C,waveFlag)

    implicit none

    integer,          intent(in)    :: waveFlag    ! -1=WIG, 0=Rossby, 1=EIG
    real(kind=dp),    intent(inout) :: C           ! phase speed (rad/sec)

    real(kind=dp)                   :: Cj(1:3)     ! see Eq. 7  in text
    real(kind=dp)                   :: Delta0      ! see Eq. 8a in text
    complex(kind=dp)                :: Deltaj      ! see eq. 8b in text
    real(kind=dp)                   :: Delta4      ! see Eq. 8c in text
    real(kind=dp)                   :: En          ! see Eq. 9  in text

    complex(kind=dp), parameter     :: i1 = cmplx(0.0_dp,1.0_dp,dp)
    complex(kind=dp), parameter     :: r1 = cmplx(1.0_dp,0.0_dp,dp)
    real(kind=dp),    parameter     :: r2 = 0.5_dp
    real(kind=dp),    parameter     :: r3 = 1.0_dp/3.0_dp

    integer                         :: j

    En     =  g*H0 / a**2 * ( n + sigma )**2
    Delta0 =  03.0_dp * k**2 * En
    Delta4 = -54.0_dp * k**4 * g*H0 * omega / a**2

    do j=1,3
       Deltaj = ( r1 * ( Delta4**2 - 4.0_dp * Delta0**3 ) )**r2
       Deltaj = ( r2 * ( Delta4    +          Deltaj    ) )**r3
       Deltaj = Deltaj * exp(2.0_dp*pi * i1 * j * r3)

       Cj(j)  = real( -r3 / k**2 * ( Deltaj + Delta0 / Deltaj ) )
    end do

    select case (waveFlag)
      case(0)
        C = -minval(abs(Cj))
      case(1)
        C = maxval(Cj)
      case(-1)
        C = minval(Cj)
    end select

  end subroutine getPhaseSpeed

!========================================================================
! eigenfunction
!========================================================================
  subroutine getPsi(lat,psi,dpsi)

    implicit none

    real(kind=dp), intent(in)    :: lat(:)            ! latitude (rad)
    real(kind=dp), intent(inout) :: psi(:)            ! see Eq. 14 in text
    real(kind=dp), intent(inout) :: dpsi(:)           ! meridional derivative of psi

    real(kind=dp), parameter     :: amp = 1.0e-8_dp   ! arbitrary amplitude for linear waves
    real(kind=dp)                :: C5(size(lat))     ! see Eq. 19 in text
    real(kind=dp)                :: C5p(size(lat))    ! meridional derivative of C5

    real(kind=dp), parameter     :: a3 = sigma * ( sigma + 1 ) * ( sigma + 2 )
    real(kind=dp), parameter     :: a4 = a3    * ( sigma + 3 )
    real(kind=dp), parameter     :: a5 = a4    * ( sigma + 4 )

    C5    = ( 04.0_dp * a5  * sin(lat)**4  - &
              20.0_dp * a4  * sin(lat)**2  + &
              15.0_dp * a3) * sin(lat)     / 15.0_dp

    C5p   = ( 04.0_dp * a5  * sin(lat)**4  - &
              12.0_dp * a4  * sin(lat)**2  + &
              03.0_dp * a3) * cos(lat)     / 3.0_dp

    psi   = amp * cos(lat) ** (sigma) * C5

    dpsi  = amp * cos(lat) ** (sigma) * ( - sigma * tan(lat) * C5 + C5p )

  end subroutine getPsi

!========================================================================
! latitude-dependent amplitudes
!========================================================================
  subroutine getAmplitudes(lat,uTilde,vTilde,hTilde,waveFlag)

    implicit none

    integer,       intent(in)    :: waveFlag          ! -1=WIG, 0=Rossby, 1=EIG
    real(kind=dp), intent(in)    :: lat(:)            ! latitude (rad)
    real(kind=dp), intent(inout) :: uTilde(:)         ! see Eq. 18 in text
    real(kind=dp), intent(inout) :: vTilde(:)         ! see Eqs. 12a and 15b in text
    real(kind=dp), intent(inout) :: hTilde(:)         ! see Eqs. 12b and 15a in text

    real(kind=dp)                :: Kp(size(lat))     ! see Eq. 13 in text
    real(kind=dp)                :: Km(size(lat))     ! see Eq. 13 in text
    real(kind=dp)                :: W1(size(lat))     ! see Eq. 16 in text
    real(kind=dp)                :: W2(size(lat))     ! see Eq. 17 in text
    real(kind=dp)                :: C                 ! phase speed (rad/sec)
    real(kind=dp)                :: psi(size(lat))    ! see Eq. 13 in text
    real(kind=dp)                :: dpsi(size(lat))   ! meridional derivative of psi

    real(kind=dp),    parameter  :: r2 = 0.5_dp
    real(kind=dp),    parameter  :: o2 = 2.0_dp * omega

    call getPhaseSpeed(C,waveFlag)
    call getPsi(lat,psi,dpsi)

    select case (waveFlag)
    case(0)
       Kp = (g*H0 + a**2 * C**2 * cos(lat)**2) / (C*cos(lat))
       Km = (g*H0 - a**2 * C**2 * cos(lat)**2) / (C*cos(lat))

       vTilde = ( o2    * abs(Km)  / cos(lat)**2  )**r2 * psi
       hTilde = ( o2    * abs(Km)  * a**2 * H0**2 )**r2 / Km * &
                ( dpsi  + tan(lat) * ( r2 * Kp/Km - o2 / C)  * psi )
       uTilde = ( o2    * sin(lat) / C ) * vTilde + &
                ( g / a / cos(lat) / C ) * hTilde

    case(-1,1)
       W1 = ( (k*C)**2 - o2**2 * sin(lat)**2 ) / ( C  * cos(lat) )
       W2 = ( (k*C)**2 - o2**2 * sin(lat)**2 ) / ( o2 * cos(lat)**2 )

       vTilde = ( o2    * abs(W1)  * g*H0 / cos(lat)**2 )**r2 / W1 * &
                ( dpsi  + tan(lat) * ( r2 - o2 / W2 + o2 / C ) * psi )
       hTilde = ( o2    * abs(W1)  * a**2 * H0 / g )**r2 * psi
       uTilde = ( o2    * sin(lat) / C ) * vTilde + &
                ( g / a / cos(lat) / C ) * hTilde
    end select

    vTilde = k * vTilde

  end subroutine getAmplitudes

!========================================================================
! u,v,h fields
!========================================================================
  subroutine getFields(latCell,lonCell,latEdge,lonEdge,time,u,v,h,waveFlag)

    implicit none

    integer,       intent(in)    :: waveFlag            ! -1=WIG, 0=Rossby, 1=EIG
    real(kind=dp), intent(in)    :: latCell(:)          ! latitude (rad) for mass
    real(kind=dp), intent(in)    :: lonCell(:)          ! longitude (rad) for mass
    real(kind=dp), intent(in)    :: latEdge(:)          ! latitude (rad) for wind
    real(kind=dp), intent(in)    :: lonEdge(:)          ! longitude (rad) for wind
    real(kind=dp), intent(in)    :: time(:)             ! time (sec)
    real(kind=dp), intent(inout) :: u(:,:)              ! zonal velocity field
    real(kind=dp), intent(inout) :: v(:,:)              ! meridional velocity field
    real(kind=dp), intent(inout) :: h(:,:)              ! free-surface height anomaly field

    real(kind=dp)                :: uTilde(size(latEdge))   ! see Eq. 18 in text
    real(kind=dp)                :: vTilde(size(latEdge))   ! see Eqs. 12a and 15b in text
    real(kind=dp)                :: hTilde(size(latCell))   ! see Eqs. 12b and 15a in text
    real(kind=dp)                :: C                       ! phase speed (rad/sec)
    
    real(kind=dp)                :: blankCell(size(latCell))   ! working array
    real(kind=dp)                :: blankEdge(size(latEdge))   ! working array
    integer                      :: n,iCell,iEdge

    call getPhaseSpeed(C,waveFlag)
    call getAmplitudes(latCell, blankCell, blankCell, hTilde   , waveFlag)
    call getAmplitudes(latEdge, uTilde   , vTilde   , blankEdge, waveFlag)

    do n = 1, size(time)
      do iCell = lbound(h, 1), ubound(h, 1)
        h(iCell,n) = hTilde(iCell) * cos( k * lonCell(iCell) - k * C * time(n) )
      end do
    end do
    h = h + H0
    
    do n = 1, size(time)
      do iEdge = lbound(u, 1), ubound(u, 1)
        u(iEdge,n) = uTilde(iEdge) * cos( k * lonEdge(iEdge) - k * C * time(n)               )
        v(iEdge,n) = vTilde(iEdge) * cos( k * lonEdge(iEdge) - k * C * time(n) - 0.5_dp * pi )
      end do
    end do    

  end subroutine getFields

!========================================================================
! surface pressure
!========================================================================
  subroutine getSurfacePressure(h,Ps)

    implicit none

    real(kind=dp), intent(in)    :: h(:,:,:)        ! free-surface height anomaly field
    real(kind=dp), intent(inout) :: Ps(:,:,:)       ! surface pressure, see Eq. 23 in text

    Ps  = P0*(1.0+h/H0)

  end subroutine getSurfacePressure

!========================================================================
! initial temperature
!========================================================================
  subroutine getInitialTemperature(T)

    implicit none

    real(kind=dp), intent(inout) :: T(:,:,:)        ! initial temperature field

    T = g*H0/Rd

  end subroutine getInitialTemperature

end module shallow_water_waves_test