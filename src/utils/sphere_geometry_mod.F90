module sphere_geometry_mod

  use params_mod
  use math_mod

  implicit none

  private

  public cartesian_transform
  public inverse_cartesian_transform
  public rotation_transform
  public inverse_rotation_transform
  public calc_distance
  public cross_product
  public norm_vector
  public calc_sphere_angle
  public calc_plane_angle
  public calc_arc_length
  public point_type

  type point_type
    real(real_kind) lon
    real(real_kind) lat
    real(real_kind) x
    real(real_kind) y
    real(real_kind) z
  contains
    procedure :: copy_coord => point_copy_coord
  end type point_type

  interface cartesian_transform
    module procedure cartesian_transform_1
    module procedure cartesian_transform_2
  end interface cartesian_transform

  interface inverse_cartesian_transform
    module procedure inverse_cartesian_transform_1
    module procedure inverse_cartesian_transform_2
  end interface inverse_cartesian_transform

contains

  subroutine cartesian_transform_1(lon, lat, x, y, z)

    real(real_kind), intent(in)  :: lon, lat
    real(real_kind), intent(out) :: x, y, z

    real(real_kind) cos_lat

    cos_lat = cos(lat)
    x = cos_lat * cos(lon)
    y = cos_lat * sin(lon)
    z = sin(lat)

  end subroutine cartesian_transform_1

  subroutine cartesian_transform_2(point)

    class(point_type), intent(inout) :: point

    real(real_kind) cos_lat

    cos_lat = cos(point%lat)
    point%x = cos_lat * cos(point%lon)
    point%y = cos_lat * sin(point%lon)
    point%z = sin(point%lat)

  end subroutine cartesian_transform_2

  subroutine inverse_cartesian_transform_1(lon, lat, x, y, z)

    real(real_kind), intent(out) :: lon, lat
    real(real_kind), intent(in)  :: x, y, z

    lon = atan2(y, x)
    lat = asin(z)

    if (lon < 0.0d0) lon = lon + pi2

  end subroutine inverse_cartesian_transform_1

  subroutine inverse_cartesian_transform_2(point)

    class(point_type), intent(inout) :: point

    point%lon = atan2(point%y, point%x)
    point%lat = asin(point%z)

    if (point%lon < 0.0d0) point%lon = point%lon + pi2

  end subroutine inverse_cartesian_transform_2

  ! ************************************************************************** !
  ! Rotation transform                                                         !
  ! Purpose:                                                                   !
  !   Calculate the rotating transformation and its inverse of the original    !
  !   coordinate system (lonO,latO) to the rotated one (lonR, latR) with       !
  !   the north pole (lonP,latP) defined at the original coordinate system.    !
  ! ************************************************************************** !

  subroutine rotation_transform(lonP, latP, lonO, latO, lonR, latR)

    real(real_kind), intent(in) :: lonP, latP ! Rotated pole coordinate
    real(real_kind), intent(in) :: lonO, latO ! Original coordinate
    real(real_kind), intent(out), optional :: lonR, latR ! Rotated coordinate

    real(real_kind) tmp1, tmp2, tmp3, dlon

    dlon = lonO - lonP
    if (present(lonR)) then
        tmp1 = cos(latO) * sin(dlon)
        tmp2 = cos(latO) * sin(latP) * cos(dlon) - cos(latP) * sin(latO)
        lonR = atan2(tmp1, tmp2)
        if (lonR < 0.0d0) lonR = PI2 + lonR
    end if
    if (present(latR)) then
        tmp1 = sin(latO) * sin(latP)
        tmp2 = cos(latO) * cos(latP) * cos(dlon)
        tmp3 = tmp1 + tmp2
        tmp3 = min(max(tmp3, -1.0d0), 1.0d0)
        latR = asin(tmp3)
    end if

  end subroutine rotation_transform

  subroutine inverse_rotation_transform(lonP, latP, lonO, latO, lonR, latR)

      real(real_kind), intent(in)  :: lonP, latP ! Rotated pole coordinate
      real(real_kind), intent(out) :: lonO, latO ! Original coordinate
      real(real_kind), intent(in)  :: lonR, latR ! Rotated coordinate

      real(real_kind) sinLonR, cosLonR, sinLatR, cosLatR, sinLatP, cosLatP
      real(real_kind) tmp1, tmp2, tmp3

      sinLonR = sin(lonR)
      cosLonR = cos(lonR)
      sinLatR = sin(latR)
      cosLatR = cos(latR)
      sinLatP = sin(latP)
      cosLatP = cos(LatP)

      tmp1 = cosLatR * sinLonR
      tmp2 = sinLatR * cosLatP + cosLatR * cosLonR * sinLatP
      ! This trick is due to the inaccuracy of trigonometry calculation.
      if (abs(tmp2) < eps) tmp2 = 0.0d0
      lonO = atan2(tmp1, tmp2)
      lonO = lonP + lonO
      if (lonO > PI2) lonO = lonO - PI2
      tmp1 = sinLatR * sinLatP
      tmp2 = cosLatR * cosLatP * cosLonR
      tmp3 = tmp1 - tmp2
      tmp3 = min(max(tmp3, -1.0d0), 1.0d0)
      latO = asin(tmp3)

  end subroutine inverse_rotation_transform

  real(real_kind) function calc_distance(lon1, lat1, lon2, lat2) result(res)

    real(real_kind), intent(in) :: lon1
    real(real_kind), intent(in) :: lat1
    real(real_kind), intent(in) :: lon2
    real(real_kind), intent(in) :: lat2

    res = radius * acos(min(1.0d0, max(-1.0d0, sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))))

  end function calc_distance

  function norm_vector(x) result(res)

    real(real_kind), intent(in) :: x(:)
    real(real_kind) res(size(x))

    real(real_kind) n

    n = sqrt(sum(x * x))
    if (n /= 0) then
      res = x / n
    else
      res = x
    end if

  end function norm_vector

  ! Calculate the dihedra angle between plane AB and plane AC.

  real(real_kind) function calc_sphere_angle(a, b, c) result(res)

    real(real_kind), intent(in) :: a(3)
    real(real_kind), intent(in) :: b(3)
    real(real_kind), intent(in) :: c(3)

    real(real_kind) nab(3) ! Normal vector of plane AB
    real(real_kind) nac(3) ! Normal vector of plane AC

    nab = norm_vector(cross_product(a, b))
    nac = norm_vector(cross_product(a, c))
    res = acos(max(min(dot_product(nab, nac), 1.0d0), -1.0d0))

    ! Judge the cyclic direction with respect to point A to handle obtuse angle.
    if (dot_product(cross_product(b - a, c - a), a) < 0.0) res = -res

  end function calc_sphere_angle

  ! Calculate the angle between vector AB and AC on a plane.

  real(real_kind) function calc_plane_angle(a, b, c, n) result(res)

    real(real_kind), intent(in) :: a(:) ! Reference point
    real(real_kind), intent(in) :: b(:) ! Point 1
    real(real_kind), intent(in) :: c(:) ! Point 2
    real(real_kind), intent(in) :: n(:) ! Normal vector

    real(real_kind) ab(size(a)), ac(size(a))

    ab = norm_vector(b - a)
    ac = norm_vector(c - a)

    res = acos(max(min(dot_product(ab, ac), 1.0d0), -1.0d0))

    ! Judge the cyclic direction with respect to given normal vector to handle obtuse angle.
    if (dot_product(cross_product(ab, ac), n) < 0.0d0) res = -res

  end function calc_plane_angle

  ! Calculate the great circle arc length from A to B by assuming A and B are on the unit sphere surface.

  real(real_kind) function calc_arc_length(a, b) result(res)

    real(real_kind), intent(in) :: a(3)
    real(real_kind), intent(in) :: b(3)

    res = acos(max(min(dot_product(a, b), 1.0d0), -1.0d0))

  end function calc_arc_length

  subroutine point_copy_coord(a, b)

    class(point_type), intent(inout) :: a
    class(point_type), intent(in)    :: b

    a%lon = b%lon
    a%lat = b%lat
    a%x   = b%x
    a%y   = b%y
    a%z   = b%z

  end subroutine point_copy_coord

end module sphere_geometry_mod
