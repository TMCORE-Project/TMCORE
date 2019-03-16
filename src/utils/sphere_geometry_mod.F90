module sphere_geometry_mod

  use params_mod
  use math_mod

  implicit none

  private

  public calc_distance
  public cross_product
  public norm_vector
  public calc_sphere_angle
  public calc_plane_angle
  public calc_arc_length

contains

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

end module sphere_geometry_mod