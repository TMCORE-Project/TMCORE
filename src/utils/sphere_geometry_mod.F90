module sphere_geometry_mod

  use params_mod

  implicit none

contains

  real(real_kind) function distance(lon1, lat1, lon2, lat2) result(res)

    real(real_kind), intent(in) :: lon1
    real(real_kind), intent(in) :: lat1
    real(real_kind), intent(in) :: lon2
    real(real_kind), intent(in) :: lat2

    res = radius * acos(min(1.0d0, max(-1.0d0, sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))))

  end function distance

end module sphere_geometry_mod