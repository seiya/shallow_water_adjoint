module rk4_module
  use constants_module, only: dp
  use variables_module, only: nlon, nlat, dt
  use equations_module, only: rhs
  implicit none
contains
  subroutine rk4_step(h,hn,u,v,lat)
    real(dp), intent(in) :: h(nlon,nlat), u(nlon+1,nlat), v(nlon,nlat+1), lat(nlat)
    real(dp), intent(out) :: hn(nlon,nlat)
    real(dp) :: k1(nlon,nlat), k2(nlon,nlat)
    real(dp) :: k3(nlon,nlat), k4(nlon,nlat)
    real(dp) :: htmp(nlon,nlat)
    call rhs(h, k1, u, v, lat)
    htmp = h + 0.5d0*dt*k1
    call rhs(htmp, k2, u, v, lat)
    htmp = h + 0.5d0*dt*k2
    call rhs(htmp, k3, u, v, lat)
    htmp = h + dt*k3
    call rhs(htmp, k4, u, v, lat)
    hn = h + dt*(k1 + 2.d0*k2 + 2.d0*k3 + k4)/6.d0
  end subroutine rk4_step
end module rk4_module
