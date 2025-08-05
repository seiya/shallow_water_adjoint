module rk4_module
  use constants_module, only: dp
  use variables_module, only: nlon, nlat, dt
  use equations_module, only: rhs
  implicit none
contains

  !$FAD CONSTANT_VARS: lat
  subroutine rk4_step(h,u,v,hn,un,vn,lat)
    real(dp), intent(in) :: h(nlon,nlat), u(nlon,nlat), v(nlon,nlat+1)
    real(dp), intent(out) :: hn(nlon,nlat), un(nlon,nlat), vn(nlon,nlat+1)
    real(dp), intent(in) :: lat(nlat)
    real(dp) :: k1h(nlon,nlat), k2h(nlon,nlat), k3h(nlon,nlat), k4h(nlon,nlat)
    real(dp) :: k1u(nlon,nlat), k2u(nlon,nlat), k3u(nlon,nlat), k4u(nlon,nlat)
    real(dp) :: k1v(nlon,nlat+1), k2v(nlon,nlat+1), k3v(nlon,nlat+1), k4v(nlon,nlat+1)
    real(dp) :: htmp(nlon,nlat), utmp(nlon,nlat), vtmp(nlon,nlat+1)
    call rhs(h, u, v, k1h, k1u, k1v, lat)
    htmp = h + 0.5d0*dt*k1h
    utmp = u + 0.5d0*dt*k1u
    vtmp = v + 0.5d0*dt*k1v
    call rhs(htmp, utmp, vtmp, k2h, k2u, k2v, lat)
    htmp = h + 0.5d0*dt*k2h
    utmp = u + 0.5d0*dt*k2u
    vtmp = v + 0.5d0*dt*k2v
    call rhs(htmp, utmp, vtmp, k3h, k3u, k3v, lat)
    htmp = h + dt*k3h
    utmp = u + dt*k3u
    vtmp = v + dt*k3v
    call rhs(htmp, utmp, vtmp, k4h, k4u, k4v, lat)
    hn = h + dt*(k1h + 2.d0*k2h + 2.d0*k3h + k4h)/6.d0
    un = u + dt*(k1u + 2.d0*k2u + 2.d0*k3u + k4u)/6.d0
    vn = v + dt*(k1v + 2.d0*k2v + 2.d0*k3v + k4v)/6.d0
    vn(:,1) = 0.d0
    vn(:,nlat+1) = 0.d0
  end subroutine rk4_step

end module rk4_module
