module rk4_module
  use constants_module, only: dp
  use variables_module, only: nx, ny, dt, ihalo, exchange_halo_x
  use equations_module, only: rhs
  implicit none
contains

  !$FAD CONSTANT_VARS: no_momentum_tendency
  subroutine rk4_step(h,u,v,hn,un,vn,no_momentum_tendency)
    real(dp), intent(inout) :: h(1-ihalo:nx+ihalo,ny), u(1-ihalo:nx+ihalo,ny), v(1-ihalo:nx+ihalo,ny+1)
    real(dp), intent(out) :: hn(1-ihalo:nx+ihalo,ny), un(1-ihalo:nx+ihalo,ny), vn(1-ihalo:nx+ihalo,ny+1)
    logical, intent(in), optional :: no_momentum_tendency
    real(dp) :: k1h(1-ihalo:nx+ihalo,ny), k2h(1-ihalo:nx+ihalo,ny), &
                 k3h(1-ihalo:nx+ihalo,ny), k4h(1-ihalo:nx+ihalo,ny)
    real(dp) :: k1u(1-ihalo:nx+ihalo,ny), k2u(1-ihalo:nx+ihalo,ny), &
                 k3u(1-ihalo:nx+ihalo,ny), k4u(1-ihalo:nx+ihalo,ny)
    real(dp) :: k1v(1-ihalo:nx+ihalo,ny+1), k2v(1-ihalo:nx+ihalo,ny+1), &
                 k3v(1-ihalo:nx+ihalo,ny+1), k4v(1-ihalo:nx+ihalo,ny+1)
    real(dp) :: htmp(1-ihalo:nx+ihalo,ny), utmp(1-ihalo:nx+ihalo,ny), &
                 vtmp(1-ihalo:nx+ihalo,ny+1)
    logical :: skip_momentum
    skip_momentum = .false.
    if (present(no_momentum_tendency)) skip_momentum = no_momentum_tendency

    if (skip_momentum) then
       call rhs(h, u, v, k1h, k1u, k1v, no_momentum_tendency=.true.)
       htmp = h + 0.5d0*dt*k1h
       call rhs(htmp, u, v, k2h, k2u, k2v, no_momentum_tendency=.true.)
       htmp = h + 0.5d0*dt*k2h
       call rhs(htmp, u, v, k3h, k3u, k3v, no_momentum_tendency=.true.)
       htmp = h + dt*k3h
       call rhs(htmp, u, v, k4h, k4u, k4v, no_momentum_tendency=.true.)
       hn = h + dt*(k1h + 2.d0*k2h + 2.d0*k3h + k4h)/6.d0
       un = u
       vn = v
       vn(:,1) = 0.d0
       vn(:,ny+1) = 0.d0
       call exchange_halo_x(hn)
       call exchange_halo_x(un)
       call exchange_halo_x(vn)
       return
    end if

    call rhs(h, u, v, k1h, k1u, k1v, no_momentum_tendency=.false.)
    htmp = h + 0.5d0*dt*k1h
    utmp = u + 0.5d0*dt*k1u
    vtmp = v + 0.5d0*dt*k1v
    call rhs(htmp, utmp, vtmp, k2h, k2u, k2v, no_momentum_tendency=.false.)
    htmp = h + 0.5d0*dt*k2h
    utmp = u + 0.5d0*dt*k2u
    vtmp = v + 0.5d0*dt*k2v
    call rhs(htmp, utmp, vtmp, k3h, k3u, k3v, no_momentum_tendency=.false.)
    htmp = h + dt*k3h
    utmp = u + dt*k3u
    vtmp = v + dt*k3v
    call rhs(htmp, utmp, vtmp, k4h, k4u, k4v, no_momentum_tendency=.false.)
    hn = h + dt*(k1h + 2.d0*k2h + 2.d0*k3h + k4h)/6.d0
    un = u + dt*(k1u + 2.d0*k2u + 2.d0*k3u + k4u)/6.d0
    vn = v + dt*(k1v + 2.d0*k2v + 2.d0*k3v + k4v)/6.d0
    vn(:,1) = 0.d0
    vn(:,ny+1) = 0.d0
    call exchange_halo_x(hn)
    call exchange_halo_x(un)
    call exchange_halo_x(vn)
  end subroutine rk4_step

end module rk4_module
