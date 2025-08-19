module rk4_module
  use constants_module, only: dp
  use variables_module, only: nx, ny, dt, ihalo, is, ie, js, je, exchange_halo
  use equations_module, only: rhs
  use mpi_decomp_module, only: istart, iend, jstart, jend
  implicit none
contains

  subroutine rk4_step(h,u,v,hn,un,vn,no_momentum_tendency)
    real(dp), intent(in) :: h(is:ie,js:je), u(is:ie,js:je), v(is:ie,js:jend+1)
    real(dp), intent(out) :: hn(is:ie,js:je), un(is:ie,js:je), vn(is:ie,js:jend+1)
    logical, intent(in), optional :: no_momentum_tendency
    real(dp), allocatable :: k1h(:,:), k2h(:,:), k3h(:,:), k4h(:,:)
    real(dp), allocatable :: k1u(:,:), k2u(:,:), k3u(:,:), k4u(:,:)
    real(dp), allocatable :: k1v(:,:), k2v(:,:), k3v(:,:), k4v(:,:)
    real(dp), allocatable :: htmp(:,:), utmp(:,:), vtmp(:,:)
    logical :: skip_momentum
    skip_momentum = .false.
    if (present(no_momentum_tendency)) skip_momentum = no_momentum_tendency

    allocate(k1h(is:ie,js:je), k2h(is:ie,js:je), k3h(is:ie,js:je), k4h(is:ie,js:je))
    allocate(k1u(is:ie,js:je), k2u(is:ie,js:je), k3u(is:ie,js:je), k4u(is:ie,js:je))
    allocate(k1v(is:ie,js:jend+1), k2v(is:ie,js:jend+1), k3v(is:ie,js:jend+1), k4v(is:ie,js:jend+1))
    allocate(htmp(is:ie,js:je), utmp(is:ie,js:je), vtmp(is:ie,js:jend+1))

    if (skip_momentum) then
       call rhs(h, u, v, k1h, k1u, k1v, no_momentum_tendency=.true.)
       htmp(istart:iend,jstart:jend) = h(istart:iend,jstart:jend) + 0.5d0*dt*k1h(istart:iend,jstart:jend)
       call exchange_halo(htmp)
       call rhs(htmp, u, v, k2h, k2u, k2v, no_momentum_tendency=.true.)
       htmp(istart:iend,jstart:jend) = h(istart:iend,jstart:jend) + 0.5d0*dt*k2h(istart:iend,jstart:jend)
       call exchange_halo(htmp)
       call rhs(htmp, u, v, k3h, k3u, k3v, no_momentum_tendency=.true.)
       htmp(istart:iend,jstart:jend) = h(istart:iend,jstart:jend) + dt*k3h(istart:iend,jstart:jend)
       call exchange_halo(htmp)
       call rhs(htmp, u, v, k4h, k4u, k4v, no_momentum_tendency=.true.)
       hn(istart:iend,jstart:jend) = h(istart:iend,jstart:jend) + dt*(k1h(istart:iend,jstart:jend) + 2.d0*k2h(istart:iend,jstart:jend) + 2.d0*k3h(istart:iend,jstart:jend) + k4h(istart:iend,jstart:jend))/6.d0
       un(istart:iend,jstart:jend) = u(istart:iend,jstart:jend)
       vn(istart:iend,jstart:jend) = v(istart:iend,jstart:jend)
       call exchange_halo(hn)
       call exchange_halo(un)
       call exchange_halo(vn)
       return
    end if

    call rhs(h, u, v, k1h, k1u, k1v, no_momentum_tendency=.false.)
    htmp(istart:iend,jstart:jend) = h(istart:iend,jstart:jend) + 0.5d0*dt*k1h(istart:iend,jstart:jend)
    utmp(istart:iend,jstart:jend) = u(istart:iend,jstart:jend) + 0.5d0*dt*k1u(istart:iend,jstart:jend)
    vtmp(istart:iend,jstart:jend) = v(istart:iend,jstart:jend) + 0.5d0*dt*k1v(istart:iend,jstart:jend)
    call exchange_halo(htmp)
    call exchange_halo(utmp)
    call exchange_halo(vtmp)
    call rhs(htmp, utmp, vtmp, k2h, k2u, k2v, no_momentum_tendency=.false.)
    htmp(istart:iend,jstart:jend) = h(istart:iend,jstart:jend) + 0.5d0*dt*k2h(istart:iend,jstart:jend)
    utmp(istart:iend,jstart:jend) = u(istart:iend,jstart:jend) + 0.5d0*dt*k2u(istart:iend,jstart:jend)
    vtmp(istart:iend,jstart:jend) = v(istart:iend,jstart:jend) + 0.5d0*dt*k2v(istart:iend,jstart:jend)
    call exchange_halo(htmp)
    call exchange_halo(utmp)
    call exchange_halo(vtmp)
    call rhs(htmp, utmp, vtmp, k3h, k3u, k3v, no_momentum_tendency=.false.)
    htmp(istart:iend,jstart:jend) = h(istart:iend,jstart:jend) + dt*k3h(istart:iend,jstart:jend)
    utmp(istart:iend,jstart:jend) = u(istart:iend,jstart:jend) + dt*k3u(istart:iend,jstart:jend)
    vtmp(istart:iend,jstart:jend) = v(istart:iend,jstart:jend) + dt*k3v(istart:iend,jstart:jend)
    call exchange_halo(htmp)
    call exchange_halo(utmp)
    call exchange_halo(vtmp)
    call rhs(htmp, utmp, vtmp, k4h, k4u, k4v, no_momentum_tendency=.false.)
    hn(istart:iend,jstart:jend) = h(istart:iend,jstart:jend) + dt*(k1h(istart:iend,jstart:jend) + 2.d0*k2h(istart:iend,jstart:jend) + 2.d0*k3h(istart:iend,jstart:jend) + k4h(istart:iend,jstart:jend))/6.d0
    un(istart:iend,jstart:jend) = u(istart:iend,jstart:jend) + dt*(k1u(istart:iend,jstart:jend) + 2.d0*k2u(istart:iend,jstart:jend) + 2.d0*k3u(istart:iend,jstart:jend) + k4u(istart:iend,jstart:jend))/6.d0
    vn(istart:iend,jstart:jend) = v(istart:iend,jstart:jend) + dt*(k1v(istart:iend,jstart:jend) + 2.d0*k2v(istart:iend,jstart:jend) + 2.d0*k3v(istart:iend,jstart:jend) + k4v(istart:iend,jstart:jend))/6.d0
    call exchange_halo(hn)
    call exchange_halo(un)
    call exchange_halo(vn)
  end subroutine rk4_step

end module rk4_module
