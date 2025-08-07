module equations_module
  use constants_module, only: dp
  use variables_module, only: nlon, nlat, dx, dy, g, f0, h0, h1, pi
  implicit none
contains

  subroutine init_height(h, x, y)
    real(dp), intent(out) :: h(nlon,nlat)
    real(dp), intent(in) :: x(nlon), y(nlat)
    real(dp) :: x0, y0, r0, dist
    integer :: i,j
    x0 = 1.5d0*dx*nlon
    y0 = 0.5d0*dy*nlat
    r0 = dy*nlat/3.d0
    do j=1,nlat
       do i=1,nlon
          h(i,j) = h0
          dist = sqrt( (x(i)-x0)**2 + (y(j)-y0)**2 )
          if (dist < r0) then
             h(i,j) = h0 + 0.5d0*h1*(1.d0+cos(pi*dist/r0))
          end if
       end do
    end do
  end subroutine init_height

  subroutine velocity_field(u,v,x,y,alpha)
    real(dp), intent(out) :: u(nlon,nlat), v(nlon,nlat+1)
    real(dp), intent(in) :: x(nlon), y(nlat), alpha
    integer :: i,j
    u = 0.d0
    v = 0.d0
    v(:,1) = 0.d0
    v(:,nlat+1) = 0.d0
  end subroutine velocity_field

  subroutine analytic_height(ha, x, y, t, alpha)
    real(dp), intent(out) :: ha(nlon,nlat)
    real(dp), intent(in) :: x(nlon), y(nlat), t, alpha
    call init_height(ha, x, y)
  end subroutine analytic_height

  subroutine rhs(h, u, v, dhdt, dudt, dvdt, no_momentum_tendency)
    real(dp), intent(in) :: h(nlon,nlat), u(nlon,nlat), v(nlon,nlat+1)
    real(dp), intent(out) :: dhdt(nlon,nlat)
    real(dp), intent(out) :: dudt(nlon,nlat), dvdt(nlon,nlat+1)
    logical, intent(in), optional :: no_momentum_tendency
    integer :: i,j,ip1,im1,jp1,jm1
    real(dp) :: fe,fw,fn,fs,ue,uw,vn,vs
    real(dp) :: h_e, h_w, h_n, h_s, v_avg, u_avg
    real(dp) :: fcor

    ! continuity equation
    do j=1,nlat
       jp1 = min(j+1,nlat)
       jm1 = max(j-1,1)
       do i=1,nlon
          ip1 = mod(i,nlon)+1
          im1 = mod(i-2+nlon,nlon)+1
          ue = u(ip1,j)
          uw = u(i,j)
          if (ue > 0.d0) then
             fe = ue * h(i,j)
          else
             fe = ue * h(ip1,j)
          end if
          if (uw > 0.d0) then
             fw = uw * h(im1,j)
          else
             fw = uw * h(i,j)
          end if
          vn = v(i,j+1)
          vs = v(i,j)
          if (vn > 0.d0) then
             fn = vn * h(i,j)
          else
             fn = vn * h(i,jp1)
          end if
          if (vs > 0.d0) then
             fs = vs * h(i,jm1)
          else
             fs = vs * h(i,j)
          end if
          dhdt(i,j) = -((fe - fw)/dx + (fn - fs)/dy)
       end do
    end do

    if (present(no_momentum_tendency)) then
       if (no_momentum_tendency) then
          dudt = 0.d0
          dvdt = 0.d0
          return
       end if
    end if

    fcor = f0

    ! zonal momentum
    do j=1,nlat
       jp1 = min(j+1,nlat)
       jm1 = max(j-1,1)
       do i=1,nlon
          ip1 = mod(i,nlon)+1
          im1 = mod(i-2+nlon,nlon)+1
          h_e = h(i,j)
          h_w = h(im1,j)
          v_avg = 0.25d0*(v(im1,j) + v(im1,j+1) + v(i,j) + v(i,j+1))
          dudt(i,j) = - u(i,j) * (u(ip1,j) - u(im1,j)) / (2.d0*dx) &
                      - v_avg * (u(i,jp1) - u(i,jm1)) / (2.d0*dy) &
                      - g * (h_e - h_w) / dx &
                      + fcor * v_avg
       end do
    end do

    ! meridional momentum
    do j=2,nlat
       jp1 = min(j+1,nlat)
       jm1 = j-1
       do i=1,nlon
          ip1 = mod(i,nlon)+1
          im1 = mod(i-2+nlon,nlon)+1
          h_n = h(i,j)
          h_s = h(i,jm1)
          u_avg = 0.25d0*(u(i,jm1) + u(ip1,jm1) + u(i,j) + u(ip1,j))
          dvdt(i,j) = - u_avg * (v(ip1,j) - v(im1,j)) / (2.d0*dx) &
                      - v(i,j) * (v(i,jp1) - v(i,jm1)) / (2.d0*dy) &
                      - g * (h_n - h_s) / dy &
                      - fcor * u_avg
       end do
    end do
    dvdt(:,1) = 0.d0
    dvdt(:,nlat+1) = 0.d0
  end subroutine rhs

end module equations_module
