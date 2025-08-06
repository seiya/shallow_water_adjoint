module equations_module
  use constants_module, only: dp
  use variables_module, only: nlon, nlat, pi, radius, dlon, dlat, h0, h1, omega, g
  implicit none

contains

  !$FAD CONSTANT_VARS: lon, lat
  subroutine init_height(h, lon, lat)
    real(dp), intent(out) :: h(nlon,nlat)
    real(dp), intent(in) :: lon(nlon), lat(nlat)
    real(dp) :: lambda0, phi0, r0, dist
    integer :: i,j
    lambda0 = 3.d0*pi/2.d0
    phi0 = 0.d0
    r0 = radius/3.d0
    do j=1,nlat
       do i=1,nlon
          h(i,j) = h0
          call gc_distance(lon(i),lat(j),lambda0,phi0,dist)
          if (dist < r0) then
             h(i,j) = h0 + 0.5d0*h1*(1.d0+cos(pi*dist/r0))
          end if
       end do
    end do
  end subroutine init_height

  !$FAD CONSTANT_VARS: lon, lat, alpha
  subroutine velocity_field(u,v,lon,lat,alpha)
    real(dp), intent(out) :: u(nlon,nlat), v(nlon,nlat+1)
    real(dp), intent(in) :: lon(nlon), lat(nlat), alpha
    real(dp) :: u0, lon_edge
    integer :: i,j
    u0 = omega*radius
    do j=1,nlat
       do i=1,nlon
          lon_edge = (i-1)*dlon
          u(i,j) = u0*(cos(lat(j))*cos(alpha) + sin(lat(j))*cos(lon_edge)*sin(alpha))
       end do
    end do
    do j=1,nlat+1
       do i=1,nlon
          v(i,j) = -u0*sin(lon(i))*sin(alpha)
       end do
    end do
    v(:,1) = 0.d0
    v(:,nlat+1) = 0.d0
  end subroutine velocity_field

  !$FAD SKIP
  subroutine analytic_height(ha, lon, lat, t, alpha)
    real(dp), intent(out) :: ha(nlon,nlat)
    real(dp), intent(in) :: lon(nlon), lat(nlat), t, alpha
    real(dp) :: lonr, latr, theta
    integer :: i,j
    theta = -omega*t
    do j=1,nlat
       do i=1,nlon
          call rotate_point(lon(i),lat(j),alpha,theta,lonr,latr)
          ha(i,j) = initial_at_point(lonr,latr)
       end do
    end do
  end subroutine analytic_height

  !$FAD CONSTANT_VARS: lon, lat
  function initial_at_point(lon,lat) result(hp)
    real(dp), intent(in) :: lon, lat
    real(dp) :: hp, lambda0, phi0, r0, dist
    lambda0 = 3.d0*pi/2.d0
    phi0 = 0.d0
    r0 = radius/3.d0
    hp = h0
    call gc_distance(lon,lat,lambda0,phi0,dist)
    if (dist < r0) hp = h0 + 0.5d0*h1*(1.d0+cos(pi*dist/r0))
  end function initial_at_point

  !$FAD SKIP
  subroutine gc_distance(lon1,lat1,lon2,lat2,dist)
    real(dp), intent(in) :: lon1,lat1,lon2,lat2
    real(dp), intent(out) :: dist
    dist = radius*acos( sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon1-lon2) )
  end subroutine gc_distance

  !$FAD SKIP
  subroutine rotate_point(lon,lat,alpha,theta,lonp,latp)
    real(dp), intent(in) :: lon, lat, alpha, theta
    real(dp), intent(out) :: lonp, latp
    real(dp) :: x,y,z,xp,yp,zp
    real(dp) :: n1,n2,n3,dot,c1,c2,c3,ct,st
    x = cos(lat)*cos(lon)
    y = cos(lat)*sin(lon)
    z = sin(lat)
    n1 = -sin(alpha)
    n2 = 0.d0
    n3 = cos(alpha)
    ct = cos(theta)
    st = sin(theta)
    dot = n1*x + n2*y + n3*z
    c1 = n2*z - n3*y
    c2 = n3*x - n1*z
    c3 = n1*y - n2*x
    xp = x*ct + c1*st + n1*dot*(1.d0-ct)
    yp = y*ct + c2*st + n2*dot*(1.d0-ct)
    zp = z*ct + c3*st + n3*dot*(1.d0-ct)
    lonp = atan2(yp,xp)
    if (lonp < 0.d0) lonp = lonp + 2.d0*pi
    latp = asin(zp)
  end subroutine rotate_point

  !$FAD CONSTANT_VARS: lat, no_momentum_tendency
  subroutine rhs(h, u, v, dhdt, dudt, dvdt, lat, no_momentum_tendency)
    real(dp), intent(in) :: h(nlon,nlat), u(nlon,nlat), v(nlon,nlat+1)
    real(dp), intent(out) :: dhdt(nlon,nlat)
    real(dp), intent(out) :: dudt(nlon,nlat), dvdt(nlon,nlat+1)
    real(dp), intent(in) :: lat(nlat)
    logical, intent(in), optional :: no_momentum_tendency
    integer :: i,j,ip1,im1,jp1,jm1
    real(dp) :: fe,fw,fn,fs,ue,uw,vn,vs
    real(dp) :: fcor, h_e, h_w, h_n, h_s, v_avg, u_avg
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
          vn = v(i,j+1) * cos((lat(j) + lat(jp1)) * 0.5d0)
          vs = v(i,j) * cos((lat(jm1) + lat(j)) * 0.5d0)
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
          dhdt(i,j) = -( (fe - fw)/dlon + (fn - fs)/dlat ) / (radius * cos(lat(j)))
       end do
    end do

    if (present(no_momentum_tendency)) then
       if (no_momentum_tendency) then
          dudt = 0.d0
          dvdt = 0.d0
          return
       end if
    end if

    ! zonal momentum
    do j=1,nlat
       jp1 = mod(j,nlat)+1
       jm1 = mod(j-2+nlat,nlat)+1
       fcor = 2.d0*omega*sin(lat(j))
       do i=1,nlon
         ip1 = mod(i,nlon)+1
         im1 = mod(i-2+nlon,nlon)+1
          h_e = h(i,j)
          h_w = h(im1,j)
          v_avg = 0.25d0*(v(im1,j) + v(im1,j+1) + v(i,j) + v(i,j+1))
          dudt(i,j) = - u(i,j) * (u(ip1,j) - u(im1,j)) / (2.0d0 * dlon * radius * cos(lat(j))) &
                      - v_avg * (u(i,jp1) - u(i,jm1)) / (2.0d0 * dlat * radius) &
                      + u(i,j) * v_avg * tan(lat(j)) / radius &
                      - g * (h_e - h_w) / (dlon * radius * cos(lat(j))) &
                      + fcor * v_avg
       end do
    end do

    ! meridional momentum
    do j=2,nlat
        jp1 = j+1
        jm1 = j-1
        fcor = 2.d0*omega*sin(-pi/2.d0 + (j-1)*dlat)
        do i=1,nlon
          ip1 = mod(i,nlon)+1
          im1 = mod(i-2+nlon,nlon)+1
          h_n = h(i,j)
          h_s = h(i,jm1)
          u_avg = 0.25d0*(u(i,jm1) + u(ip1,jm1) + u(i,j) + u(ip1,j))
          dvdt(i,j) = - u_avg * (v(ip1,j) - v(im1,j)) / (2.0d0 * dlon * radius * cos(lat(j))) &
                      - v(i,j) * (v(i,jp1) - v(i,jm1)) / (2.0d0 * dlat * radius) &
                      - u_avg**2 * tan(lat(j)) / radius &
                      - g * (h_n - h_s) / (dlat * radius) &
                      - fcor * u_avg
       end do
    end do
    dvdt(:,1) = 0.d0
    dvdt(:,nlat+1) = 0.d0
  end subroutine rhs

end module equations_module
