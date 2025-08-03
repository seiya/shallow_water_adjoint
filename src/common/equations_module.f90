module equations_module
  use constants_module, only: dp
  use variables_module, only: nlon, nlat, pi, radius, dlon, dlat, h0, h1, omega
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
    real(dp), intent(out) :: u(nlon+1,nlat), v(nlon,nlat+1)
    real(dp), intent(in) :: lon(nlon), lat(nlat), alpha
    real(dp) :: u0, lon_edge
    integer :: i,j
    u0 = omega*radius
    do j=1,nlat
       do i=1,nlon+1
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

  !$FAD CONSTANT_VARS: lon, lat, t, alpha
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

  !$FAD CONSTANT_VARS: lat
  subroutine rhs(h,dhdt,u,v,lat)
    real(dp), intent(in) :: h(nlon,nlat), u(nlon+1,nlat), v(nlon,nlat+1), lat(nlat)
    real(dp), intent(out) :: dhdt(nlon,nlat)
    integer :: i,j,ip1,im1,jp1,jm1
    real(dp) :: fe,fw,fn,fs,ue,uw,vn,vs
    do j=1,nlat
       jp1 = min(j+1,nlat)
       jm1 = max(j-1,1)
       do i=1,nlon
          ip1 = mod(i,nlon)+1
          im1 = mod(i-2+nlon,nlon)+1
          ue = u(i+1,j)
          uw = u(i,j)
          if (ue > 0.d0) then
             fe = ue*h(i,j)
          else
             fe = ue*h(ip1,j)
          end if
          if (uw > 0.d0) then
             fw = uw*h(im1,j)
          else
             fw = uw*h(i,j)
          end if
          vn = v(i,j+1)
          vs = v(i,j)
          if (vn > 0.d0) then
             fn = vn*h(i,j)
          else
             fn = vn*h(i,jp1)
          end if
          if (vs > 0.d0) then
             fs = vs*h(i,jm1)
          else
             fs = vs*h(i,j)
          end if
          dhdt(i,j) = -((fe - fw)/(dlon*radius*cos(lat(j))) + (fn - fs)/(dlat*radius))
       end do
    end do
  end subroutine rhs

end module equations_module
