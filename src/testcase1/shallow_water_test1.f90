program shallow_water_test1
  implicit none
  ! Shallow water equation solver for cosine bell advection test case
  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: nlon=64, nlat=32
  real(dp), parameter :: pi=3.14159265358979323846d0
  real(dp), parameter :: radius=6371220.d0, g=9.80616d0
  real(dp), parameter :: day=86400.d0
  real(dp), parameter :: dlon=2.d0*pi/nlon, dlat=pi/nlat
  real(dp), parameter :: h0=10000.d0, h1=2000.d0
  real(dp), parameter :: omega=2.d0*pi/(12.d0*day)
  real(dp), parameter :: dt=600.d0
  integer, parameter :: nsteps=nint(12.d0*day/dt)
  real(dp) :: lon(nlon), lat(nlat)
  real(dp) :: h(nlon,nlat), hn(nlon,nlat)
  real(dp) :: ha(nlon,nlat)
  real(dp) :: u(nlon,nlat), v(nlon,nlat)
  real(dp) :: t, maxerr, alpha
  integer :: i,j,n
  character(len=32) :: carg

  ! Read solid body rotation angle alpha in degrees from command line
  if (command_argument_count() >= 1) then
     call get_command_argument(1,carg)
     read(carg,*) alpha
  else
     alpha = 0.d0
  end if
  alpha = alpha*pi/180.d0

  do i=1,nlon
     lon(i) = (i-1)*dlon
  end do
  do j=1,nlat
     lat(j) = -pi/2.d0 + (j-0.5d0)*dlat
  end do

  call init_height(h, lon, lat)
  call velocity_field(u, v, lon, lat, alpha)

  open(unit=10,file='error.dat',status='replace')
  do n=0,nsteps
     t = n*dt
     call analytic_height(ha, lon, lat, t, alpha)
     maxerr = 0.d0
     do j=1,nlat
        do i=1,nlon
           maxerr = max(maxerr, abs(h(i,j)-ha(i,j)))
        end do
     end do
     write(10,'(f10.4,1x,e14.6)') t/day, maxerr
     if (n == nsteps) exit
     call step(h, hn, u, v, lat)
     h = hn
  end do
  close(10)

contains

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

  subroutine velocity_field(u,v,lon,lat,alpha)
    real(dp), intent(out) :: u(nlon,nlat), v(nlon,nlat)
    real(dp), intent(in) :: lon(nlon), lat(nlat), alpha
    real(dp) :: u0
    integer :: i,j
    u0 = omega*radius
    do j=1,nlat
       do i=1,nlon
          u(i,j) = u0*(cos(lat(j))*cos(alpha) + sin(lat(j))*cos(lon(i))*sin(alpha))
          v(i,j) = -u0*sin(lon(i))*sin(alpha)
       end do
    end do
  end subroutine velocity_field

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

  subroutine gc_distance(lon1,lat1,lon2,lat2,dist)
    real(dp), intent(in) :: lon1,lat1,lon2,lat2
    real(dp), intent(out) :: dist
    dist = radius*acos( sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon1-lon2) )
  end subroutine gc_distance

  subroutine rotate_point(lon,lat,alpha,theta,lonp,latp)
    real(dp), intent(in) :: lon, lat, alpha, theta
    real(dp), intent(out) :: lonp, latp
    real(dp) :: x,y,z,xp,yp,zp
    real(dp) :: n1,n2,n3,dot,c1,c2,c3,ct,st
    x = cos(lat)*cos(lon)
    y = cos(lat)*sin(lon)
    z = sin(lat)
    n1 = sin(alpha)
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

  subroutine rhs(h,dhdt,u,v,lat)
    real(dp), intent(in) :: h(nlon,nlat), u(nlon,nlat), v(nlon,nlat), lat(nlat)
    real(dp), intent(out) :: dhdt(nlon,nlat)
    integer :: i,j,ip,im,jp,jm
    real(dp) :: hx, hy
    do j=1,nlat
       jp = min(j+1,nlat)
       jm = max(j-1,1)
       do i=1,nlon
          ip = mod(i,nlon)+1
          im = mod(i-2+nlon,nlon)+1
          if (u(i,j) > 0.d0) then
             hx = (h(i,j)-h(im,j))/(dlon*radius*cos(lat(j)))
          else
             hx = (h(ip,j)-h(i,j))/(dlon*radius*cos(lat(j)))
          end if
          if (v(i,j) > 0.d0) then
             hy = (h(i,j)-h(i,jm))/(dlat*radius)
          else
             hy = (h(i,jp)-h(i,j))/(dlat*radius)
          end if
          dhdt(i,j) = -(u(i,j)*hx + v(i,j)*hy)
       end do
    end do
  end subroutine rhs

  subroutine step(h,hn,u,v,lat)
    real(dp), intent(in) :: h(nlon,nlat), u(nlon,nlat), v(nlon,nlat), lat(nlat)
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
  end subroutine step

end program shallow_water_test1
