program shallow_water_test1
  use cost_module, only: calc_mse, calc_mass_residual, dp, sp
  implicit none
  ! Shallow water equation solver for cosine bell advection test case
  integer, parameter :: nlon=128, nlat=64
  real(dp), parameter :: pi=3.14159265358979323846d0
  real(dp), parameter :: radius=6371220.d0, g=9.80616d0
  real(dp), parameter :: day=86400.d0
  real(dp), parameter :: dlon=2.d0*pi/nlon, dlat=pi/nlat
  real(dp), parameter :: h0=10000.d0, h1=2000.d0
  real(dp), parameter :: omega=2.d0*pi/(12.d0*day)
  real(dp), parameter :: dt=600.d0
  integer, parameter :: nsteps=nint(12.d0*day/dt)
  integer, parameter :: output_interval=48
  real(dp) :: lon(nlon), lat(nlat)
  real(dp) :: h(nlon,nlat), hn(nlon,nlat)
  real(dp) :: ha(nlon,nlat)
  ! Arakawa C-grid staggering: u on zonal cell edges, v on meridional edges
  real(dp) :: u(nlon+1,nlat), v(nlon,nlat+1)
  real(dp) :: t, maxerr, l1err, l2err, alpha, wtsum, err, w, mse, mass_res
  real(sp) :: hsp(nlon,nlat), usp(nlon,nlat), vsp(nlon,nlat)
  integer :: i,j,n
  character(len=32) :: carg, filename

  ! Read solid body rotation angle alpha in degrees from command line
  if (command_argument_count() >= 1) then
     call get_command_argument(1,carg)
     read(carg,*) alpha
  else
     alpha = 0.d0
  end if
  alpha = alpha*pi/180.d0

  ! Write grid dimensions for use by plotting script
  open(unit=30,file='grid_params.txt',status='replace')
  write(30,*) nlon, nlat
  close(30)

  do i=1,nlon
     lon(i) = (i-0.5d0)*dlon
  end do
  do j=1,nlat
     lat(j) = -pi/2.d0 + (j-0.5d0)*dlat
  end do

  call init_height(h, lon, lat)
  ! initialize reference mass for conservation check
  mass_res = calc_mass_residual(h)
  call velocity_field(u, v, lon, lat, alpha)

  open(unit=10,file='error.dat',status='replace')  ! t(days), L1, L2, Linf
  do n=0,nsteps
     t = n*dt
     call analytic_height(ha, lon, lat, t, alpha)
     maxerr = 0.d0
     l1err = 0.d0
     l2err = 0.d0
     wtsum = 0.d0
     do j=1,nlat
        w = cos(lat(j))
        wtsum = wtsum + w
        do i=1,nlon
           err = h(i,j) - ha(i,j)
           maxerr = max(maxerr, abs(err))
           l1err = l1err + abs(err)*w
           l2err = l2err + err*err*w
        end do
     end do
     l1err = l1err/(nlon*wtsum)
     l2err = sqrt(l2err/(nlon*wtsum))
     write(10,'(f10.4,3(1x,e14.6))') t/day, l1err, l2err, maxerr

     if (mod(n,output_interval) == 0) then
        hsp = real(h,sp)
        ! Output cell-centered velocities averaged from C-grid edges
        usp = real(0.5d0*(u(1:nlon,:) + u(2:nlon+1,:)), sp)
        vsp = real(0.5d0*(v(:,1:nlat) + v(:,2:nlat+1)), sp)
        write(filename,'("snapshot_",i4.4,".bin")') n
        open(unit=20,file=filename,form='unformatted',access='stream',status='replace')
        write(20) hsp, usp, vsp
        close(20)
     end if

     if (n == nsteps) exit
     call step(h, hn, u, v, lat)
     h = hn
  end do
  close(10)

  ! Compute cost function diagnostics at final time
  open(unit=40,file='cost.log',status='replace')
  mse = calc_mse(h, ha)
  mass_res = calc_mass_residual(h)
  write(40,'(a,1x,e16.8)') 'MSE', mse
  write(40,'(a,1x,e16.8)') 'MassResidual', mass_res
  close(40)

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
    ! Compute velocity components on an Arakawa C-grid
    real(dp), intent(out) :: u(nlon+1,nlat), v(nlon,nlat+1)
    real(dp), intent(in) :: lon(nlon), lat(nlat), alpha
    real(dp) :: u0, lon_edge
    integer :: i,j
    u0 = omega*radius
    ! Zonal velocity on longitudinal cell edges
    do j=1,nlat
       do i=1,nlon+1
          lon_edge = (i-1)*dlon
          u(i,j) = u0*(cos(lat(j))*cos(alpha) + sin(lat(j))*cos(lon_edge)*sin(alpha))
       end do
    end do
    ! Meridional velocity on latitudinal cell edges
    do j=1,nlat+1
       do i=1,nlon
          v(i,j) = -u0*sin(lon(i))*sin(alpha)
       end do
    end do
    v(:,1) = 0.d0
    v(:,nlat+1) = 0.d0
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

  subroutine rhs(h,dhdt,u,v,lat)
    real(dp), intent(in) :: h(nlon,nlat), u(nlon+1,nlat), v(nlon,nlat+1), lat(nlat)
    real(dp), intent(out) :: dhdt(nlon,nlat)
    ! First-order upwind fluxes on a C-grid
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

  subroutine step(h,hn,u,v,lat)
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
  end subroutine step

end program shallow_water_test1
