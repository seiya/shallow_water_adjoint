module variables_module
  use constants_module, only: dp, sp
  implicit none
  integer, parameter :: nlon=128, nlat=64
  real(dp), parameter :: pi=3.14159265358979323846d0
  real(dp), parameter :: radius=6371220.d0
  real(dp) :: g=9.80616d0
  real(dp), parameter :: day=86400.d0
  real(dp), parameter :: dlon=2.d0*pi/nlon, dlat=pi/nlat
  real(dp), parameter :: h0=10000.d0, h1=2000.d0
  real(dp) :: omega=2.d0*pi/(12.d0*day)
  real(dp), parameter :: dt=600.d0
  integer, parameter :: nsteps=nint(12.d0*day/dt)
  integer :: output_interval=48
  real(dp), allocatable :: lon(:), lat(:)
  real(dp), allocatable :: h(:,:), hn(:,:), ha(:,:)
  real(dp), allocatable :: u(:,:), v(:,:)
  real(sp), allocatable :: hsp(:,:), usp(:,:), vsp(:,:)

  !$FAD CONSTANT_VARS: lon, lat
  !$FAD CONSTANT_VARS: ha
  !$FAD CONSTANT_VARS: hsp, usp, vsp
  !$FAD CONSTANT_VARS: g, omega

contains

  subroutine init_variables()
    integer :: i, j
    allocate(lon(nlon), lat(nlat))
    allocate(h(nlon,nlat), hn(nlon,nlat), ha(nlon,nlat))
    allocate(u(nlon,nlat), v(nlon,nlat+1))
    allocate(hsp(nlon,nlat), usp(nlon,nlat), vsp(nlon,nlat))
    do i=1,nlon
       lon(i) = (i-0.5d0)*dlon
    end do
    do j=1,nlat
       lat(j) = -pi/2.d0 + (j-0.5d0)*dlat
    end do
  end subroutine init_variables

  subroutine finalize_variables()
    if (allocated(lon)) deallocate(lon)
    if (allocated(lat)) deallocate(lat)
    if (allocated(h)) deallocate(h)
    if (allocated(hn)) deallocate(hn)
    if (allocated(ha)) deallocate(ha)
    if (allocated(u)) deallocate(u)
    if (allocated(v)) deallocate(v)
    if (allocated(hsp)) deallocate(hsp)
    if (allocated(usp)) deallocate(usp)
    if (allocated(vsp)) deallocate(vsp)
  end subroutine finalize_variables

end module variables_module
