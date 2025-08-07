module variables_module
  use constants_module, only: dp, sp
  implicit none
  integer, parameter :: nlon=128, nlat=64
  real(dp), parameter :: pi=3.14159265358979323846d0
  real(dp), parameter :: g=9.80616d0
  real(dp), parameter :: day=86400.d0
  real(dp), parameter :: Lx=2.d0*pi, Ly=1.d0
  real(dp), parameter :: dx=Lx/nlon, dy=Ly/nlat
  real(dp), parameter :: h0=10000.d0, h1=2000.d0
  real(dp), parameter :: f0=1.0d-4
  real(dp), parameter :: dt=600.d0
  integer, parameter :: nsteps=nint(12.d0*day/dt)
  integer :: output_interval=48
  real(dp), allocatable :: x(:), y(:)
  real(dp), allocatable :: h(:,:), hn(:,:), ha(:,:)
  real(dp), allocatable :: u(:,:), v(:,:)
  real(sp), allocatable :: hsp(:,:), usp(:,:), vsp(:,:)

  !$FAD CONSTANT_VARS: x, y
  !$FAD CONSTANT_VARS: ha
  !$FAD CONSTANT_VARS: hsp, usp, vsp
  !$FAD CONSTANT_VARS: f0

contains

  subroutine init_variables()
    integer :: i, j
    allocate(x(nlon), y(nlat))
    allocate(h(nlon,nlat), hn(nlon,nlat), ha(nlon,nlat))
    allocate(u(nlon,nlat), v(nlon,nlat+1))
    allocate(hsp(nlon,nlat), usp(nlon,nlat), vsp(nlon,nlat))
    do i=1,nlon
       x(i) = (i-0.5d0)*dx
    end do
    do j=1,nlat
       y(j) = (j-0.5d0)*dy
    end do
  end subroutine init_variables

  subroutine finalize_variables()
    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)
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
