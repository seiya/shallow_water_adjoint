module variables_module
  use constants_module, only: dp, sp
  implicit none
  integer, parameter :: nx=128, ny=64
  integer, parameter :: ihalo=1
  real(dp), parameter :: pi=3.14159265358979323846d0
  real(dp), parameter :: radius=6371220.d0
  real(dp), parameter :: g=9.80616d0
  real(dp), parameter :: day=86400.d0
  real(dp), parameter :: Lx=2.d0*pi*radius, Ly=pi*radius
  real(dp), parameter :: dx=Lx/nx, dy=Ly/ny
  real(dp), parameter :: h0=10000.d0, h1=2000.d0
  real(dp), parameter :: omega=2.0d0*pi/(12.d0*day)
  real(dp), parameter :: u0 = omega * radius
  real(dp), parameter :: f0=1.0d-4
  real(dp), parameter :: dt=600.d0
  integer, parameter :: nsteps=nint(12.d0*day/dt)
  integer :: output_interval=48
  real(dp), allocatable :: x(:), y(:)
  real(dp), allocatable :: h(:,:), hn(:,:), ha(:,:)
  real(dp), allocatable :: u(:,:), v(:,:)
  real(dp), allocatable :: b(:,:)
    real(sp), allocatable :: hsp(:,:), usp(:,:), vsp(:,:)

    !$FAD CONSTANT_VARS: x, y, b
    !$FAD CONSTANT_VARS: ha
    !$FAD CONSTANT_VARS: hsp, usp, vsp

  contains

  subroutine init_variables()
    integer :: i, j
    allocate(x(1-ihalo:nx+ihalo), y(ny))
    allocate(h(1-ihalo:nx+ihalo,ny), hn(1-ihalo:nx+ihalo,ny), &
             ha(1-ihalo:nx+ihalo,ny))
    allocate(u(1-ihalo:nx+ihalo,ny), v(1-ihalo:nx+ihalo,ny+1))
    allocate(b(1-ihalo:nx+ihalo,ny))
    allocate(hsp(1-ihalo:nx+ihalo,ny), usp(1-ihalo:nx+ihalo,ny), &
             vsp(1-ihalo:nx+ihalo,ny))
      do i=1,nx
         x(i) = (i-0.5d0)*dx
      end do
      call exchange_halo_x_1d(x)
    do j=1,ny
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
    if (allocated(b)) deallocate(b)
    if (allocated(hsp)) deallocate(hsp)
    if (allocated(usp)) deallocate(usp)
    if (allocated(vsp)) deallocate(vsp)
  end subroutine finalize_variables

    subroutine exchange_halo_x(field)
      real(dp), intent(inout) :: field(1-ihalo:,:)
      field(1-ihalo:0, :) = field(nx+1-ihalo:nx, :)
      field(nx+1:nx+ihalo, :) = field(1:ihalo, :)
    end subroutine exchange_halo_x

    subroutine exchange_halo_x_1d(field)
      real(dp), intent(inout) :: field(1-ihalo:)
      field(1-ihalo:0) = field(nx+1-ihalo:nx)
      field(nx+1:nx+ihalo) = field(1:ihalo)
    end subroutine exchange_halo_x_1d

end module variables_module
