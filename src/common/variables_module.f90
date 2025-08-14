module variables_module
  use constants_module, only: dp, sp
  use mpi
  implicit none
  integer, parameter :: nx=128, ny=64
  integer, parameter :: ihalo=1
  integer, parameter :: is=1-ihalo
  integer, parameter :: ie=nx+ihalo
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

  !$FAD CONSTANT_VARS: x, y, b
  !$FAD CONSTANT_VARS: ha

contains

  subroutine init_variables()
    integer :: i, j
    allocate(x(is:ie), y(ny))
    allocate(h(is:ie,ny), hn(is:ie,ny), ha(is:ie,ny))
    allocate(u(is:ie,ny), v(is:ie,ny+1))
    allocate(b(is:ie,ny))

    !$omp parallel do
    do i=1,nx
        x(i) = (i-0.5d0)*dx
    end do
    !$omp end parallel do
    call exchange_halo_x_1d(x)
    !$omp parallel do
    do j=1,ny
       y(j) = (j-0.5d0)*dy
    end do
    !$omp end parallel do
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
  end subroutine finalize_variables

  subroutine exchange_halo_x(field)
    use mpi_decomp_module, only : comm_cart, nbr_west, nbr_east
    use mpi
    real(dp), intent(inout) :: field(is:,:)
    integer :: ierr, ny_local, cnt
    integer :: requests(4)

    ny_local = size(field,2)
    cnt = ny_local

    call MPI_Irecv(field(is,:),   cnt, MPI_DOUBLE_PRECISION, nbr_west, 0, comm_cart, requests(1), ierr)
    call MPI_Irecv(field(nx+1,:), cnt, MPI_DOUBLE_PRECISION, nbr_east, 1, comm_cart, requests(2), ierr)
    call MPI_Isend(field(nx,:),   cnt, MPI_DOUBLE_PRECISION, nbr_east, 0, comm_cart, requests(3), ierr)
    call MPI_Isend(field(1,:),    cnt, MPI_DOUBLE_PRECISION, nbr_west, 1, comm_cart, requests(4), ierr)
    call MPI_Waitall(4, requests, MPI_STATUSES_IGNORE, ierr)
  end subroutine exchange_halo_x

  subroutine exchange_halo_x_1d(field)
    use mpi_decomp_module, only : comm_cart, nbr_west, nbr_east
    use mpi
    real(dp), intent(inout) :: field(is:)
    integer :: ierr, cnt
    integer :: requests(4)

    cnt = ihalo

    call MPI_Irecv(field(is:0),   cnt, MPI_DOUBLE_PRECISION, nbr_west, 0, comm_cart, requests(1), ierr)
    call MPI_Irecv(field(nx+1:ie), cnt, MPI_DOUBLE_PRECISION, nbr_east, 1, comm_cart, requests(2), ierr)
    call MPI_Isend(field(nx+1-ihalo:nx), cnt, MPI_DOUBLE_PRECISION, nbr_east, 0, comm_cart, requests(3), ierr)
    call MPI_Isend(field(1:ihalo),       cnt, MPI_DOUBLE_PRECISION, nbr_west, 1, comm_cart, requests(4), ierr)
    call MPI_Waitall(4, requests, MPI_STATUSES_IGNORE, ierr)
  end subroutine exchange_halo_x_1d

end module variables_module
