module variables_module
  use constants_module, only: dp, sp
  use mpi
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
  integer :: is, ie, js, je

  !$FAD CONSTANT_VARS: x, y, b
  !$FAD CONSTANT_VARS: ha

contains

  subroutine init_variables()
    use mpi_decomp_module, only : &
      nbr_west, nbr_east, nbr_south, nbr_north, &
      istart, iend, jstart, jend
    integer :: i, j

    if (nbr_west == MPI_PROC_NULL) then
      is = istart
    else
      is = istart - ihalo
    end if
    if (nbr_east == MPI_PROC_NULL) then
      ie = iend
    else
      ie = iend + ihalo
    end if
    if (nbr_south == MPI_PROC_NULL) then
      js = jstart
    else
      js = jstart - ihalo
    end if
    if (nbr_north == MPI_PROC_NULL) then
      je = jend
    else
      je = jend + ihalo
    end if

    allocate(x(is:ie), y(js:je))
    allocate(h(is:ie,js:je), hn(is:ie,js:je), ha(is:ie,js:je))
    allocate(u(is:ie,js:je), v(is:ie,js:jend+1))
    allocate(b(is:ie,js:je))

    !$omp parallel do
    do i = is, ie
        x(i) = (i-0.5d0)*dx
    end do
    !$omp end parallel do
    call exchange_halo_x_1d(x)
    !$omp parallel do
    do j = js, je
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
    use mpi_decomp_module, only : &
      comm_cart, nbr_west, nbr_east, nbr_south, nbr_north, &
      istart, iend, jstart, jend
    use mpi
    real(dp), intent(inout) :: field(is:,:)
    integer :: ierr, cntx, cnty
    integer :: requests(8)
    real(dp) :: recvbuf(ihalo,jstart:jend,2)

    cntx = (iend - istart + 1) * ihalo
    cnty = (jend - jstart + 1) * ihalo

    call MPI_Irecv(recvbuf(:,:,1),                 cnty, MPI_DOUBLE_PRECISION, nbr_west,  0, comm_cart, requests(1), ierr)
    call MPI_Irecv(recvbuf(:,:,2),                 cnty, MPI_DOUBLE_PRECISION, nbr_east,  1, comm_cart, requests(2), ierr)
    call MPI_Irecv(field(istart:iend,js:jstart-1), cntx, MPI_DOUBLE_PRECISION, nbr_south, 2, comm_cart, requests(3), ierr)
    call MPI_Irecv(field(istart:iend,jend+1:je),   cntx, MPI_DOUBLE_PRECISION, nbr_north, 3, comm_cart, requests(4), ierr)
    call MPI_Isend(field(iend-ihalo+1:iend,jstart:jend),     cnty, MPI_DOUBLE_PRECISION, nbr_east,  0, comm_cart, requests(5), ierr)
    call MPI_Isend(field(istart:istart+ihalo-1,jstart:jend), cnty, MPI_DOUBLE_PRECISION, nbr_west,  1, comm_cart, requests(6), ierr)
    call MPI_Isend(field(istart:iend,jend-ihalo+1:jend),     cntx, MPI_DOUBLE_PRECISION, nbr_north, 2, comm_cart, requests(7), ierr)
    call MPI_Isend(field(istart:iend,jstart:jstart+ihalo-1), cntx, MPI_DOUBLE_PRECISION, nbr_south, 3, comm_cart, requests(8), ierr)
    call MPI_Waitall(8, requests, MPI_STATUSES_IGNORE, ierr)
    field(is:istart-1,jstart:jend) = recvbuf(:,:,1)
    field(iend+1:ie,jstart:jend) = recvbuf(:,:,2)
  end subroutine exchange_halo_x

  subroutine exchange_halo_x_1d(field)
    use mpi_decomp_module, only : &
      comm_cart, nbr_west, nbr_east, &
      istart, iend
    use mpi
    real(dp), intent(inout) :: field(is:)
    integer :: ierr, cnt
    integer :: requests(4)

    cnt = ihalo

    call MPI_Irecv(field(is:istart-1),   cnt, MPI_DOUBLE_PRECISION, nbr_west, 0, comm_cart, requests(1), ierr)
    call MPI_Irecv(field(iend+1:ie),     cnt, MPI_DOUBLE_PRECISION, nbr_east, 1, comm_cart, requests(2), ierr)
    call MPI_Isend(field(iend-ihalo+1:iend),     cnt, MPI_DOUBLE_PRECISION, nbr_east, 0, comm_cart, requests(3), ierr)
    call MPI_Isend(field(istart:istart+ihalo-1), cnt, MPI_DOUBLE_PRECISION, nbr_west, 1, comm_cart, requests(4), ierr)
    call MPI_Waitall(4, requests, MPI_STATUSES_IGNORE, ierr)
  end subroutine exchange_halo_x_1d

end module variables_module
