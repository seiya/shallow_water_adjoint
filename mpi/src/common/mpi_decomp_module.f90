module mpi_decomp_module
  use mpi
  implicit none
  integer :: comm_cart = MPI_COMM_NULL
  integer :: mpi_rank = -1, mpi_size = 0
  integer :: dims(2) = 0, coords(2) = 0
  logical :: periods(2) = (/ .false., .false. /)
  integer :: nbr_west, nbr_east, nbr_south, nbr_north
  integer :: istart, iend, jstart, jend
contains

  !$FAD SKIP
  subroutine init_decomp(nx_global, ny_global)
    integer, intent(in) :: nx_global, ny_global
    integer :: cart_rank
    integer :: ierr
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierr)
    dims = 0
    call MPI_Dims_create(mpi_size, 2, dims, ierr)
    periods = (/ .true., .false. /)
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .false., comm_cart, ierr)
    call MPI_Comm_rank(comm_cart, cart_rank, ierr)
    call MPI_Cart_coords(comm_cart, cart_rank, 2, coords, ierr)
    istart = block_start(coords(1), dims(1), nx_global)
    iend   = block_end(coords(1), dims(1), nx_global)
    jstart = block_start(coords(2), dims(2), ny_global)
    jend   = block_end(coords(2), dims(2), ny_global)
    call MPI_Cart_shift(comm_cart, 0, 1, nbr_west, nbr_east, ierr)
    call MPI_Cart_shift(comm_cart, 1, 1, nbr_south, nbr_north, ierr)
  end subroutine init_decomp

  !$FAD SKIP
  function block_start(coord, nprocs, n_global) result(start)
    integer, intent(in) :: coord, nprocs, n_global
    integer :: start, base, rem
    base = n_global / nprocs
    rem  = mod(n_global, nprocs)
    start = coord * base + min(coord, rem) + 1
  end function block_start

  !$FAD SKIP
  function block_end(coord, nprocs, n_global) result(endp)
    integer, intent(in) :: coord, nprocs, n_global
    integer :: endp, base, rem
    base = n_global / nprocs
    rem  = mod(n_global, nprocs)
    endp = block_start(coord, nprocs, n_global) + base - 1
    if (coord < rem) endp = endp + 1
  end function block_end

  !$FAD SKIP
  subroutine finalize_decomp()
    integer :: ierr
    if (comm_cart /= MPI_COMM_NULL) then
      call MPI_Comm_free(comm_cart, ierr)
      comm_cart = MPI_COMM_NULL
    end if
    call MPI_Finalize(ierr)
  end subroutine finalize_decomp

end module mpi_decomp_module
