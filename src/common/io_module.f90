module io_module
  use constants_module, only: dp, sp
  use variables_module, only: nx, ny, day, pi, h, u, v, &
                              ihalo, is, ie, js, je, exchange_halo
  implicit none
contains

  !$FAD SKIP
  subroutine read_output_interval(interval)
    integer, intent(inout) :: interval
    character(len=32) :: carg
    if (command_argument_count() >= 1) then
       call get_command_argument(1,carg)
       read(carg,*) interval
    end if
  end subroutine read_output_interval

  !$FAD SKIP
  subroutine read_field(field, filename)
    !! Read a two-dimensional field from a binary file.
    use mpi
    use mpi_decomp_module, only: istart, iend, jstart, jend
    real(dp), intent(out) :: field(is:ie,js:je)
    character(len=*), intent(in) :: filename
    integer :: fh, ierr
    integer :: ni, nj, j
    integer(kind=MPI_Offset_kind) :: offset
    real(dp), allocatable :: buf(:)

    ni = iend - istart + 1
    nj = jend - jstart + 1

    allocate(buf(ni))
    call MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    if (ierr /= MPI_SUCCESS) stop 'MPI_File_open failed'

    do j = 0, nj-1
       offset = int(((jstart - 1 + j) * nx + (istart - 1)), MPI_OFFSET_KIND) * 8_MPI_OFFSET_KIND
       call MPI_File_read_at(fh, offset, buf, ni, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
       if (ierr /= MPI_SUCCESS) stop 'MPI_File_read_at failed'
       field(istart:iend, jstart + j) = buf
    end do

    call MPI_File_close(fh, ierr)
    deallocate(buf)
    call exchange_halo(field)
  end subroutine read_field

  !$FAD SKIP
  subroutine write_grid_params()
    open(unit=30,file='grid_params.txt',status='replace')
    write(30,*) nx, ny
    close(30)
  end subroutine write_grid_params

  !$FAD SKIP
  subroutine open_error_file()
    open(unit=10,file='error.dat',status='replace')
  end subroutine open_error_file

  !$FAD SKIP
  subroutine write_error(t,l1err,l2err,maxerr)
    real(dp), intent(in) :: t, l1err, l2err, maxerr
    write(10,'(f10.4,3(1x,e14.6))') t/day, l1err, l2err, maxerr
  end subroutine write_error

  !$FAD SKIP
  subroutine close_error_file()
    close(10)
  end subroutine close_error_file

  !$FAD SKIP
  subroutine write_snapshot(n, h, u, v)
    use mpi
    use mpi_decomp_module, only: istart, iend, jstart, jend
    integer, intent(in) :: n
    real(dp), intent(in) :: h(is:ie,js:je)
    real(dp), intent(in) :: u(is:ie,js:je), v(is:ie,js:jend+1)
    real(sp) :: hsp(istart:iend,jstart:jend)
    real(sp) :: usp(istart:iend,jstart:jend)
    real(sp) :: vsp(istart:iend,jstart:jend)
    integer :: ni, nj
    character(len=32) :: filename
    integer :: fh, newtype, ierr
    integer(kind=MPI_Offset_kind) :: disp
    integer :: i, j
    do j = jstart, jend
      do i = istart, iend
        hsp(i,j) = real(h(i,j), sp)
        usp(i,j) = real(0.5d0 * (u(i,j) + u(i+1,j)), sp)
        vsp(i,j) = real(0.5d0 * (v(i,j) + v(i,j+1)), sp)
      end do
    end do
    ni = iend - istart + 1
    nj = jend - jstart + 1
    write(filename,'("snapshot_",i4.4,".bin")') n
    call MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)
    if (ierr /= MPI_SUCCESS) stop 'MPI_File_open failed'
    call MPI_File_set_size(fh, 0_MPI_OFFSET_KIND, ierr)
    call MPI_Type_create_subarray(2, [nx, ny], [ni, nj], [istart-1, jstart-1], MPI_ORDER_FORTRAN, MPI_REAL, newtype, ierr)
    call MPI_Type_commit(newtype, ierr)

    disp = 0_MPI_OFFSET_KIND
    call MPI_File_set_view(fh, disp, MPI_REAL, newtype, 'native', MPI_INFO_NULL, ierr)
    call MPI_File_write_all(fh, hsp, ni * nj, MPI_REAL, MPI_STATUS_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) stop 'MPI_File_write_all (hsp) failed'

    disp = 4_MPI_OFFSET_KIND * (nx * ny) ! 4 byte * nx * ny
    call MPI_File_set_view(fh, disp, MPI_REAL, newtype, 'native', MPI_INFO_NULL, ierr)
    call MPI_File_write_all(fh, usp, ni * nj, MPI_REAL, MPI_STATUS_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) stop 'MPI_File_write_all (usp) failed'

    disp = 8_MPI_OFFSET_KIND * (nx * ny) ! 2 * 4 byte * nx * ny
    call MPI_File_set_view(fh, disp, MPI_REAL, newtype, 'native', MPI_INFO_NULL, ierr)
    call MPI_File_write_all(fh, vsp, ni * nj, MPI_REAL, MPI_STATUS_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) stop 'MPI_File_write_all (vsp) failed'

    call MPI_Type_free(newtype, ierr)
    call MPI_File_close(fh, ierr)
    !open(unit=20,file=filename,form='unformatted',access='stream',status='replace')
    !write(20) hsp, usp, vsp
    !close(20)
  end subroutine write_snapshot

  !$FAD SKIP
  subroutine write_cost_log(mse,mass_res)
    real(dp), intent(in) :: mse, mass_res
    open(unit=40,file='cost.log',status='replace')
    write(40,'(a,1x,e25.16)') 'MSE', mse
    write(40,'(a,1x,e25.16)') 'MassResidual', mass_res
    close(40)
  end subroutine write_cost_log

  !$FAD SKIP
  subroutine write_cost_log2(mass_res, energy_res, wave)
    real(dp), intent(in) :: mass_res, energy_res, wave
    open(unit=40,file='cost.log',status='replace')
    write(40,'(a,1x,e25.16)') 'MassResidual', mass_res
    write(40,'(a,1x,e25.16)') 'EnergyResidual', energy_res
    write(40,'(a,1x,e25.16)') 'WavePattern', wave
    close(40)
  end subroutine write_cost_log2
end module io_module
