module io_module
  use constants_module, only: dp, sp
  use variables_module, only: nx, ny, day, pi, h, u, v, &
                              ihalo, is, ie, exchange_halo_x
  implicit none
contains
  !$FAD SKIP
  subroutine read_alpha(alpha)
    real(dp), intent(out) :: alpha
    character(len=32) :: carg
    if (command_argument_count() >= 1) then
       call get_command_argument(1,carg)
       read(carg,*) alpha
    else
       alpha = 0.d0
    end if
    alpha = alpha*pi/180.d0
  end subroutine read_alpha

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
    real(dp), intent(out) :: field(is:ie,ny)
    character(len=*), intent(in) :: filename
    open(unit=50,file=filename,form='unformatted',access='stream',status='old')
    read(50) field(1:nx,1:ny)
    close(50)
    call exchange_halo_x(field)
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
    integer, intent(in) :: n
    real(dp), intent(in) :: h(is:ie,ny)
    real(dp), intent(in) :: u(is:ie,ny), v(is:ie,ny+1)
    real(sp) :: hsp(nx,ny), usp(nx,ny), vsp(nx,ny)
    character(len=32) :: filename
    hsp = real(h(1:nx,:),sp)
    usp = real(0.5d0*(u(1:nx,:) + u(2:nx+1,:)), sp)
    vsp = real(0.5d0*(v(1:nx,1:ny) + v(1:nx,2:ny+1)), sp)
    write(filename,'("snapshot_",i4.4,".bin")') n
    open(unit=20,file=filename,form='unformatted',access='stream',status='replace')
    write(20) hsp, usp, vsp
    close(20)
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
