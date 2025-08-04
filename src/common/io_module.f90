module io_module
  use constants_module, only: dp, sp
  use variables_module, only: nlon, nlat, day, pi, h, u, v, hsp, usp, vsp
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
  subroutine read_snapshot_flag(flag)
    logical, intent(out) :: flag
    character(len=32) :: carg
    integer :: inum
    if (command_argument_count() >= 2) then
       call get_command_argument(2,carg)
       read(carg,*) inum
       flag = (inum /= 0)
    else
       flag = .true.
    end if
  end subroutine read_snapshot_flag

  !$FAD SKIP
  subroutine read_field(field, filename)
    !! Read a two-dimensional field from a binary file.
    real(dp), intent(out) :: field(nlon,nlat)
    character(len=*), intent(in) :: filename
    open(unit=50,file=filename,form='unformatted',access='stream',status='old')
    read(50) field
    close(50)
  end subroutine read_field

  !$FAD SKIP
  subroutine write_grid_params()
    open(unit=30,file='grid_params.txt',status='replace')
    write(30,*) nlon, nlat
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
    real(dp), intent(in) :: h(nlon,nlat)
    real(dp), intent(in) :: u(nlon+1,nlat), v(nlon,nlat+1)
    character(len=32) :: filename
    hsp = real(h,sp)
    usp = real(0.5d0*(u(1:nlon,:) + u(2:nlon+1,:)), sp)
    vsp = real(0.5d0*(v(:,1:nlat) + v(:,2:nlat+1)), sp)
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
end module io_module
