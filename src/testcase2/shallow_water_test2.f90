program shallow_water_test2
  use constants_module, only: dp
  use cost_module, only: calc_mse, calc_mass_residual, calc_error_norms
  use variables_module
  use rk4_module
  use io_module
  implicit none

  real(dp) :: t, maxerr, l1err, l2err, mse, mass_res
  integer :: n
  character(len=256) :: carg
  real(dp) :: un(nlon,nlat), vn(nlon,nlat+1)

  call init_variables()
  call read_output_interval(output_interval)
  call write_grid_params()

  if (command_argument_count() >= 3) then
     call get_command_argument(3, carg)
     call read_field(h, trim(carg))
     ha = h
  else
     call init_geostrophic_height(h, lon, lat)
     ha = h
  end if

  call geostrophic_velocity(u, v, lat)
  mass_res = calc_mass_residual(h)
  call open_error_file()
  do n = 0, nsteps
     t = n*dt
     call calc_error_norms(h, ha, lat, l1err, l2err, maxerr)
     call write_error(t, l1err, l2err, maxerr)
     if (output_interval /= -1) then
        if (output_interval == 0) then
           if (n == nsteps) call write_snapshot(n, h, u, v)
        else if (mod(n, output_interval) == 0) then
           call write_snapshot(n, h, u, v)
        end if
     end if
     if (n == nsteps) exit
     call rk4_step(h, u, v, hn, un, vn, lat)
     h = hn
     u = un
     v = vn
  end do
  call close_error_file()
  mse = calc_mse(h, ha)
  mass_res = calc_mass_residual(h)
  call write_cost_log(mse, mass_res)
  call finalize_variables()

contains

  !$FAD CONSTANT_VARS: lon, lat
  subroutine init_geostrophic_height(h, lon, lat)
    real(dp), intent(out) :: h(nlon,nlat)
    real(dp), intent(in) :: lon(nlon), lat(nlat)
    integer :: i, j
    real(dp) :: coeff
    real(dp), parameter :: u0 = 20.d0
    coeff = radius*omega*u0/g
    do j = 1, nlat
       do i = 1, nlon
          h(i,j) = h0 - coeff * sin(lat(j))**2
       end do
    end do
  end subroutine init_geostrophic_height

  !$FAD CONSTANT_VARS: lat
  subroutine geostrophic_velocity(u, v, lat)
    real(dp), intent(out) :: u(nlon,nlat), v(nlon,nlat+1)
    real(dp), intent(in) :: lat(nlat)
    integer :: i, j
    real(dp), parameter :: u0 = 20.d0
    do j = 1, nlat
       do i = 1, nlon
          u(i,j) = u0 * cos(lat(j))
       end do
    end do
    v = 0.d0
  end subroutine geostrophic_velocity

end program shallow_water_test2
