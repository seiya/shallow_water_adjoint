program shallow_water_test1
  use constants_module, only: dp
  use cost_module, only: calc_mse, calc_mass_residual, calc_error_norms
  use variables_module
  use equations_module
  use rk4_module
  use io_module
  implicit none

  real(dp) :: t, maxerr, l1err, l2err, alpha, mse, mass_res
  integer :: n
  real(dp) :: un(nlon,nlat), vn(nlon,nlat+1)
  character(len=256) :: carg

  call init_variables()
  call read_alpha(alpha)
  call read_output_interval(output_interval)
  call write_grid_params()
  if (command_argument_count() >= 3) then
     call get_command_argument(3, carg)
     call read_field(h, trim(carg))
  else
     call init_height(h, lon, lat)
  end if
  mass_res = calc_mass_residual(h)
  call velocity_field(u, v, lon, lat, alpha)
  call open_error_file()
  do n = 0, nsteps
     t = n*dt
     call analytic_height(ha, lon, lat, t, alpha)
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
     call rk4_step(h, u, v, hn, un, vn, lat, no_momentum_tendency=.true.)
     h = hn
  end do
  call close_error_file()
  mse = calc_mse(h, ha)
  mass_res = calc_mass_residual(h)
  call write_cost_log(mse, mass_res)
  call finalize_variables()
end program shallow_water_test1
