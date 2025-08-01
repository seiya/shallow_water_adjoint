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
  logical :: snapshot_flag
  call init_variables()
  call read_alpha(alpha)
  call read_snapshot_flag(snapshot_flag)
  call write_grid_params()
  call init_height(h, lon, lat)
  mass_res = calc_mass_residual(h)
  call velocity_field(u, v, lon, lat, alpha)
  call open_error_file()
  do n=0,nsteps
     t = n*dt
     call analytic_height(ha, lon, lat, t, alpha)
     call calc_error_norms(h, ha, lat, l1err, l2err, maxerr)
     call write_error(t, l1err, l2err, maxerr)
    if (snapshot_flag .and. mod(n,output_interval) == 0) then
       call write_snapshot(n)
    end if
     if (n == nsteps) exit
     call rk4_step(h, hn, u, v, lat)
     h = hn
  end do
  call close_error_file()
  mse = calc_mse(h, ha)
  mass_res = calc_mass_residual(h)
  call write_cost_log(mse, mass_res)
  call finalize_variables()
end program shallow_water_test1
