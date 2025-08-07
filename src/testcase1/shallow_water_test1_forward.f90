program shallow_water_test1_forward
  use constants_module, only: dp
  use cost_module_ad, only: calc_mse_fwd_ad, calc_mass_residual_fwd_ad, calc_error_norms_fwd_ad
  use variables_module
  use variables_module_ad
  use equations_module
  use equations_module_ad
  use rk4_module
  use rk4_module_ad
  use io_module
  implicit none

  real(dp) :: t, maxerr, l1err, l2err, mse, mass_res
  real(dp) :: t_ad, maxerr_ad, l1err_ad, l2err_ad, mse_ad, mass_res_ad
  integer :: n
  character(len=256) :: carg
  real(dp) :: un(is:ie,ny), vn(is:ie,ny+1)
  real(dp) :: un_ad(is:ie,ny), vn_ad(is:ie,ny+1)

  call init_variables()
  call read_output_interval(output_interval)
  call write_grid_params()
  call init_variables_fwd_ad()
  if (command_argument_count() >= 2) then
     call get_command_argument(2, carg)
     call read_field(h, trim(carg))
  else
     call init_height(h, x, y)
  end if
  if (command_argument_count() >= 3) then
     call get_command_argument(3, carg)
     call read_field(h_ad, trim(carg))
  else
     h_ad(:,:) = 0.0_dp
     h_ad(nx/2-1:nx/2+2,ny/2-1:ny/2+2) = 0.5_dp
     h_ad(nx/2:nx/2+1,ny/2:ny/2+1) = 1.0_dp
  end if
  u_ad = 0.0_dp
  v_ad = 0.0_dp
  call velocity_field(u, v, x, y)
  do n = 0, nsteps
     t = n*dt
     if (output_interval > 0) then
        if (mod(n, output_interval) == 0) then
           call write_snapshot(n, h_ad, u, v)
        end if
     end if
     if (n == nsteps) exit
     call rk4_step_fwd_ad(h, h_ad, u, u_ad, v, v_ad, hn, hn_ad, un, un_ad, vn, vn_ad, no_momentum_tendency=.true.)
     h_ad = hn_ad
     h = hn
  end do
  if (output_interval == 0) then
     call write_snapshot(nsteps, h_ad, u, v)
  end if
  t = nsteps*dt
  call analytic_height(ha, x, y, t)
  call calc_mse_fwd_ad(h, h_ad, ha, mse, mse_ad)
  call calc_mass_residual_fwd_ad(h, h_ad, mass_res, mass_res_ad)
  call write_cost_log(mse, mass_res)
  print *, mse_ad, mass_res_ad
  call finalize_variables_fwd_ad()

end program shallow_water_test1_forward
