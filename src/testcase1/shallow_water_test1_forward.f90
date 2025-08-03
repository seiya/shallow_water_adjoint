program shallow_water_test1_forward
  use constants_module, only: dp
  use cost_module, only: calc_mse, calc_mass_residual, calc_error_norms
  use cost_module_ad, only: calc_mse_fwd_ad, calc_mass_residual_fwd_ad, calc_error_norms_fwd_ad
  use variables_module
  use variables_module_ad
  use equations_module
  use equations_module_ad
  use rk4_module
  use rk4_module_ad
  use io_module
  use io_module_ad
  implicit none

  real(dp) :: t, maxerr, l1err, l2err, alpha, mse, mass_res
  real(dp) :: t_ad, maxerr_ad, l1err_ad, l2err_ad, mse_ad, mass_res_ad
  integer :: n

  call init_variables()
  call init_variables_fwd_ad()
  call read_alpha(alpha)
  call init_height_fwd_ad(h, h_ad, lon, lat)
  call calc_mass_residual_fwd_ad(h, h_ad, mass_res, mass_res_ad)
  call velocity_field_fwd_ad(u, u_ad, v, v_ad, lon, lat, alpha)
  call open_error_file()
  do n = 0, nsteps
     t = n*dt
     call analytic_height_fwd_ad(ha, ha_ad, lon, lat, t, alpha)
     call calc_error_norms_fwd_ad(h, h_ad, ha, lat, l1err, l1err_ad, l2err, l2err_ad, maxerr, maxerr_ad)
     call write_error(t, l1err, l2err, maxerr)
     if (n == nsteps) exit
     call rk4_step_fwd_ad(h, h_ad, hn, hn_ad, u, u_ad, v, v_ad, lat)
     h_ad = hn_ad
     h = hn
  end do
  call close_error_file()
  call calc_mse_fwd_ad(h, h_ad, ha, mse, mse_ad)
  call calc_mass_residual_fwd_ad(h, h_ad, mass_res, mass_res_ad)
  call write_cost_log(mse, mass_res)
  call finalize_variables_fwd_ad()
end program shallow_water_test1_forward
