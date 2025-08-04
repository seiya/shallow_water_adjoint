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

  real(dp) :: t, maxerr, l1err, l2err, alpha, mse, mass_res
  real(dp) :: t_ad, maxerr_ad, l1err_ad, l2err_ad, mse_ad, mass_res_ad
  integer :: n
  logical :: snapshot_flag

  call init_variables()
  call read_alpha(alpha)
  call read_snapshot_flag(snapshot_flag)
  call write_grid_params()
  call init_variables_fwd_ad()
  call read_alpha(alpha)
  call init_height_fwd_ad(h, h_ad, lon, lat)
  h_ad(nlon/2,nlat/2) = 1.0_dp
  call velocity_field_fwd_ad(u, u_ad, v, v_ad, lon, lat, alpha)
  do n = 0, nsteps
     t = n*dt
     if (snapshot_flag .and. mod(n,output_interval) == 0) then
        call write_snapshot(n, h_ad, u, v)
     end if
     if (n == nsteps) exit
     call rk4_step_fwd_ad(h, h_ad, hn, hn_ad, u, u_ad, v, v_ad, lat)
     h_ad = hn_ad
     h = hn
  end do
  call calc_mse_fwd_ad(h, h_ad, ha, mse, mse_ad)
  call calc_mass_residual_fwd_ad(h, h_ad, mass_res, mass_res_ad)
  call write_cost_log(mse, mass_res)
  print *, mse_ad, mass_res_ad
  call finalize_variables_fwd_ad()
end program shallow_water_test1_forward
