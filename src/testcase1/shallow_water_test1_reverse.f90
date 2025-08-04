program shallow_water_test1_reverse
  use constants_module, only: dp
  use cost_module, only: calc_mse, calc_mass_residual
  use cost_module_ad, only: calc_mse_rev_ad, calc_mass_residual_rev_ad
  use variables_module
  use variables_module_ad
  use equations_module
  use equations_module_ad
  use rk4_module
  use rk4_module_ad
  use io_module
  use io_module_ad
  use fautodiff_stack
  implicit none

  real(dp) :: t, maxerr, l1err, l2err, alpha, mse, mass_res
  real(dp) :: maxerr_ad, l1err_ad, l2err_ad, mse_ad, mass_res_ad
  real(dp) :: grad_dot_d
  integer :: n
  logical :: snapshot_flag
  character(len=256) :: carg
  real(dp), allocatable :: d(:,:)

  call init_variables()
  call read_alpha(alpha)
  call read_snapshot_flag(snapshot_flag)
  call write_grid_params()
  if (command_argument_count() >= 3) then
     call get_command_argument(3, carg)
     call read_field(h, trim(carg))
  else
     call init_height(h, lon, lat)
  end if
  allocate(d(nlon,nlat))
  if (command_argument_count() >= 4) then
     call get_command_argument(4, carg)
     call read_field(d, trim(carg))
  else
     d = 0.0_dp
  end if
  mass_res = calc_mass_residual(h)
  call velocity_field(u, v, lon, lat, alpha)
  do n = 0, nsteps
     call fautodiff_stack_push_r(h)
     if (n == nsteps) exit
     call rk4_step(h, hn, u, v, lat)
     h = hn
  end do
  t = nsteps*dt
  call analytic_height(ha, lon, lat, t, alpha)
  mse = calc_mse(h, ha)
  mass_res = calc_mass_residual(h)

  call finalize_variables_rev_ad()

  mse_ad = 1.0_dp
  mass_res_ad = 0.0_dp
  h_ad = 0.0_dp

  call calc_mass_residual_rev_ad(h, h_ad, mass_res_ad)
  call calc_mse_rev_ad(h, h_ad, ha, mse_ad)
  do n = nsteps, 0, -1
     call fautodiff_stack_pop_r(h)
     if (n .ne. nsteps) then
        hn_ad = h_ad
        h_ad = 0.0_dp
        call rk4_step_rev_ad(h, h_ad, hn_ad, u, u_ad, v, v_ad, lat)
     end if
     if (snapshot_flag .and. mod(n,output_interval) == 0) then
        call write_snapshot(n, h_ad, u, v)
     end if
  end do
  !call velocity_field_rev_ad(u, u_ad, v, v_ad, lon, lat, alpha)
  !call init_height_rev_ad(h, h_ad, lon, lat)
  grad_dot_d = sum(h_ad*d)
  print *, sum(h_ad), minval(h_ad), maxval(h_ad)
  print *, grad_dot_d
  call init_variables_rev_ad()

  call finalize_variables()

end program shallow_water_test1_reverse
