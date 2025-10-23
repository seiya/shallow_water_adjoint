program shallow_water_test2_forward
  use constants_module_ad, only: dp
  use cost_module_ad, only: calc_mse, calc_mass_residual, calc_mse_fwd_ad, calc_mass_residual_fwd_ad
  use variables_module_ad
  use equations_module_ad
  use rk4_module_ad
  use io_module_ad
  implicit none

  real(dp) :: t, mse, mass_res
  real(dp) :: t_ad, mse_ad, mass_res_ad
  integer :: n, i1, i2, j1, j2
  character(len=256) :: carg
  real(dp), allocatable :: un(:,:), vn(:,:)
  real(dp), allocatable :: un_ad(:,:), vn_ad(:,:)

  call init_variables()

  allocate(un(nx,ny), vn(nx,ny+1))
  allocate(un_ad(nx,ny), vn_ad(nx,ny+1))

  call read_output_interval(output_interval)
  call write_grid_params()
  call init_variables_fwd_ad()
  ha = 0.0_dp
  if (command_argument_count() >= 2) then
     call get_command_argument(2, carg)
     call read_field(h, trim(carg))
     ha = h
  else
     call init_geostrophic_height_fwd_ad(h, h_ad, y)
     ha = h
  end if
  if (command_argument_count() >= 3) then
     call get_command_argument(3, carg)
     call read_field(h_ad, trim(carg))
  else
     h_ad = 0.0_dp
     i1 = nx/2-1
     i2 = nx/2+2
     j1 = ny/2-1
     j2 = ny/2+2
     h_ad(i1:i2, j1:j2) = 0.5_dp
     i1 = nx/2
     i2 = nx/2+1
     j1 = ny/2
     j2 = ny/2+1
     h_ad(i1:i2, j1:j2) = 1.0_dp
  end if
  call geostrophic_velocity_fwd_ad(u, u_ad, v, v_ad, h, h_ad)
  do n = 0, nsteps
     t = n*dt
     if (output_interval > 0) then
        if (mod(n, output_interval) == 0) then
           call write_snapshot(n, h_ad, u_ad, v_ad)
        end if
     end if
     if (n == nsteps) exit
     call rk4_step_fwd_ad(h, h_ad, u, u_ad, v, v_ad, hn, hn_ad, un, un_ad, vn, vn_ad)
     h = hn
     u = un
     v = vn
     h_ad = hn_ad
     u_ad = un_ad
     v_ad = vn_ad
  end do
  if (output_interval == 0) then
     call write_snapshot(nsteps, h_ad, u_ad, v_ad)
  end if
  mse = calc_mse(h, ha)
  mass_res = calc_mass_residual(h)
  call calc_mse_fwd_ad(h, h_ad, ha, mse, mse_ad)
  call calc_mass_residual_fwd_ad(h, h_ad, mass_res, mass_res_ad)
  call write_cost_log(mse, mass_res)
  print *, mse_ad, mass_res_ad
  call finalize_variables_fwd_ad()
end program shallow_water_test2_forward
