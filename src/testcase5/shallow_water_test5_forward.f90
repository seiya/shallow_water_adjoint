program shallow_water_test5_forward
  use constants_module, only: dp
  use cost_module, only: calc_wave_pattern, calc_mass_residual, calc_energy_residual
  use cost_module_ad, only: calc_mass_residual_fwd_ad, calc_energy_residual_fwd_ad
  use variables_module
  use variables_module_ad
  use equations_module
  use equations_module_ad
  use rk4_module
  use rk4_module_ad
  use io_module
  implicit none

  real(dp) :: t, mass_res, energy_res, wave
  real(dp) :: mass_res_ad, energy_res_ad
  integer :: n
  character(len=256) :: carg
  real(dp) :: un(nx,ny), vn(nx,ny+1)
  real(dp) :: un_ad(nx,ny), vn_ad(nx,ny+1)

  call init_variables()
  call init_variables_fwd_ad()
  call read_output_interval(output_interval)
  call write_grid_params()
  call init_topography(b, x, y)
  if (command_argument_count() >= 2) then
     call get_command_argument(2, carg)
     call read_field(h, trim(carg))
  else
     h = h0 - b
  end if
  if (command_argument_count() >= 3) then
     call get_command_argument(3, carg)
     call read_field(h_ad, trim(carg))
  else
     h_ad = 0.d0
     h_ad(nx/2, ny/2) = 1.d0
  end if
  u_ad = 0.d0
  v_ad = 0.d0
  call velocity_field(u, v, x, y)
  mass_res = calc_mass_residual(h)
  energy_res = calc_energy_residual(h, u, v)
  do n = 0, nsteps
     t = n*dt
     if (output_interval /= -1) then
        if (output_interval == 0) then
           if (n == nsteps) call write_snapshot(n, h_ad, u, v)
        else if (mod(n, output_interval) == 0) then
           call write_snapshot(n, h_ad, u, v)
        end if
     end if
     if (n == nsteps) exit
     call rk4_step_fwd_ad(h, h_ad, u, u_ad, v, v_ad, hn, hn_ad, un, un_ad, vn, vn_ad)
     h = hn; h_ad = hn_ad
     u = un; u_ad = un_ad
     v = vn; v_ad = vn_ad
  end do
  mass_res = calc_mass_residual(h)
  call calc_mass_residual_fwd_ad(h, h_ad, mass_res, mass_res_ad)
  call calc_energy_residual_fwd_ad(h, h_ad, u, u_ad, v, v_ad, energy_res, energy_res_ad)
  wave = calc_wave_pattern(h)
  call write_cost_log2(mass_res, energy_res, wave)
  print *, energy_res_ad, mass_res_ad, 0.d0
  call finalize_variables_fwd_ad()
end program shallow_water_test5_forward
