program shallow_water_test5
  use constants_module, only: dp
  use cost_module, only: calc_mass_residual, calc_energy_residual, calc_wave_pattern
  use variables_module
  use equations_module
  use rk4_module
  use io_module
  implicit none

  real(dp) :: t, mass_res, energy_res, wave
  integer :: n
  real(dp) :: un(nx,ny), vn(nx,ny+1)
  character(len=256) :: carg
  real(dp) :: hgeo(nx,ny)

  call init_variables()
  call read_output_interval(output_interval)
  call write_grid_params()
  call init_topography(b, x, y)
  if (command_argument_count() >= 2) then
     call get_command_argument(2, carg)
     call read_field(hgeo, trim(carg))
  else
     call init_geostrophic_height(hgeo, y)
  end if
  call geostrophic_velocity(u, v, hgeo)
  h = hgeo - b
  mass_res = calc_mass_residual(h)
  energy_res = calc_energy_residual(h, u, v)
  do n = 0, nsteps
     t = n*dt
     if (output_interval /= -1) then
        if (output_interval == 0) then
           if (n == nsteps) call write_snapshot(n, h, u, v)
        else if (mod(n, output_interval) == 0) then
           call write_snapshot(n, h, u, v)
        end if
     end if
     if (n == nsteps) exit
     call rk4_step(h, u, v, hn, un, vn)
     h = hn
     u = un
     v = vn
  end do
  mass_res = calc_mass_residual(h)
  energy_res = calc_energy_residual(h, u, v)
  wave = calc_wave_pattern(h)
  call write_cost_log2(mass_res, energy_res, wave)
  call finalize_variables()
end program shallow_water_test5
