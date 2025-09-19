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
  use mpi_decomp_module, only: init_decomp, finalize_decomp, mpi_rank, jend
  implicit none

  real(dp) :: t, mass_res, energy_res, wave
  real(dp) :: mass_res_ad, energy_res_ad
  integer :: n, i1, i2, j1, j2
  character(len=256) :: carg
  real(dp), allocatable :: un(:,:), vn(:,:)
  real(dp), allocatable :: un_ad(:,:), vn_ad(:,:)
  real(dp), allocatable :: hgeo(:,:)

  call init_decomp(nx, ny)
  call init_variables()

  allocate(un(is:ie,js:je), vn(is:ie,js:jend+1))
  allocate(un_ad(is:ie,js:je), vn_ad(is:ie,js:jend+1))
  allocate(hgeo(is:ie,js:je))

  call init_variables_fwd_ad()
  call read_output_interval(output_interval)
  call write_grid_params()
  call init_topography(b, x, y)
  if (command_argument_count() >= 2) then
     call get_command_argument(2, carg)
     call read_field(hgeo, trim(carg))
  else
     call init_geostrophic_height(hgeo, y)
  end if
  if (command_argument_count() >= 3) then
     call get_command_argument(3, carg)
     call read_field(h_ad, trim(carg))
  else
     h_ad = 0.0_dp
     i1 = max(nx/4-1, is)
     i2 = min(nx/4+2, ie)
     j1 = max(ny*3/4-1, js)
     j2 = min(ny*3/4+2, je)
     if (i1 <= i2 .and. j1 <= j2) h_ad(i1:i2, j1:j2) = 0.5_dp
     i1 = max(nx/4, is)
     i2 = min(nx/4+1, ie)
     j1 = max(ny*3/4, js)
     j2 = min(ny*3/4+1, je)
     if (i1 <= i2 .and. j1 <= j2) h_ad(i1:i2, j1:j2) = 1.0_dp
  end if
  call geostrophic_velocity_fwd_ad(u, u_ad, v, v_ad, hgeo, h_ad)
  h = hgeo - b
  mass_res = calc_mass_residual(h)
  energy_res = calc_energy_residual(h, u, v)
  do n = 0, nsteps
     t = n*dt
     if (output_interval > 0) then
        if (mod(n, output_interval) == 0) then
           call write_snapshot(n, h_ad, u_ad, v_ad)
        end if
     end if
     if (n == nsteps) exit
     call rk4_step_fwd_ad(h, h_ad, u, u_ad, v, v_ad, hn, hn_ad, un, un_ad, vn, vn_ad)
     h = hn; h_ad = hn_ad
     u = un; u_ad = un_ad
     v = vn; v_ad = vn_ad
  end do
  if (output_interval == 0) then
     call write_snapshot(nsteps, h_ad, u_ad, v_ad)
  end if
  mass_res = calc_mass_residual(h)
  call calc_mass_residual_fwd_ad(h, h_ad, mass_res, mass_res_ad)
  call calc_energy_residual_fwd_ad(h, h_ad, u, u_ad, v, v_ad, energy_res, energy_res_ad)
  wave = calc_wave_pattern(h)
  call write_cost_log2(mass_res, energy_res, wave)
  if (mpi_rank == 0) then
     print *, energy_res_ad, mass_res_ad, 0.d0
  end if
  call finalize_variables_fwd_ad()
  call finalize_decomp()
end program shallow_water_test5_forward
