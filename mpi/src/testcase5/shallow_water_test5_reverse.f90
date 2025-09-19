program shallow_water_test5_reverse
  use mpi
  use constants_module, only: dp
  use cost_module, only: calc_mass_residual, calc_energy_residual, calc_wave_pattern
  use cost_module_ad, only: calc_mass_residual_rev_ad, calc_energy_residual_rev_ad, calc_mass_residual_fwd_rev_ad, calc_energy_residual_fwd_rev_ad
  use variables_module
  use variables_module_ad
  use equations_module
  use equations_module_ad
  use rk4_module
  use rk4_module_ad
  use io_module
  use io_module_ad
  use fautodiff_stack
  use mpi_decomp_module, only: init_decomp, finalize_decomp, mpi_rank, jend
  implicit none

  real(dp) :: mass_res, energy_res, wave
  real(dp) :: mass_res_ad, energy_res_ad
  real(dp) :: grad_dot_d
  integer :: n
  character(len=256) :: carg
  real(dp), allocatable :: un(:,:), vn(:,:)
  real(dp), allocatable :: un_ad(:,:), vn_ad(:,:)
  real(dp), allocatable :: hgeo(:,:), d(:,:)

  real(dp) :: h_ad_sum, h_ad_min, h_ad_max
  integer :: ierr

  call init_decomp(nx, ny)
  call init_variables()

  allocate(un(is:ie,js:je), vn(is:ie,js:jend+1))
  allocate(un_ad(is:ie,js:je), vn_ad(is:ie,js:jend+1))
  allocate(hgeo(is:ie,js:je))
  allocate(d(is:ie,js:je))

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
     call read_field(d, trim(carg))
  else
     d = 0.d0
  end if

  call geostrophic_velocity(u, v, hgeo)

  h = hgeo - b
  mass_res = calc_mass_residual(h)
  energy_res = calc_energy_residual(h, u, v)
  do n = 0, nsteps
     call fautodiff_stack_push_r(h)
     call fautodiff_stack_push_r(u)
     call fautodiff_stack_push_r(v)
     if (n == nsteps) exit
     call rk4_step(h, u, v, hn, un, vn)
     h = hn
     u = un
     v = vn
  end do
  wave = calc_wave_pattern(h)
  call calc_mass_residual_fwd_rev_ad()
  mass_res = calc_mass_residual(h)
  call calc_energy_residual_fwd_rev_ad()
  energy_res = calc_energy_residual(h, u, v)

  call finalize_variables_rev_ad()

  if (mpi_rank == 0) then
    energy_res_ad = 1.d0
  else
    energy_res_ad = 0.d0
  end if
  mass_res_ad = 0.d0
  h_ad = 0.d0
  u_ad = 0.d0
  v_ad = 0.d0

  call calc_energy_residual_rev_ad(h, h_ad, u, u_ad, v, v_ad, energy_res_ad)
  call calc_mass_residual_rev_ad(h, h_ad, mass_res_ad)

  do n = nsteps, 0, -1
     call fautodiff_stack_pop_r(v)
     call fautodiff_stack_pop_r(u)
     call fautodiff_stack_pop_r(h)
     if (n /= nsteps) then
        hn_ad = h_ad
        un_ad = u_ad
        vn_ad = v_ad
        h_ad = 0.d0
        u_ad = 0.d0
        v_ad = 0.d0
        call rk4_step_rev_ad(h, h_ad, u, u_ad, v, v_ad, hn_ad, un_ad, vn_ad)
     end if
     if (output_interval > 0) then
        if (mod(n, output_interval) == 0) then
           call exchange_halo_rev_ad(h_ad)
           call exchange_halo_rev_ad(u_ad)
           call exchange_halo_rev_ad(v_ad)
           call write_snapshot(n, h_ad, u_ad, v_ad)
        end if
     end if
  end do
  call geostrophic_velocity_rev_ad(u_ad, v_ad, h_ad)
  call exchange_halo_rev_ad(h_ad)
  call exchange_halo_rev_ad(u_ad)
  call exchange_halo_rev_ad(v_ad)
  if (output_interval == 0) then
     call write_snapshot(0, h_ad, u_ad, v_ad)
  end if
  grad_dot_d = sum(h_ad(istart:iend,jstart:jend)*d(istart:iend,jstart:jend))
  h_ad_sum = sum(h_ad(istart:iend,jstart:jend))
  h_ad_min = minval(h_ad(istart:iend,jstart:jend))
  h_ad_max = maxval(h_ad(istart:iend,jstart:jend))
  if (mpi_rank == 0) then
    call MPI_Reduce(MPI_IN_PLACE, h_ad_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(MPI_IN_PLACE, h_ad_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(MPI_IN_PLACE, h_ad_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(MPI_IN_PLACE, grad_dot_d, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    print *, h_ad_sum, h_ad_min, h_ad_max
    print *, grad_dot_d
  else
    call MPI_Reduce(h_ad_sum, 0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(h_ad_min, 0, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(h_ad_max, 0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(grad_dot_d, 0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  end if
  call init_variables_rev_ad()
  call finalize_variables()
  call finalize_decomp()
end program shallow_water_test5_reverse
