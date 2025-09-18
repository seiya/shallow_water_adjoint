program shallow_water_test1
  use constants_module, only: dp
  use cost_module, only: calc_mse, calc_mass_residual, calc_error_norms
  use variables_module
  use equations_module
  use rk4_module
  use io_module
  use mpi_decomp_module, only: init_decomp, finalize_decomp
  implicit none

  real(dp) :: t, maxerr, l1err, l2err, mse, mass_res
  integer :: n
  real(dp), allocatable :: un(:,:), vn(:,:)
  character(len=256) :: carg

  call init_decomp(nx, ny)
  call init_variables()

  allocate(un(is:ie,js:je), vn(is:ie,js:jend+1))

  call read_output_interval(output_interval)
  call write_grid_params()
  if (command_argument_count() >= 2) then
     call get_command_argument(2, carg)
     call read_field(h, trim(carg))
  else
     call init_height(h, x, y)
  end if
  mass_res = calc_mass_residual(h)
  call velocity_field(u, v, x, y)
  call open_error_file()
  do n = 0, nsteps
     t = n*dt
     call analytic_height(ha, x, y, t)
     call calc_error_norms(h, ha, l1err, l2err, maxerr)
     call write_error(t, l1err, l2err, maxerr)
     if (output_interval > 0) then
        if (mod(n, output_interval) == 0) then
           call write_snapshot(n, h, u, v)
        end if
     end if
     if (n == nsteps) exit
     call rk4_step(h, u, v, hn, un, vn, no_momentum_tendency=.true.)
     h = hn
  end do
  if (output_interval == 0) then
     call write_snapshot(nsteps, h, u, v)
  end if
  call close_error_file()
  mse = calc_mse(h, ha)
  mass_res = calc_mass_residual(h)
  call write_cost_log(mse, mass_res)
  call finalize_variables()
  call finalize_decomp()
end program shallow_water_test1
