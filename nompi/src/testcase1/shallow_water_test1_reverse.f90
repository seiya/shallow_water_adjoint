program shallow_water_test1_reverse
  use constants_module, only: dp
  use cost_module, only: calc_mse, calc_mass_residual
  use cost_module_ad, only: calc_mse_rev_ad, calc_mass_residual_rev_ad, calc_mass_residual_fwd_rev_ad
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

  real(dp) :: t, maxerr, l1err, l2err, mse, mass_res
  real(dp) :: maxerr_ad, l1err_ad, l2err_ad, mse_ad, mass_res_ad
  real(dp) :: grad_dot_d
  integer :: n
  character(len=256) :: carg
  real(dp), allocatable :: un(:,:), vn(:,:)
  real(dp), allocatable :: un_ad(:,:), vn_ad(:,:)
  real(dp), allocatable :: d(:,:)

  real(dp) :: h_ad_sum, h_ad_min, h_ad_max
  integer :: ierr

  call init_variables()

  allocate(un(nx,ny), vn(nx,ny+1))
  allocate(un_ad(nx,ny), vn_ad(nx,ny+1))
  allocate(d(nx,ny))

  call read_output_interval(output_interval)
  call write_grid_params()
  if (command_argument_count() >= 2) then
     call get_command_argument(2, carg)
     call read_field(h, trim(carg))
  else
     call init_height(h, x, y)
  end if
  if (command_argument_count() >= 3) then
     call get_command_argument(3, carg)
     call read_field(d, trim(carg))
  else
     d = 0.0_dp
  end if

  call velocity_field(u, v, x, y)
  do n = 0, nsteps
     call fautodiff_stack_push_r(h)
     if (n == nsteps) exit
     call rk4_step(h, u, v, hn, un, vn, no_momentum_tendency=.true.)
     h = hn
  end do
  t = nsteps*dt
  call analytic_height(ha, x, y, t)
  mse = calc_mse(h, ha)
  call calc_mass_residual_fwd_rev_ad()
  mass_res = calc_mass_residual(h)

  call finalize_variables_rev_ad()

  mse_ad = 1.0_dp
  mass_res_ad = 0.0_dp
  h_ad = 0.0_dp
  u_ad = 0.0_dp
  v_ad = 0.0_dp

  call calc_mass_residual_rev_ad(h_ad, mass_res_ad)
  call calc_mse_rev_ad(h, h_ad, ha, mse_ad)
  do n = nsteps, 0, -1
     call fautodiff_stack_pop_r(h)
     if (n .ne. nsteps) then
        hn_ad = h_ad
        h_ad = 0.0_dp
        call rk4_step_rev_ad(h, h_ad, u, u_ad, v, v_ad, hn_ad, un_ad, vn_ad, no_momentum_tendency=.true.)
     end if
     if (output_interval > 0) then
        if (mod(n, output_interval) == 0) then
           u_ad = 0.0
           v_ad = 0.0
           call write_snapshot(n, h_ad, u_ad, v_ad)
        end if
     end if
  end do
  if (output_interval == 0) then
     u_ad = 0.0
     v_ad = 0.0
     call write_snapshot(0, h_ad, u_ad, v_ad)
  end if
  h_ad_sum = sum(h_ad(:,:))
  h_ad_min = minval(h_ad(:,:))
  h_ad_max = maxval(h_ad(:,:))
  grad_dot_d = sum(h_ad(:,:)*d(:,:))
  print *, h_ad_sum, h_ad_min, h_ad_max
  print *, grad_dot_d
  call init_variables_rev_ad()

  call finalize_variables()

end program shallow_water_test1_reverse
