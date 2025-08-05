program shallow_water_test2_reverse
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

  real(dp) :: t, mse, mass_res
  real(dp) :: mse_ad, mass_res_ad, grad_dot_d
  integer :: n
  character(len=256) :: carg
  real(dp), allocatable :: d(:,:)
  real(dp) :: un(nlon,nlat), vn(nlon,nlat+1)
  real(dp) :: un_ad(nlon,nlat), vn_ad(nlon,nlat+1)

  call init_variables()
  call read_output_interval(output_interval)
  call write_grid_params()
  if (command_argument_count() >= 3) then
     call get_command_argument(3, carg)
     call read_field(h, trim(carg))
  else
     call init_geostrophic_height(h, lon, lat)
  end if
  ha = h
  allocate(d(nlon,nlat))
  if (command_argument_count() >= 4) then
     call get_command_argument(4, carg)
     call read_field(d, trim(carg))
  else
     d = 0.0_dp
  end if
  call geostrophic_velocity(u, v, lat)
  do n = 0, nsteps
     call fautodiff_stack_push_r(h)
     call fautodiff_stack_push_r(u)
     call fautodiff_stack_push_r(v)
     if (n == nsteps) exit
     call rk4_step(h, u, v, hn, un, vn, lat)
     h = hn
     u = un
     v = vn
  end do
  mse = calc_mse(h, ha)
  mass_res = calc_mass_residual(h)

  call finalize_variables_rev_ad()

  mse_ad = 1.0_dp
  mass_res_ad = 0.0_dp
  h_ad = 0.0_dp

  call calc_mass_residual_rev_ad(h, h_ad, mass_res_ad)
  call calc_mse_rev_ad(h, h_ad, ha, mse_ad)
  do n = nsteps, 0, -1
     call fautodiff_stack_pop_r(v)
     call fautodiff_stack_pop_r(u)
     call fautodiff_stack_pop_r(h)
     if (n /= nsteps) then
        hn_ad = h_ad
        un_ad = u_ad
        vn_ad = v_ad
        h_ad = 0.0_dp
        u_ad = 0.0_dp
        v_ad = 0.0_dp
        call rk4_step_rev_ad(h, h_ad, u, u_ad, v, v_ad, hn_ad, un_ad, vn_ad, lat)
     end if
     if (output_interval /= -1) then
        if (output_interval == 0) then
           if (n == 0) call write_snapshot(n, h_ad, u, v)
        else if (mod(n, output_interval) == 0) then
           call write_snapshot(n, h_ad, u, v)
        end if
     end if
  end do
  grad_dot_d = sum(h_ad*d)
  print *, sum(h_ad), minval(h_ad), maxval(h_ad)
  print *, grad_dot_d
  call init_variables_rev_ad()
  call finalize_variables()

contains

  subroutine init_geostrophic_height(h, lon, lat)
    real(dp), intent(out) :: h(nlon,nlat)
    real(dp), intent(in) :: lon(nlon), lat(nlat)
    real(dp), parameter :: u0 = 20.d0
    real(dp) :: coeff
    integer :: i, j
    coeff = radius*omega*u0/g
    do j = 1, nlat
       do i = 1, nlon
          h(i,j) = h0 - coeff * sin(lat(j))**2
       end do
    end do
  end subroutine init_geostrophic_height

  subroutine geostrophic_velocity(u, v, lat)
    real(dp), intent(out) :: u(nlon,nlat), v(nlon,nlat+1)
    real(dp), intent(in) :: lat(nlat)
    real(dp), parameter :: u0 = 20.d0
    integer :: i, j
    do j = 1, nlat
       do i = 1, nlon
          u(i,j) = u0 * cos(lat(j))
       end do
    end do
    v = 0.d0
  end subroutine geostrophic_velocity

end program shallow_water_test2_reverse
