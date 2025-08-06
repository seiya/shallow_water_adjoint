program shallow_water_test2_forward
  use constants_module, only: dp
  use cost_module, only: calc_mse, calc_mass_residual
  use cost_module_ad, only: calc_mse_fwd_ad, calc_mass_residual_fwd_ad
  use variables_module
  use variables_module_ad
  use equations_module
  use equations_module_ad
  use rk4_module
  use rk4_module_ad
  use io_module
  implicit none

  real(dp) :: t, mse, mass_res
  real(dp) :: t_ad, mse_ad, mass_res_ad
  integer :: n
  character(len=256) :: carg
  real(dp) :: un(nlon,nlat), vn(nlon,nlat+1)
  real(dp) :: un_ad(nlon,nlat), vn_ad(nlon,nlat+1)

  call init_variables()
  call read_output_interval(output_interval)
  call write_grid_params()
  call init_variables_fwd_ad()
  ha = 0.0_dp
  if (command_argument_count() >= 3) then
     call get_command_argument(3, carg)
     call read_field(h, trim(carg))
     ha = h
  else
     call init_geostrophic_height_fwd_ad(h, h_ad, lon, lat)
     ha = h
  end if
  if (command_argument_count() >= 4) then
     call get_command_argument(4, carg)
     call read_field(h_ad, trim(carg))
  else
     h_ad(nlon/2, nlat/2) = 1.0_dp
  end if
  call geostrophic_velocity_fwd_ad(u, u_ad, v, v_ad, lat)
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
     call rk4_step_fwd_ad(h, h_ad, u, u_ad, v, v_ad, hn, hn_ad, un, un_ad, vn, vn_ad, lat)
     !print *, minval(hn), maxval(hn), minval(un), maxval(un), minval(vn), maxval(vn)
     print *, minval(hn_ad), maxval(hn_ad), minval(un_ad), maxval(un_ad), minval(vn_ad), maxval(vn_ad)
     h = hn
     u = un
     v = vn
     h_ad = hn_ad
     u_ad = un_ad
     v_ad = vn_ad
  end do
  mse = calc_mse(h, ha)
  mass_res = calc_mass_residual(h)
  call calc_mse_fwd_ad(h, h_ad, ha, mse, mse_ad)
  call calc_mass_residual_fwd_ad(h, h_ad, mass_res, mass_res_ad)
  call write_cost_log(mse, mass_res)
  print *, mse_ad, mass_res_ad
  call finalize_variables_fwd_ad()

contains

  subroutine init_geostrophic_height_fwd_ad(h, h_ad, lon, lat)
    real(dp), intent(out) :: h(nlon,nlat)
    real(dp), intent(out) :: h_ad(nlon,nlat)
    real(dp), intent(in)  :: lon(nlon)
    real(dp), intent(in)  :: lat(nlat)
    real(dp), parameter :: u0 = 20.d0
    real(dp) :: coeff
    integer :: i, j
    coeff = radius*omega*u0/g
    do j = 1, nlat
       do i = 1, nlon
          h_ad(i,j) = 0.0_dp
          h(i,j) = h0 - coeff * sin(lat(j))**2
       end do
    end do
  end subroutine init_geostrophic_height_fwd_ad

  subroutine geostrophic_velocity_fwd_ad(u, u_ad, v, v_ad, lat)
    real(dp), intent(out) :: u(nlon,nlat)
    real(dp), intent(out) :: u_ad(nlon,nlat)
    real(dp), intent(out) :: v(nlon,nlat+1)
    real(dp), intent(out) :: v_ad(nlon,nlat+1)
    real(dp), intent(in)  :: lat(nlat)
    real(dp), parameter :: u0 = 20.d0
    integer :: i, j
    do j = 1, nlat
       do i = 1, nlon
          u_ad(i,j) = 0.0_dp
          u(i,j) = u0 * cos(lat(j))
       end do
    end do
    v_ad = 0.0_dp
    v = 0.d0
  end subroutine geostrophic_velocity_fwd_ad

end program shallow_water_test2_forward
