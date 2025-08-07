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
  real(dp) :: un(nx,ny), vn(nx,ny+1)
  real(dp) :: un_ad(nx,ny), vn_ad(nx,ny+1)

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
     call init_geostrophic_height_fwd_ad(h, h_ad, y)
     ha = h
  end if
  if (command_argument_count() >= 4) then
     call get_command_argument(4, carg)
     call read_field(h_ad, trim(carg))
  else
     h_ad(nx/2-1:nx/2+2, ny/2-1:ny/2+2) = 0.5_dp
     h_ad(nx/2:nx/2+1, ny/2:ny/2+1) = 1.0_dp
  end if
  call geostrophic_velocity_fwd_ad(u, u_ad, v, v_ad, h, h_ad)
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

  subroutine init_geostrophic_height_fwd_ad(h, h_ad, y)
    real(dp), intent(out) :: h(nx,ny)
    real(dp), intent(out) :: h_ad(nx,ny)
    real(dp), intent(in)  :: y(ny)
    integer :: i, j
    real(dp), parameter :: coeff = f0 * u0 * radius / g
   do j = 1, ny
       do i = 1, nx
          h_ad(i,j) = 0.0_dp
          h(i,j) = h0 - coeff * sin(y(j)/radius)
       end do
    end do
  end subroutine init_geostrophic_height_fwd_ad

    subroutine geostrophic_velocity_fwd_ad(u, u_ad, v, v_ad, h, h_ad)
    real(dp), intent(out) :: u(nx,ny)
    real(dp), intent(out) :: u_ad(nx,ny)
    real(dp), intent(out) :: v(nx,ny+1)
    real(dp), intent(out) :: v_ad(nx,ny+1)
    real(dp), intent(in) :: h(nx,ny)
    real(dp), intent(in) :: h_ad(nx,ny)
    integer :: i, j
    integer :: ip1, im1, jp1, jm1

    do j = 1, ny
       jp1 = min(j+1, ny)
       jm1 = max(j-1, 1)
       do i = 1, nx
          im1 = mod(i-2+nx, nx) + 1
          u_ad(i,j) = - g / f0 * ((h_ad(im1,jp1) + h_ad(i,jp1)) - (h_ad(im1,jm1) + h_ad(i,jm1))) / (4.0d0 * dy)
          u(i,j) = - g / f0 * ((h(im1,jp1) + h(i,jp1)) - (h(im1,jm1) + h(i,jm1))) / (4.0d0 * dy)
       end do
    end do
    do j = 2, ny
       jm1 = j - 1
       do i = 1, nx
          ip1 = mod(i, nx) + 1
          im1 = mod(i-2+nx, nx) + 1
          v_ad(i,j) = g / f0 * ((h_ad(ip1,jm1) + h_ad(ip1,j)) - (h_ad(im1,jm1) + h_ad(im1,j))) / (4.0d0 * dx)
          v(i,j) = g / f0 * ((h(ip1,jm1) + h(ip1,j)) - (h(im1,jm1) + h(im1,j))) / (4.0d0 * dx)
       end do
    end do
    v_ad(:,1) = 0.0d0
    v_ad(:,ny+1) = 0.0d0
    v(:,1) = 0.0d0
    v(:,ny+1) = 0.0d0
  end subroutine geostrophic_velocity_fwd_ad

end program shallow_water_test2_forward
