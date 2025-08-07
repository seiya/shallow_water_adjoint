program shallow_water_test2
  use constants_module, only: dp
  use cost_module, only: calc_mse, calc_mass_residual, calc_error_norms
  use variables_module
  use rk4_module
  use io_module
  implicit none

  real(dp) :: t, maxerr, l1err, l2err, mse, mass_res
  integer :: n
  character(len=256) :: carg
  real(dp) :: un(nx,ny), vn(nx,ny+1)

  call init_variables()
  call read_output_interval(output_interval)
  call write_grid_params()

  if (command_argument_count() >= 2) then
     call get_command_argument(2, carg)
     call read_field(h, trim(carg))
     ha = h
  else
     call init_geostrophic_height(h, y)
     ha = h
  end if

  call geostrophic_velocity(u, v, h)
  mass_res = calc_mass_residual(h)
  call open_error_file()
  do n = 0, nsteps
     t = n*dt
     call calc_error_norms(h, ha, l1err, l2err, maxerr)
     call write_error(t, l1err, l2err, maxerr)
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
  call close_error_file()
  mse = calc_mse(h, ha)
  mass_res = calc_mass_residual(h)
  call write_cost_log(mse, mass_res)
  call finalize_variables()

contains

  subroutine init_geostrophic_height(h, y)
    real(dp), intent(out) :: h(nx,ny)
    real(dp), intent(in) :: y(ny)
    integer :: i, j
    real(dp), parameter :: coeff = f0 * u0 * radius / g
    do j = 1, ny
       do i = 1, nx
          h(i,j) = h0 - coeff * sin(y(j)/radius)
       end do
    end do
  end subroutine init_geostrophic_height

  subroutine geostrophic_velocity(u, v, h)
    real(dp), intent(out) :: u(nx,ny), v(nx,ny+1)
    real(dp), intent(in) :: h(nx,ny)
    integer :: i, j
    integer :: ip1, im1, jp1, jm1
    do j = 1, ny
       jp1 = min(j+1, ny)
       jm1 = max(j-1, 1)
       do i = 1, nx
          im1 = mod(i-2+nx, nx) + 1
          u(i,j) = - g / f0 * ((h(im1,jp1) + h(i,jp1)) - (h(im1,jm1) + h(i,jm1))) / (4.0d0 * dy)
       end do
    end do
    do j = 2, ny
       jm1 = j - 1
       do i = 1, nx
          ip1 = mod(i, nx) + 1
          im1 = mod(i-2+nx, nx) + 1
          v(i,j) = g / f0 * ((h(ip1,jm1) + h(ip1,j)) - (h(im1,jm1) + h(im1,j))) / (4.0d0 * dx)
       end do
    end do
    v(:,1) = 0.0d0
    v(:,ny+1) = 0.0d0
  end subroutine geostrophic_velocity

end program shallow_water_test2
