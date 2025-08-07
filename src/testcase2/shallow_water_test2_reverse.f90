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
  real(dp) :: un(nx,ny), vn(nx,ny+1)
  real(dp) :: un_ad(nx,ny), vn_ad(nx,ny+1)

  call init_variables()
  call read_output_interval(output_interval)
  call write_grid_params()
  if (command_argument_count() >= 2) then
     call get_command_argument(2, carg)
     call read_field(h, trim(carg))
  else
     call init_geostrophic_height(h, y)
  end if
  ha = h
  allocate(d(nx,ny))
  if (command_argument_count() >= 3) then
     call get_command_argument(3, carg)
     call read_field(d, trim(carg))
  else
     d = 0.0_dp
  end if
  call geostrophic_velocity(u, v, h)
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
        call rk4_step_rev_ad(h, h_ad, u, u_ad, v, v_ad, hn_ad, un_ad, vn_ad)
     end if
     if (output_interval /= -1) then
        if (output_interval == 0) then
           if (n == 0) then
              call geostrophic_velocity_rev_ad(u, u_ad, v, v_ad, h, h_ad)
              call write_snapshot(n, h_ad, u, v)
           end if
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

  subroutine geostrophic_velocity_rev_ad(u, u_ad, v, v_ad, h, h_ad)
    real(dp), intent(in)    :: u(nx,ny), v(nx,ny+1)
    real(dp), intent(inout) :: u_ad(nx,ny), v_ad(nx,ny+1)
    real(dp), intent(in)    :: h(nx,ny)
    real(dp), intent(inout) :: h_ad(nx,ny)
    integer :: i, j, ip1, im1, jp1, jm1
    real(dp) :: coeff

    coeff = - g / f0 / (4.0d0 * dy)
    do j = 1, ny
       jp1 = min(j+1, ny)
       jm1 = max(j-1, 1)
       do i = 1, nx
          im1 = mod(i-2+nx, nx) + 1
          h_ad(im1, jp1) = h_ad(im1, jp1) + coeff * u_ad(i,j)
          h_ad(i,   jp1) = h_ad(i,   jp1) + coeff * u_ad(i,j)
          h_ad(im1, jm1) = h_ad(im1, jm1) - coeff * u_ad(i,j)
          h_ad(i,   jm1) = h_ad(i,   jm1) - coeff * u_ad(i,j)
       end do
    end do
    u_ad(:,:) = 0.0d0

    v_ad(:,1) = 0.0d0
    v_ad(:,ny+1) = 0.0d0
    coeff = g / f0 / (4.0d0 * dx)
    do j = 2, ny
       jm1 = j - 1
       do i = 1, nx
          ip1 = mod(i, nx) + 1
          im1 = mod(i-2+nx, nx) + 1
          h_ad(ip1, jm1) = h_ad(ip1, jm1) + coeff * v_ad(i,j)
          h_ad(ip1, j)   = h_ad(ip1, j)   + coeff * v_ad(i,j)
          h_ad(im1, jm1) = h_ad(im1, jm1) - coeff * v_ad(i,j)
          h_ad(im1, j)   = h_ad(im1, j)   - coeff * v_ad(i,j)
       end do
    end do
    v_ad(:,:) = 0.0d0
  end subroutine geostrophic_velocity_rev_ad

end program shallow_water_test2_reverse
