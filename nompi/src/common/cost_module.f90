module cost_module
  use constants_module, only: dp
  use variables_module, only: g, nx, ny, dx, dy
  implicit none
  real(dp), save :: reference_mass = -1.d0
  real(dp), save :: reference_energy = -1.d0
contains

  !> Compute sum of squared errors between numerical and analytic heights
  !$FAD CONSTANT_VARS: height_ana
  function calc_mse(height_num, height_ana) result(mse)
    real(dp), intent(in) :: height_num(nx,ny)
    real(dp), intent(in) :: height_ana(nx,ny)
    real(dp) :: mse

    mse = sum((height_num(:,:) - height_ana(:,:))**2)
    mse = mse / (nx * ny)
  end function calc_mse

  !> Compute deviation from the initial total mass
  function calc_mass_residual(height) result(residual)
    real(dp), intent(in) :: height(nx,ny)
    real(dp) :: residual, current_mass

    current_mass = sum(height(:,:))
    if (reference_mass < 0.d0) then
       reference_mass = current_mass
       residual = 0.d0
    else
       residual = current_mass - reference_mass
    end if
  end function calc_mass_residual

  !> Compute deviation from the initial total energy
  function calc_energy_residual(height, u, v) result(residual)
    real(dp), intent(in) :: height(nx,ny), u(nx,ny), v(nx,ny+1)
    real(dp) :: residual, current_energy
    real(dp) :: um, dudx, vm, dvdy
    integer :: ip1
    integer :: i, j

    current_energy = 0.0d0
    !$omp parallel do reduction(+:current_energy) private(ip1, um, dudx, vm, dvdy)
    do j = 1, ny
      do i = 1, nx
        ip1 = modulo(i,nx) + 1
        um = (u(i,j) + u(ip1,j)) * 0.5d0
        dudx = (u(ip1,j) - u(i,j)) / dx
        vm = (v(i,j) + v(i,j+1)) * 0.5d0
        dvdy = (v(i,j+1) - v(i,j)) / dy
        current_energy = current_energy &
                       + 0.5d0 * g * height(i,j)**2 &
                       + 0.5d0 * height(i,j) * ((um**2 + dudx**2 * dx**2 / 12.0d0) + (vm**2 + dvdy**2 * dy**2 / 12.0d0))
      end do
    end do
    if (reference_energy < 0.d0) then
       reference_energy = current_energy
       residual = 0.d0
    else
       residual = current_energy - reference_energy
    end if
  end function calc_energy_residual

  !> Measure wave pattern energy as variance from zonal mean
  !$FAD SKIP
  function calc_wave_pattern(height) result(pattern)
    real(dp), intent(in) :: height(nx,ny)
    real(dp) :: pattern
    real(dp) :: zonal_mean(ny)
    integer :: cnt
    integer :: i, j

    !$omp parallel workshare
    zonal_mean(:) = sum(height(:,:),dim=1)
    !$omp end parallel workshare
    !$omp parallel workshare
    zonal_mean(:) = zonal_mean(:) / nx
    !$omp end parallel workshare
    pattern = 0.d0
    !$omp parallel do reduction(+:pattern)
    do j = 1, ny
       do i = 1, nx
          pattern = pattern + (height(i,j) - zonal_mean(j))**2
       end do
    end do
    !$omp end parallel do
    pattern = pattern / (nx*ny)
  end function calc_wave_pattern

  !> Compute L1, L2, and Linf error norms
  !$FAD CONSTANT_VARS: height_ana
  subroutine calc_error_norms(height_num, height_ana, l1err, l2err, linf)
    real(dp), intent(in) :: height_num(nx,ny), height_ana(nx,ny)
    real(dp), intent(out) :: l1err, l2err, linf
    integer :: i, j
    real(dp) :: err

    linf = 0.d0
    l1err = 0.d0
    l2err = 0.d0
    !$omp parallel do private(err) reduction(max:linf) reduction(+:l1err,l2err)
    do j = 1, ny
       do i = 1, nx
          err = height_num(i,j) - height_ana(i,j)
          linf = max(linf, abs(err))
          l1err = l1err + abs(err)
          l2err = l2err + err*err
       end do
    end do
    !$omp end parallel do

    l1err = l1err / (nx * ny)
    l2err = sqrt(l2err / (nx * ny))
  end subroutine calc_error_norms

  !> Evaluate inner product of gradients with perturbation directions
  !$FAD SKIP
  function evaluate_gradient(grad_ic, grad_param, dir_ic, dir_param) result(ip)
    real(dp), intent(in) :: grad_ic(:,:), grad_param(:)
    real(dp), intent(in) :: dir_ic(size(grad_ic,1),size(grad_ic,2))
    real(dp), intent(in) :: dir_param(size(grad_param))
    real(dp) :: ip
    ip = sum(grad_ic*dir_ic) + sum(grad_param*dir_param)
  end function evaluate_gradient

end module cost_module
