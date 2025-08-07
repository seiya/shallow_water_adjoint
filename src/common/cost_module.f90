module cost_module
  use constants_module, only: dp
  use variables_module, only: g
  implicit none
  real(dp), save :: reference_mass = -1.d0
  real(dp), save :: reference_energy = -1.d0
contains

  !> Compute sum of squared errors between numerical and analytic heights
  !$FAD CONSTANT_VARS: height_ana
  function calc_mse(height_num, height_ana) result(mse)
    real(dp), intent(in) :: height_num(:,:), height_ana(:,:)
    real(dp) :: mse
    mse = sum((height_num - height_ana)**2) / (size(height_ana))
  end function calc_mse

  !> Compute deviation from the initial total mass
  function calc_mass_residual(height) result(residual)
    real(dp), intent(in) :: height(:,:)
    real(dp) :: residual, current_mass
    current_mass = sum(height)
    if (reference_mass < 0.d0) then
       reference_mass = current_mass
       residual = 0.d0
    else
       residual = current_mass - reference_mass
    end if
  end function calc_mass_residual

  !> Compute deviation from the initial total energy
  function calc_energy_residual(height, u, v) result(residual)
    real(dp), intent(in) :: height(:,:), u(:,:), v(:,:)
    real(dp) :: residual, current_energy
    current_energy = sum(0.5d0*g*height**2 + 0.5d0*height*(u**2 + v(:,1:size(height,2))**2))
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
    real(dp), intent(in) :: height(:,:)
    real(dp) :: pattern
    real(dp), allocatable :: zonal_mean(:)
    integer :: nx, ny, i, j
    nx = size(height,1)
    ny = size(height,2)
    allocate(zonal_mean(ny))
    zonal_mean = sum(height,dim=1)/nx
    pattern = 0.d0
    do j = 1, ny
       do i = 1, nx
          pattern = pattern + (height(i,j) - zonal_mean(j))**2
       end do
    end do
    pattern = pattern / (nx*ny)
    deallocate(zonal_mean)
  end function calc_wave_pattern

  !> Compute L1, L2, and Linf error norms
  !$FAD CONSTANT_VARS: height_ana
  subroutine calc_error_norms(height_num, height_ana, l1err, l2err, maxerr)
    real(dp), intent(in) :: height_num(:,:), height_ana(:,:)
    real(dp), intent(out) :: l1err, l2err, maxerr
    integer :: i, j, nx, ny
    real(dp) :: err

    nx = size(height_num, 1)
    ny = size(height_num, 2)
    maxerr = 0.d0
    l1err = 0.d0
    l2err = 0.d0

    do j = 1, ny
       do i = 1, nx
          err = height_num(i,j) - height_ana(i,j)
          maxerr = max(maxerr, abs(err))
          l1err = l1err + abs(err)
          l2err = l2err + err*err
       end do
    end do

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
