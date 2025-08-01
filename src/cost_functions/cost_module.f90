module cost_module
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: sp=kind(1.0)
  real(dp), save :: reference_mass = -1.d0
contains

  !> Compute sum of squared errors between numerical and analytic heights
  function calc_mse(height_num, height_ana) result(mse)
    real(dp), intent(in) :: height_num(:,:), height_ana(:,:)
    real(dp) :: mse
    mse = sum((height_num - height_ana)**2)
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

  !> Evaluate inner product of gradients with perturbation directions
  function evaluate_gradient(grad_ic, grad_param, dir_ic, dir_param) result(ip)
    real(dp), intent(in) :: grad_ic(:,:), grad_param(:)
    real(dp), intent(in) :: dir_ic(size(grad_ic,1),size(grad_ic,2))
    real(dp), intent(in) :: dir_param(size(grad_param))
    real(dp) :: ip
    ip = sum(grad_ic*dir_ic) + sum(grad_param*dir_param)
  end function evaluate_gradient

end module cost_module
