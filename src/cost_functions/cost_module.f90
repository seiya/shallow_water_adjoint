module cost_module
  use constants_module, only: dp
  implicit none
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

  !> Compute L1, L2, and Linf error norms
  subroutine calc_error_norms(height_num, height_ana, lat, l1err, l2err, maxerr)
    real(dp), intent(in) :: height_num(:,:), height_ana(:,:), lat(:)
    real(dp), intent(out) :: l1err, l2err, maxerr
    integer :: i, j, nlon, nlat
    real(dp) :: err, w, wtsum

    nlon = size(height_num, 1)
    nlat = size(height_num, 2)
    maxerr = 0.d0
    l1err = 0.d0
    l2err = 0.d0
    wtsum = 0.d0

    do j = 1, nlat
       w = cos(lat(j))
       wtsum = wtsum + w
       do i = 1, nlon
          err = height_num(i,j) - height_ana(i,j)
          maxerr = max(maxerr, abs(err))
          l1err = l1err + abs(err)*w
          l2err = l2err + err*err*w
       end do
    end do

    l1err = l1err/(nlon*wtsum)
    l2err = sqrt(l2err/(nlon*wtsum))
  end subroutine calc_error_norms

  !> Evaluate inner product of gradients with perturbation directions
  function evaluate_gradient(grad_ic, grad_param, dir_ic, dir_param) result(ip)
    real(dp), intent(in) :: grad_ic(:,:), grad_param(:)
    real(dp), intent(in) :: dir_ic(size(grad_ic,1),size(grad_ic,2))
    real(dp), intent(in) :: dir_param(size(grad_param))
    real(dp) :: ip
    ip = sum(grad_ic*dir_ic) + sum(grad_param*dir_param)
  end function evaluate_gradient

end module cost_module
