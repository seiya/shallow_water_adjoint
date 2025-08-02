program shallow_water_forward
  use constants_module_ad, only: dp
  use cost_module_ad, only: calc_mse_fwd_ad
  implicit none
  real(dp) :: a(2,2), a_ad(2,2)
  real(dp) :: b(2,2), b_ad(2,2)
  real(dp) :: mse, mse_ad
  a = reshape([1.d0,2.d0,3.d0,4.d0],[2,2])
  b = reshape([1.5d0,2.5d0,3.5d0,4.5d0],[2,2])
  a_ad = 1.d0
  b_ad = 0.d0
  call calc_mse_fwd_ad(a, a_ad, b, b_ad, mse, mse_ad)
  print *, 'mse =', mse
  print *, 'mse_ad =', mse_ad
end program shallow_water_forward
