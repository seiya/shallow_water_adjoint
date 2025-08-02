program shallow_water_reverse
  use constants_module_ad, only: dp
  use cost_module_ad, only: calc_mse_rev_ad
  implicit none
  real(dp) :: a(2,2), a_ad(2,2)
  real(dp) :: b(2,2), b_ad(2,2)
  real(dp) :: mse_ad
  a = reshape([1.d0,2.d0,3.d0,4.d0],[2,2])
  b = reshape([1.5d0,2.5d0,3.5d0,4.5d0],[2,2])
  mse_ad = 1.d0
  call calc_mse_rev_ad(a, a_ad, b, b_ad, mse_ad)
  print *, 'a_ad =', a_ad
  print *, 'b_ad =', b_ad
end program shallow_water_reverse
