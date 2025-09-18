module equations_module
  use constants_module, only: dp
  use variables_module, only: nx, ny, Lx, Ly, dx, dy, g, radius, u0, f0, h0, h1, pi, b, &
                               ihalo, is, ie, js, je, exchange_halo
  implicit none
contains

  !$FAD SKIP
  subroutine init_height(h, x, y, xoffset)
    use mpi_decomp_module, only: istart, iend, jstart, jend
    real(dp), intent(out) :: h(is:ie,js:je)
    real(dp), intent(in) :: x(is:ie), y(js:je)
    real(dp), intent(in), optional :: xoffset
    real(dp) :: x0, y0, r0, dist
    integer :: i,j
    x0 = 0.75d0*Lx
    if (present(xoffset)) then
      x0 = modulo(x0 + xoffset, Lx)
    end if
    y0 = 0.5d0*Ly
    r0 = radius/3.d0
    !$omp parallel do private(dist)
    do j = jstart, jend
       do i = istart, iend
          h(i,j) = h0
          dist = sqrt( (x(i)-x0)**2 + (y(j)-y0)**2 )
          if (dist < r0) then
             h(i,j) = h0 + 0.5d0*h1*(1.d0+cos(pi*dist/r0))
          end if
       end do
    end do
    !$omp end parallel do
    call exchange_halo(h)
  end subroutine init_height

  !$FAD SKIP
  subroutine velocity_field(u, v, x, y)
    use mpi_decomp_module, only: istart, iend, jstart, jend
    real(dp), intent(out) :: u(is:ie,js:je), v(is:ie,js:jend+1)
    real(dp), intent(in) :: x(is:ie), y(js:je)
    integer :: i,j
    !$omp parallel do
    do j = jstart, jend
      do i = istart, iend
         u(i,j) = u0
      end do
    end do
    !$omp end parallel do
    !$omp parallel workshare
    v(:,:) = 0.d0
    !$omp end parallel workshare
    call exchange_halo(u)
    call exchange_halo(v)
  end subroutine velocity_field

  !$FAD SKIP
  subroutine init_topography(b, x, y)
    use mpi_decomp_module, only: istart, iend, jstart, jend
    real(dp), intent(out) :: b(is:ie,js:je)
    real(dp), intent(in) :: x(is:ie), y(js:je)
    real(dp) :: x0, y0, r0, dist
    integer :: i,j
    x0 = 0.5d0*Lx
    y0 = 0.75d0*Ly
    r0 = radius/4.d0
    !$omp parallel do private(dist)
    do j = jstart, jend
       do i = istart, iend
          dist = sqrt((x(i)-x0)**2 + (y(j)-y0)**2)
          if (dist < r0) then
             b(i,j) = h1 * (1.d0 - dist/r0)
          else
             b(i,j) = 0.d0
          end if
       end do
    end do
    !$omp end parallel do
    call exchange_halo(b)
  end subroutine init_topography

  !$FAD CONSTANT_VARS: y
  subroutine init_geostrophic_height(h, y)
    use mpi_decomp_module, only: istart, iend, jstart, jend
    real(dp), intent(out) :: h(is:ie,js:je)
    real(dp), intent(in) :: y(js:je)
    integer :: i, j
    real(dp), parameter :: coeff = f0 * u0 * radius / g
    !$omp parallel do
    do j = jstart, jend
       do i = istart, iend
          h(i,j) = h0 + coeff * sin(y(j)/radius)**2
       end do
    end do
    !$omp end parallel do
    call exchange_halo(h)
  end subroutine init_geostrophic_height

  subroutine geostrophic_velocity(u, v, h)
    use mpi_decomp_module, only: istart, iend, jstart, jend
    real(dp), intent(out) :: u(is:ie,js:je), v(is:ie,js:jend+1)
    real(dp), intent(in)  :: h(is:ie,js:je)
    integer :: i, j
    integer :: jp1, jm1
    !$omp parallel do private(jp1,jm1)
    do j = jstart, jend
       jp1 = min(j+1, ny)
       jm1 = max(j-1, 1)
       do i = istart, iend
          u(i,j) = - g / f0 * ((h(i-1,jp1) + h(i,jp1)) - (h(i-1,jm1) + h(i,jm1))) / (4.0d0 * dy)
       end do
    end do
    !$omp end parallel do
    !$omp parallel do private(jm1)
    do j = max(jstart, 2), jend
       jm1 = j - 1
       do i = istart, iend
          v(i,j) = g / f0 * ((h(i+1,jm1) + h(i+1,j)) - (h(i-1,jm1) + h(i-1,j))) / (4.0d0 * dx)
       end do
    end do
    !$omp end parallel do
    call exchange_halo(u)
    call exchange_halo(v)
    if (jstart == 1) then
       v(:,1) = 0.0_dp
    end if
    if (jend == ny) then
       v(:,ny+1) = 0.0_dp
    end if
  end subroutine geostrophic_velocity

  !$FAD SKIP
  subroutine analytic_height(ha, x, y, t)
    real(dp), intent(out) :: ha(is:ie,js:je)
    real(dp), intent(in) :: x(is:ie), y(ny), t
    call init_height(ha, x, y, t*u0)
  end subroutine analytic_height

  subroutine rhs(h, u, v, dhdt, dudt, dvdt, no_momentum_tendency)
    use mpi_decomp_module, only: istart, iend, jstart, jend
    real(dp), intent(in) :: h(is:ie,js:je), u(is:ie,js:je), v(is:ie,js:jend+1)
    real(dp), intent(out) :: dhdt(is:ie,js:je)
    real(dp), intent(out) :: dudt(is:ie,js:je), dvdt(is:ie,js:jend+1)
    logical, intent(in), optional :: no_momentum_tendency
    integer :: i,j,jp1,jm1
    real(dp) :: fe,fw,fn,fs,ue,uw,vn,vs
    real(dp) :: h_e, h_w, h_n, h_s, v_avg, u_avg

    ! continuity equation
    !$omp parallel do private(jp1,jm1,ue,uw,vn,vs,fe,fw,fn,fs)
    do j = jstart, jend
       jp1 = min(j+1,ny)
       jm1 = max(j-1,1)
       do i = istart, iend
          ue = u(i+1,j)
          uw = u(i,j)
          if (ue > 0.d0) then
             fe = ue * h(i,j)
          else
             fe = ue * h(i+1,j)
          end if
          if (uw > 0.d0) then
             fw = uw * h(i-1,j)
          else
             fw = uw * h(i,j)
          end if
          vn = v(i,j+1)
          vs = v(i,j)
          if (vn > 0.d0) then
             fn = vn * h(i,j)
          else
             fn = vn * h(i,jp1)
          end if
          if (vs > 0.d0) then
             fs = vs * h(i,jm1)
          else
             fs = vs * h(i,j)
          end if
          dhdt(i,j) = -((fe - fw)/dx + (fn - fs)/dy)
       end do
    end do
    !$omp end parallel do

    if (present(no_momentum_tendency)) then
       if (no_momentum_tendency) then
          !$omp parallel workshare
          dudt = 0.d0
          dvdt = 0.d0
          !$omp end parallel workshare
          return
       end if
    end if

    ! zonal momentum
    !$omp parallel do private(jp1,jm1,h_e,h_w,v_avg)
    do j = jstart, jend
       jp1 = min(j+1,ny)
       jm1 = max(j-1,1)
       do i = istart, iend
          h_e = h(i,j) + b(i,j)
          h_w = h(i-1,j) + b(i-1,j)
          v_avg = 0.25d0*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1))
          dudt(i,j) = - u(i,j) * (u(i+1,j) - u(i-1,j)) / (2.d0*dx) &
                      - v_avg * (u(i,jp1) - u(i,jm1)) / (2.d0*dy) &
                      - g * (h_e - h_w) / dx &
                      + f0 * v_avg
       end do
    end do
    !$omp end parallel do

    ! meridional momentum
    !$omp parallel do private(jp1,jm1,h_n,h_s,u_avg)
    do j = max(jstart, 2), jend
       jp1 = j+1
       jm1 = j-1
       do i = istart, iend
          h_n = h(i,j) + b(i,j)
          h_s = h(i,jm1) + b(i,jm1)
          u_avg = 0.25d0*(u(i,jm1) + u(i+1,jm1) + u(i,j) + u(i+1,j))
          dvdt(i,j) = - u_avg * (v(i+1,j) - v(i-1,j)) / (2.d0*dx) &
                      - v(i,j) * (v(i,jp1) - v(i,jm1)) / (2.d0*dy) &
                      - g * (h_n - h_s) / dy &
                      - f0 * u_avg
       end do
    end do
    !$omp end parallel do
  end subroutine rhs

end module equations_module
