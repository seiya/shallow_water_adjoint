module equations_module
  use constants_module, only: dp
  use variables_module, only: nx, ny, Lx, Ly, dx, dy, g, radius, u0, f0, h0, h1, pi, b
  implicit none
contains

  !$FAD SKIP
  subroutine init_height(h, x, y, xoffset)
    real(dp), intent(out) :: h(nx,ny)
    real(dp), intent(in) :: x(nx), y(ny)
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
    do j = 1, ny
       do i = 1, nx
          h(i,j) = h0
          dist = sqrt( (x(i)-x0)**2 + (y(j)-y0)**2 )
          if (dist < r0) then
             h(i,j) = h0 + 0.5d0*h1*(1.d0+cos(pi*dist/r0))
          end if
       end do
    end do
    !$omp end parallel do
  end subroutine init_height

  !$FAD SKIP
  subroutine velocity_field(u, v, x, y)
    real(dp), intent(out) :: u(nx,ny), v(nx,ny+1)
    real(dp), intent(in) :: x(nx), y(ny)
    integer :: i,j
    !$omp parallel do
    do j = 1, ny
      do i = 1, nx
         u(i,j) = u0
      end do
    end do
    !$omp end parallel do
    !$omp parallel workshare
    v(:,:) = 0.d0
    !$omp end parallel workshare
  end subroutine velocity_field

  !$FAD SKIP
  subroutine init_topography(b, x, y)
    real(dp), intent(out) :: b(nx,ny)
    real(dp), intent(in) :: x(nx), y(ny)
    real(dp) :: x0, y0, r0, dist
    integer :: i,j
    x0 = 0.5d0*Lx
    y0 = 0.75d0*Ly
    r0 = radius/4.d0
    !$omp parallel do private(dist)
    do j = 1, ny
       do i = 1, nx
          dist = sqrt((x(i)-x0)**2 + (y(j)-y0)**2)
          if (dist < r0) then
             b(i,j) = h1 * (1.d0 - dist/r0)
          else
             b(i,j) = 0.d0
          end if
       end do
    end do
    !$omp end parallel do
  end subroutine init_topography

  !$FAD CONSTANT_VARS: y
  subroutine init_geostrophic_height(h, y)
    real(dp), intent(out) :: h(nx,ny)
    real(dp), intent(in) :: y(ny)
    integer :: i, j
    real(dp), parameter :: coeff = f0 * u0 * radius / g
    !$omp parallel do
    do j = 1, ny
       do i = 1, nx
          h(i,j) = h0 + coeff * sin(y(j)/radius)**2
       end do
    end do
    !$omp end parallel do
  end subroutine init_geostrophic_height

  subroutine geostrophic_velocity(u, v, h)
    real(dp), intent(out) :: u(nx,ny), v(nx,ny+1)
    real(dp), intent(in)  :: h(nx,ny)
    integer :: i, j
    integer :: ip1, im1
    integer :: jp1, jm1
    !$omp parallel do private(im1,jp1,jm1)
    do j = 1, ny
       jp1 = min(j+1, ny)
       jm1 = max(j-1, 1)
       do i = 1, nx
          im1 = modulo(i-2, nx) + 1
          u(i,j) = - g / f0 * ((h(im1,jp1) + h(i,jp1)) - (h(im1,jm1) + h(i,jm1))) / (4.0d0 * dy)
       end do
    end do
    !$omp end parallel do
    !$omp parallel do private(ip1,im1,jm1)
    do j = 2, ny - 1
       jm1 = j - 1
       do i = 1, nx
          ip1 = modulo(i, nx) + 1
          im1 = modulo(i-2, nx) + 1
          v(i,j) = g / f0 * ((h(ip1,jm1) + h(ip1,j)) - (h(im1,jm1) + h(im1,j))) / (4.0d0 * dx)
       end do
    end do
    !$omp end parallel do
    v(:,1) = 0.0_dp
    v(:,ny+1) = 0.0_dp
  end subroutine geostrophic_velocity

  !$FAD SKIP
  subroutine analytic_height(ha, x, y, t)
    real(dp), intent(out) :: ha(nx,ny)
    real(dp), intent(in) :: x(nx), y(ny), t
    call init_height(ha, x, y, t*u0)
  end subroutine analytic_height

  subroutine rhs(h, u, v, dhdt, dudt, dvdt, no_momentum_tendency)
    real(dp), intent(in) :: h(nx,ny), u(nx,ny), v(nx,ny+1)
    real(dp), intent(out) :: dhdt(nx,ny)
    real(dp), intent(out) :: dudt(nx,ny), dvdt(nx,ny+1)
    logical, intent(in), optional :: no_momentum_tendency
    integer :: i, j
    integer :: ip1, im1
    integer :: jp1, jm1
    real(dp) :: fe,fw,fn,fs,ue,uw,vn,vs
    real(dp) :: h_e, h_w, h_n, h_s, v_avg, u_avg

    ! continuity equation
    !$omp parallel do private(ip1,im1,jp1,jm1,ue,uw,vn,vs,fe,fw,fn,fs)
    do j = 1, ny
       jp1 = min(j+1,ny)
       jm1 = max(j-1,1)
       do i = 1, nx
          ip1 = modulo(i,nx) + 1
          im1 = modulo(i-2,nx) + 1
          ue = u(ip1,j)
          uw = u(i,j)
          if (ue > 0.d0) then
             fe = ue * h(i,j)
          else
             fe = ue * h(ip1,j)
          end if
          if (uw > 0.d0) then
             fw = uw * h(im1,j)
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
    !$omp parallel do private(ip1,im1,jp1,jm1,h_e,h_w,v_avg)
    do j = 1, ny
       jp1 = min(j+1,ny)
       jm1 = max(j-1,1)
       do i = 1, nx
          ip1 = modulo(i,nx) + 1
          im1 = modulo(i-2,nx) + 1
          h_e = h(i,j) + b(i,j)
          h_w = h(im1,j) + b(im1,j)
          v_avg = 0.25d0*(v(im1,j) + v(im1,j+1) + v(i,j) + v(i,j+1))
          dudt(i,j) = - u(i,j) * (u(ip1,j) - u(im1,j)) / (2.d0*dx) &
                      - v_avg * (u(i,jp1) - u(i,jm1)) / (2.d0*dy) &
                      - g * (h_e - h_w) / dx &
                      + f0 * v_avg
       end do
    end do
    !$omp end parallel do

    ! meridional momentum
    !$omp parallel do private(ip1,im1,jp1,jm1,h_n,h_s,u_avg)
    do j = 2, ny - 1
       jp1 = j+1
       jm1 = j-1
       do i = 1, nx
          ip1 = modulo(i,nx) + 1
          im1 = modulo(i-2,nx) + 1
          h_n = h(i,j) + b(i,j)
          h_s = h(i,jm1) + b(i,jm1)
          u_avg = 0.25d0*(u(i,jm1) + u(ip1,jm1) + u(i,j) + u(ip1,j))
          dvdt(i,j) = - u_avg * (v(ip1,j) - v(im1,j)) / (2.d0*dx) &
                      - v(i,j) * (v(i,jp1) - v(i,jm1)) / (2.d0*dy) &
                      - g * (h_n - h_s) / dy &
                      - f0 * u_avg
       end do
    end do
    !$omp end parallel do
  end subroutine rhs

end module equations_module
