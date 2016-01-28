subroutine calc_velocity(nx, ny, x, y, alpha, u)
  use omp_lib
  implicit none
  integer :: nx, ny
  real(kind=8), dimension(nx,ny), intent(in) :: x, y, alpha
  real(kind=8), dimension(2,nx,ny), intent(out) :: u

  integer :: i, j, ii, jj, maxthreads
  real(kind=8) :: delta, dx, dy, r2, rt2, u_tmp_1, u_tmp_2, pi
  pi = 4.*atan(1.)
  delta = 1.0
  !$OMP PARALLEL
  !$OMP DO !PRIVATE(u_tmp_1,u_tmp_2, dx, dy, r2, rt2)
  do j=1,ny  
     do i=1,nx
        u_tmp_1 = 0.0
        u_tmp_2 = 0.0
        do jj=1,ny
           do ii=1,nx
              dx = x(i,j) - x(ii,jj)
              dy = y(i,j) - y(ii,jj)
              r2 = dx**2 + dy**2
              rt2 = r2/delta**2
              u_tmp_1 = u_tmp_1 + dy*alpha(ii,jj)/(2.0*pi*(r2 + 1e-16))*(1.0 - exp(-rt2))
              u_tmp_2 = u_tmp_2 -dx*alpha(ii,jj)/(2.0*pi*(r2 + 1e-16))*(1.0 - exp(-rt2))
           end do
        end do
        u(1,i,j) = u_tmp_1
        u(2,i,j) = u_tmp_2
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine calc_velocity
