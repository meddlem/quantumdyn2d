module simulation
  use constants
  use plotroutines
  use omp_lib
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(psi, x, y, n, M_x, M_y, A_x, A_y)
    complex(dp), intent(inout) :: psi(:,:)
    complex(dp), intent(in)    :: A_x(:,:,:), A_y(:,:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:)
    integer, intent(in)        :: n, M_x, M_y

    integer :: i

    do i=1,n
      call solve_nxt(psi, M_x, M_y, A_x, A_y)
      
      if (mod(i,100)==0) call plot_wavef(psi, x, y, M_x, M_y)
    enddo
  end subroutine

  subroutine solve_nxt(psi, M_x, M_y, A_x, A_y)
    complex(dp), intent(inout) :: psi(:,:)
    complex(dp), intent(in)    :: A_x(:,:,:), A_y(:,:,:)
    integer, intent(in)        :: M_x, M_y

    complex(dp), allocatable :: g_x(:), g_y(:), A_x_d(:), A_x_l(:), &
                                A_x_u(:), A_y_d(:), A_y_l(:), &
                                A_y_u(:)
    integer                  :: i, info

    allocate(A_x_d(M_x), A_x_l(M_x-1), A_x_u(M_x-1), g_x(M_x), A_y_d(M_y), &
      A_y_l(M_y-1), A_y_u(M_y-1), g_y(M_y))

    ! horizontal sweep
    !$omp parallel do private(A_x_d,A_x_u,A_x_l,g_x)
    do i=1,M_y
      ! init temp arrays
      g_x = A_x(2,:,i)
      A_x_d = A_x(2,:,i)
      A_x_u = A_x(1,1:M_x-1,i)
      A_x_l = A_x(1,1:M_x-1,i)

      ! explicit part of calculation, mat-vec multiplication
      call zgbmv('N', M_x, M_x, 1, 1, one, conjg(A_x(:,:,i)), 3, psi(:,i), 1, &
        zero, g_x, 1)

      ! solve for psi at t=n+1/2
      call zgtsv(M_x, 1, A_x_l, A_x_d, A_x_u, g_x, M_x, info)
      psi(:,i) = g_x
    enddo
    !$omp end parallel do

    ! vertical sweep
    !$omp parallel do private(A_y_d,A_y_u,A_y_l,g_y)
    do i=1,M_x
      ! init temp arrays
      g_y = A_y(2,:,i)
      A_y_d = A_y(2,:,i)
      A_y_u = A_y(1,1:M_y-1,i)
      A_y_l = A_y(1,1:M_y-1,i)

      ! explicit part of calculation, mat-vec multiplication
      call zgbmv('N', M_y, M_y, 1, 1, one, conjg(A_y(:,:,i)), 3, psi(i,:), 1, &
        zero, g_y, 1)

      ! solve for psi at t=n+1
      call zgtsv(M_y, 1, A_y_l, A_y_d, A_y_u, g_y, M_y, info)
      psi(i,:) = g_y
    enddo
    !$omp end parallel do

    deallocate(A_x_d, A_x_l, A_x_u, A_y_d, A_y_l, A_y_u, g_x, g_y)
  end subroutine
end module
