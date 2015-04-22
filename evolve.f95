module evolve
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
      call inc_time(psi, M_x, M_y, A_x, A_y)
      
      if (mod(i,30)==0) call plot_wavef(psi, x, y, M_x, M_y)
    enddo
  end subroutine

  subroutine inc_time(psi, M_x, M_y, A_x, A_y)
    complex(dp), intent(inout) :: psi(:,:)
    complex(dp), intent(in)    :: A_x(:,:,:), A_y(:,:,:)
    integer, intent(in)        :: M_x, M_y

    complex(dp), allocatable :: g_x(:), g_y(:), A_x_d_tmp(:), A_x_l_tmp(:), &
                                A_x_u_tmp(:), A_y_d_tmp(:), A_y_l_tmp(:), &
                                A_y_u_tmp(:)
    integer                  :: i, info

    allocate(A_x_d_tmp(M_x), A_x_l_tmp(M_x-1), A_x_u_tmp(M_x-1), g_x(M_x), &
      A_y_d_tmp(M_y), A_y_l_tmp(M_y-1), A_y_u_tmp(M_y-1), g_y(M_y))

    ! horizontal sweep
    !$omp parallel do private(A_x_d_tmp,A_x_u_tmp,A_x_l_tmp,g_x)
    do i=1,M_y
      ! init temp arrays
      g_x = A_x(2,:,i)
      A_x_d_tmp = A_x(2,:,i)
      A_x_u_tmp = A_x(1,1:M_x-1,i)
      A_x_l_tmp = A_x(1,1:M_x-1,i)

      ! explicit part of calculation, mat-vec multiplication
      call zgbmv('N',M_x,M_x,1,1,one,conjg(A_x(:,:,i)),3,psi(:,i),1,zero,g_x,1)

      ! solve for psi at t=n+1/2
      call zgtsv(M_x,1,A_x_l_tmp,A_x_d_tmp,A_x_u_tmp,g_x,M_x,info)
      psi(:,i) = g_x
    enddo
    !$omp end parallel do

    ! vertical sweep
    !$omp parallel do private(A_y_d_tmp,A_y_u_tmp,A_y_l_tmp,g_y)
    do i=1,M_x
      ! init temp arrays
      g_y = A_y(2,:,i)
      A_y_d_tmp = A_y(2,:,i)
      A_y_u_tmp = A_y(1,1:M_y-1,i)
      A_y_l_tmp = A_y(1,1:M_y-1,i)

      ! explicit part of calculation, mat-vec multiplication
      call zgbmv('N',M_y,M_y,1,1,one,conjg(A_y(:,:,i)),3,psi(i,:),1,zero,g_y,1)

      ! solve for psi at t=n+1
      call zgtsv(M_y,1,A_y_l_tmp,A_y_d_tmp,A_y_u_tmp,g_y,M_y,info)
      psi(i,:) = g_y
    enddo
    !$omp end parallel do

    deallocate(A_x_d_tmp,A_x_l_tmp,A_x_u_tmp,A_y_d_tmp,A_y_l_tmp,A_y_u_tmp,&
      g_x,g_y)
  end subroutine
end module
