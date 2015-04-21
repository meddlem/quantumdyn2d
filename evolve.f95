module evolve
  use constants
  use plotroutines
  use omp_lib
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(psi, x, y, V, n, M, A_x_d, A_y_d, A_u, A_x_conj, &
      A_y_conj)
    complex(dp), intent(inout) :: psi(:,:)
    complex(dp), intent(in)    :: A_x_d(:,:), A_y_d(:,:), A_u(:,:), &
      A_x_conj(:,:,:), A_y_conj(:,:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:), V(:,:)
    integer, intent(in)        :: n, M

    integer :: i

    do i=1,n
      call inc_time(psi, M, A_x_d, A_y_d, A_u, A_x_conj, A_y_conj)
      
      if (mod(i,30)==0) then
        call plot_wavef(psi, x, y, V, M, M)
      endif
    enddo
  end subroutine

  subroutine inc_time(psi, M, A_x_d, A_y_d, A_u, A_x_conj, A_y_conj)
    complex(dp), intent(inout) :: psi(:,:)
    complex(dp), intent(in)    :: A_x_d(:,:), A_y_d(:,:), A_u(:,:), &
      A_x_conj(:,:,:), A_y_conj(:,:,:) 
    integer, intent(in) :: M

    complex(dp), allocatable :: g(:), A_d_tmp(:), A_l_tmp(:), A_u_tmp(:)
    integer ::  i, info

    allocate(A_d_tmp(M),A_l_tmp(M-1),A_u_tmp(M-1),g(M))
    
    ! modify lapack routine, so temp arrays are not needed? 

    ! horizontal sweep
    !$omp parallel do private(A_d_tmp,A_u_tmp,A_l_tmp,g)
    do i=1,M
      ! define needed temp arrays
      A_d_tmp = A_x_d(:,i)
      A_u_tmp = A_u(:,i)
      A_l_tmp = A_u(:,i)

      ! explicit part of calculation
      g = matmul(A_x_conj(:,:,i),psi(:,i))

      ! solve for psi at t=n+1/2
      call zgtsv(M,1,A_l_tmp,A_d_tmp,A_u_tmp,g,M,info)
      psi(:,i) = g
    enddo
    !$omp end parallel do

    ! vertical sweep
    !$omp parallel do private(A_d_tmp,A_u_tmp,A_l_tmp,g)
    do i=1,M
      ! define needed temp arrays
      A_d_tmp = A_y_d(:,i)
      A_u_tmp = A_u(:,i)
      A_l_tmp = A_u(:,i)

      ! explicit part of calculation
      g = matmul(A_y_conj(:,:,i),psi(i,:))

      ! solve for psi at t=n+1
      call zgtsv(M,1,A_l_tmp,A_d_tmp,A_u_tmp,g,M,info)
      psi(i,:) = g
    enddo
    !$omp end parallel do

    deallocate(A_d_tmp,A_l_tmp,A_u_tmp,g)
  end subroutine
end module
