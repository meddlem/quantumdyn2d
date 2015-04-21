module evolve
  use constants
  use plotroutines
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(psi, x, y, V, n, M, A_d, A_u, A_conj)
    complex(dp), intent(inout) :: psi(:,:)
    complex(dp), intent(in)    :: A_d(:,:), A_u(:,:), A_conj(:,:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:), V(:,:)
    integer, intent(in)        :: n, M

    integer :: i

    do i=1,n
      call inc_time(psi, M, A_d, A_u, A_conj)
      
      if (mod(i,10)==0) then
        call plot_wavef(psi, x, y, V, M)
      endif
    enddo
  end subroutine

  subroutine inc_time(psi, M, A_d, A_u, A_conj)
    complex(dp), intent(inout) :: psi(:,:)
    complex(dp), intent(in)    :: A_d(:,:), A_u(:,:), A_conj(:,:,:)
    integer, intent(in) :: M

    complex(dp), allocatable :: r(:), A_d_tmp(:), A_l_tmp(:), A_u_tmp(:)
    integer ::  i, info

    allocate(A_d_tmp(M),A_l_tmp(M-1),A_u_tmp(M-1),r(M))
    
    ! modify lapack routine, so temp arrays are not needed? 

    ! horizontal sweep
    do i=1,M
      ! define needed temp arrays
      A_d_tmp = A_d(:,i)
      A_u_tmp = A_u(:,i)
      A_l_tmp = A_u(:,i)

      ! explicit part of calculation
      r = matmul(A_conj(:,:,i),psi(:,i))

      ! solve for psi at t=n+1/2
      call zgtsv(M,1,A_l_tmp,A_d_tmp,A_u_tmp,r,M,info)
      psi(:,i) = r
    enddo

    ! vertical sweep
    do i=1,M
      ! define needed temp arrays
      A_d_tmp = A_d(:,i)
      A_u_tmp = A_u(:,i)
      A_l_tmp = A_u(:,i)

      ! explicit part of calculation
      r = matmul(A_conj(:,:,i),psi(i,:))

      ! solve for psi at t=n+1/2
      call zgtsv(M,1,A_l_tmp,A_d_tmp,A_u_tmp,r,M,info)
      psi(i,:) = r
    enddo

    deallocate(A_d_tmp,A_l_tmp,A_u_tmp,r)
  end subroutine
end module
