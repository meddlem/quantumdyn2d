module evolve
  use constants
  use plotroutines
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(psi, x, y, V, n, M, A_d, A_u, A_conj)
    complex(dp), intent(inout) :: psi(:,:)
    complex(dp), intent(in)    :: A_d(:), A_u(:), A_conj(:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:), V(:,:)
    integer, intent(in)        :: n, M

    integer :: i

    do i=1,n
      call inc_time(psi, M, A_d, A_u, A_conj)
      call plot_wavef(psi, x, V, M)
    enddo
  end subroutine

  subroutine inc_time(psi, M, A_d, A_u, A_conj)
    complex(dp), intent(inout) :: psi(:,:)
    complex(dp), intent(in)    :: A_d(:), A_u(:), A_conj(:,:)
    integer, intent(in) :: M

    complex(dp), allocatable :: r(:), A_d_tmp(:), A_l_tmp(:), A_u_tmp(:)
    integer ::  info

    allocate(A_d_tmp(M),A_l_tmp(M-1),A_u_tmp(M-1),r(M))
    
    ! define temp arrays
    A_d_tmp = A_d
    A_u_tmp = A_u
    A_l_tmp = A_u
    
    ! explicit part of calculation
    r = matmul(A_conj,psi)

    ! enforce fixed bcs
    r(1) = (0._dp,0._dp)
    r(M) = (0._dp,0._dp)

    ! solve for psi at next timestep
    call zgtsv(M,1,A_l_tmp,A_d_tmp,A_u_tmp,r,M,info)

    ! collect new psi 
    psi = r

    deallocate(A_d_tmp,A_l_tmp,A_u_tmp,r)
  end subroutine
end module
