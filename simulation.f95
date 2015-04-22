module simulation
  use constants
  use plotroutines
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(psi, x, V, n, M, A_x)
    complex(dp), intent(inout) :: psi(:)
    complex(dp), intent(in)    :: A_x(:,:)
    real(dp), intent(in)       :: x(:), V(:)
    integer, intent(in)        :: n, M

    integer :: i

    do i=1,n
      call solve_nxt(psi, M, A_x)
      call plot_wavef(psi, x, V, M)
    enddo
  end subroutine

  subroutine solve_nxt(psi, M, A_x)
    complex(dp), intent(inout) :: psi(:)
    complex(dp), intent(in)    :: A_x(:,:)
    integer, intent(in)        :: M

    complex(dp), allocatable :: g(:), A_x_d_tmp(:), A_x_l_tmp(:), A_x_u_tmp(:)
    integer                  ::  info

    allocate(A_x_d_tmp(M), A_x_l_tmp(M-1), A_x_u_tmp(M-1), g(M))
    
    ! init temp arrays
    A_x_l_tmp = A_x(1,1:M-1)
    A_x_d_tmp = A_x(2,:)
    A_x_u_tmp = A_x(1,1:M-1)
    
    ! explicit part of calculation, mat-vec multiplication
    call zgbmv('N',M,M,1,1,one,conjg(A_x),3,psi,1,zero,g,1)

    ! solve for wavefunction at t=n+1
    call zgtsv(M,1,A_x_l_tmp,A_x_d_tmp,A_x_u_tmp,g,M,info)

    ! collect wavefunction at t=n+1
    psi = g

    deallocate(A_x_d_tmp, A_x_l_tmp, A_x_u_tmp, g)
  end subroutine
end module
