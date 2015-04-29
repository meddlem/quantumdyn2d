module simulation
  use constants
  use plotroutines
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(psi, x, V, n, M, Ax)
    complex(dp), intent(inout) :: psi(:)
    complex(dp), intent(in)    :: Ax(:,:)
    real(dp), intent(in)       :: x(:), V(:)
    integer, intent(in)        :: n, M

    integer :: i

    do i=1,n
      call solve_nxt(psi, M, Ax)
      if(mod(i,10)==0) call plot_wavef(psi, x, V, M)
    enddo
  end subroutine

  subroutine solve_nxt(psi, M, Ax)
    complex(dp), intent(inout) :: psi(:)
    complex(dp), intent(in)    :: Ax(:,:)
    integer, intent(in)        :: M

    complex(dp), allocatable :: g(:), Ax_d_tmp(:), Ax_l_tmp(:), Ax_u_tmp(:)
    integer                  ::  info

    allocate(Ax_d_tmp(M), Ax_l_tmp(M-1), Ax_u_tmp(M-1), g(M))
    
    ! init temp arrays
    Ax_l_tmp = Ax(1,1:M-1)
    Ax_d_tmp = Ax(2,:)
    Ax_u_tmp = Ax(1,1:M-1)
    
    ! explicit part of calculation, mat-vec multiplication
    call zgbmv('N',M,M,1,1,one,conjg(Ax),3,psi,1,zero,g,1)

    ! solve for wavefunction at t=n+1
    call zgtsv(M,1,Ax_l_tmp,Ax_d_tmp,Ax_u_tmp,g,M,info)

    ! collect wavefunction at t=n+1
    psi = g

    deallocate(Ax_d_tmp, Ax_l_tmp, Ax_u_tmp, g)
  end subroutine
end module
