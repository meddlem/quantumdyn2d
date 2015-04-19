module evolve
  use constants
  use plotroutines
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(psi, x, L, n, M, opp_d, opp_dl, opp_du, opm)
    complex(dp), intent(inout) :: psi(:), opp_d(:), opp_dl(:), opp_du(:), &
      opm(:,:)
    real(dp), intent(in)       :: x(:), L
    integer, intent(in)        :: n, M

    integer :: i

    do i=1,n
      call inc_time(psi, M, opp_d, opp_dl, opp_du, opm)
!      call plot_wavef(psi, x, M)
    enddo
    call line_plot(x,abs(psi)**2,'x','P','','',2)

  end subroutine

  subroutine inc_time(psi, M, opp_d, opp_dl, opp_du, opm)
    complex(dp), intent(inout) :: psi(:), opp_d(:), opp_dl(:), opp_du(:),&
      opm(:,:)
    integer, intent(in) :: M

    complex(dp), allocatable :: psi_tmp(:), opp_d_tmp(:), opp_dl_tmp(:), &
      opp_du_tmp(:)
    integer ::  info

    allocate(opp_d_tmp(M),opp_dl_tmp(M),opp_du_tmp(M),psi_tmp(M))
    
    ! create temp arrays
    psi_tmp = matmul(opm,psi)
    opp_d_tmp = opp_d
    opp_dl_tmp = opp_dl
    opp_du_tmp = opp_du

    ! solve for psi at next timestep
    call dgtsv(M,1,opp_dl_tmp,opp_d_tmp,opp_du_tmp,psi_tmp,M,info)
    psi = psi_tmp

    deallocate(opp_d_tmp,opp_dl_tmp,opp_du_tmp,psi_tmp)
  end subroutine
end module
