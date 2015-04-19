module evolve
  use constants
  use plotroutines
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(psi, x, L, n, M, opp_d, opp_du, opm)
    complex(dp), intent(inout) :: psi(:)
    complex(dp), intent(in)    :: opp_d(:), opp_du(:), opm(:,:)
    real(dp), intent(in)       :: x(:), L
    integer, intent(in)        :: n, M

    integer :: i

    do i=1,n
      call inc_time(psi, M, opp_d, opp_du, opm)
      call plot_wavef(psi, x, M)
    enddo
  end subroutine

  subroutine inc_time(psi, M, opp_d, opp_du, opm)
    complex(dp), intent(inout) :: psi(:)
    complex(dp), intent(in)    :: opp_d(:), opp_du(:), opm(:,:)
    integer, intent(in) :: M

    complex(dp), allocatable :: r(:), opp_d_tmp(:), opp_dl_tmp(:), opp_du_tmp(:)
    integer ::  info

    allocate(opp_d_tmp(M),opp_dl_tmp(M-1),opp_du_tmp(M-1),r(M))
    
    ! create temp arrays
    r = matmul(opm,psi)
    r(1) = (0._dp,0._dp)
    r(M) = (0._dp,0._dp)

    opp_d_tmp = opp_d
    opp_du_tmp = opp_du
    opp_dl_tmp = opp_du

    ! solve for psi at next timestep
    call zgtsv(M,1,opp_dl_tmp,opp_d_tmp,opp_du_tmp,r,M,info)
    ! collect new psi 
    psi = r

    deallocate(opp_d_tmp,opp_du_tmp,r)
  end subroutine
end module
