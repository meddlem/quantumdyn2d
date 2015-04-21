module evolve
  use constants
  use plotroutines
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(psi, x, V, n, M, opp_d, opp_u, opm)
    complex(dp), intent(inout) :: psi(:)
    complex(dp), intent(in)    :: opp_d(:), opp_u(:), opm(:,:)
    real(dp), intent(in)       :: x(:), V(:)
    integer, intent(in)        :: n, M

    integer :: i

    do i=1,n
      call inc_time(psi, M, opp_d, opp_u, opm)
      call plot_wavef(psi, x, V, M)
    enddo
  end subroutine

  subroutine inc_time(psi, M, opp_d, opp_u, opm)
    complex(dp), intent(inout) :: psi(:)
    complex(dp), intent(in)    :: opp_d(:), opp_u(:), opm(:,:)
    integer, intent(in) :: M

    complex(dp), allocatable :: r(:), opp_d_tmp(:), opp_l_tmp(:), opp_u_tmp(:)
    integer ::  info

    allocate(opp_d_tmp(M),opp_l_tmp(M-1),opp_u_tmp(M-1),r(M))
    
    ! define temp arrays
    opp_d_tmp = opp_d
    opp_u_tmp = opp_u
    opp_l_tmp = opp_u
    
    ! explicit part of calculation
    r = matmul(opm,psi)

    ! enforce fixed bcs
    r(1) = (0._dp,0._dp)
    r(M) = (0._dp,0._dp)

    ! solve for psi at next timestep
    call zgtsv(M,1,opp_l_tmp,opp_d_tmp,opp_u_tmp,r,M,info)

    ! collect new psi 
    psi = r

    deallocate(opp_d_tmp,opp_l_tmp,opp_u_tmp,r)
  end subroutine
end module
