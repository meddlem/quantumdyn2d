module initialize
  use constants
  implicit none
  private
  public :: init_wavef, init_V, init_ops

contains
  
  subroutine init_wavef(psi_0,x,dx,k,L,M)
    complex(dp), intent(inout) :: psi_0(:) 
    real(dp), intent(inout)    :: x(:)
    real(dp), intent(in)       :: dx, k, L
    integer, intent(in)        :: M

    integer :: i
    
    do i = 1,M
      x(i) = i*dx
    enddo

    psi_0 = 0._dp
    where(x<L/2) psi_0 = 1/sqrt(L/2._dp)*exp(cmplx(0._dp,k*x,dp))
  end subroutine

  subroutine init_V(V,x,L)
    real(dp), intent(in) :: x(:), L
    real(dp), intent(inout) :: V(:)
    
    V = 0._dp
    ! block potential
!    where(x>3*L/4) V = 1._dp
  end subroutine
    
  subroutine init_ops(opp_d,opp_dl,opp_du,opm,V,dt,dx,M)
    complex(dp), intent(inout) :: opp_d(:), opp_dl(:), opp_du(:), opm(:,:)
    real(dp), intent(in)       :: V(:), dt, dx
    integer, intent(in)        :: M

    integer :: i

    opm = 0._dp
    ! construct operators 
    do i = 1,M
      opm(i,i) = cmplx(1._dp , (-dt*2._dp/(dx**2) - V(i))/2._dp , dp)
      opp_d(i) = cmplx(1._dp , (dt*2._dp/(dx**2) + V(i))/2._dp , dp)

      if (i>1) then
        opm(i,i-1) = cmplx(0._dp , dt*0.5_dp/(dx**2) , dp)
        opp_dl(i-1) = cmplx(0._dp , -dt*0.5_dp/(dx**2) , dp)
      endif

      if (i<M) then
        opm(i,i+1) = cmplx(0._dp , dt*0.5_dp/(dx**2) , dp)
        opp_du(i) = cmplx(0._dp , -dt*0.5_dp/(dx**2) , dp)
      endif
    enddo
  end subroutine
end module 
