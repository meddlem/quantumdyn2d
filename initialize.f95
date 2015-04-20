module initialize
  use constants
  implicit none
  private
  public :: init_wavef, init_V, init_ops

contains
  
  subroutine init_wavef(psi_0,x,dx,L,k,M)
    complex(dp), intent(inout) :: psi_0(:) 
    real(dp), intent(inout)    :: x(:)
    real(dp), intent(in)       :: dx, L, k
    integer, intent(in)        :: M
    integer :: i
    
    do i = 1,M
      x(i) = i*dx
    enddo

    ! isq well
    ! psi_0 = sin(2*pi*x/L)*exp(cmplx(0._dp,k*x,dp))
    
    ! gaussian wavepacket
    psi_0 = exp(-1.0_dp*(x-L/3)**2)*exp(cmplx(0._dp,k*x,dp))

    ! normalize wavefunction
    psi_0 = psi_0/sqrt(sum(abs(psi_0)**2*dx))
  end subroutine

  subroutine init_V(V,x,L)
    real(dp), intent(in) :: x(:), L
    real(dp), intent(inout) :: V(:)
    
    ! block/constant potential
    V = 0._dp
    where(x>78._dp .and. x<85._dp) V = 1._dp
    
    ! harmonic potential
    !V = 1._dp/32*(x-L/2)**2
  end subroutine
    
  subroutine init_ops(opp_d,opp_u,opm,V,dt,dx,M)
    complex(dp), intent(inout) :: opp_d(:), opp_u(:), opm(:,:)
    real(dp), intent(in)       :: V(:), dt, dx
    integer, intent(in)        :: M

    integer :: i

    opm = (0._dp,0._dp)
    ! construct operators 
    do i = 1,M
      opm(i,i) = cmplx(1._dp,((-dt*2._dp/(dx**2) - V(i))/2._dp), dp)
      opp_d(i) = cmplx(1._dp, (dt*2._dp/(dx**2) + V(i))/2._dp, dp)

      if (i>1) then
        opm(i,i-1) = cmplx(0._dp, dt*0.5_dp/(dx**2), dp)
      endif

      if (i<M) then
        opm(i,i+1) = cmplx(0._dp, dt*0.5_dp/(dx**2), dp)
        opp_u(i) = cmplx(0._dp, -dt*0.5_dp/(dx**2), dp)
      endif
    enddo
  end subroutine
end module 
