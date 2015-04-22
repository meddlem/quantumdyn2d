module initialize
  use constants
  implicit none
  private
  public :: init_wavef, init_V, init_ops

contains
  
  subroutine init_wavef(psi,x,dx,L,k,M)
    complex(dp), intent(inout) :: psi(:) 
    real(dp), intent(inout)    :: x(:)
    real(dp), intent(in)       :: dx, L, k
    integer, intent(in)        :: M
    integer :: i
    
    do i = 1,M
      x(i) = i*dx
    enddo

    ! isq well
    ! psi_0 = sin(2*pi*x/L)*exp(cmplx(0._dp,k*x,dp))
    
    ! gaussian wavepackets
    psi = exp(-0.5_dp*(x-L/2)**2)*exp(cmplx(0._dp,k*x,dp)) !+ &
 !     exp(-0.5_dp*(x-2*L/3)**2)*exp(cmplx(0._dp,-k*x,dp))

    ! normalize wavefunction
    psi = psi/sqrt(sum(abs(psi)**2*dx))
  end subroutine

  subroutine init_V(V,x,L)
    real(dp), intent(in) :: x(:), L
    real(dp), intent(inout) :: V(:)
    
    ! block/scattering potential
    !V = 0._dp
    !where(28._dp<x .and. x<32._dp) V = 1._dp
    
    ! harmonic potential
    V = 1._dp*(x-L/2)**2
  end subroutine
    
  subroutine init_ops(A_x,V,dt,dx,M)
    complex(dp), intent(inout) :: A_x(:,:)
    real(dp), intent(in)       :: V(:), dt, dx
    integer, intent(in)        :: M

    integer :: i

    ! construct matrix operators, band storage fmt
    do i = 1,M
      A_x(1,i) = cmplx(0._dp, -0.5_dp*dt/dx**2, dp)
      A_x(2,i) = cmplx(1._dp, 0.5_dp*dt*(2._dp/dx**2 + V(i)), dp)
      A_x(3,i) = cmplx(0._dp, -0.5_dp*dt/dx**2, dp)
    enddo
  end subroutine
end module 
