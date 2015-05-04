module initialize
  use constants
  implicit none
  private
  public :: init_param, init_wavef, init_ops

contains
  subroutine init_param(dx, dt, L, M, n)
    real(dp), intent(out) :: dx, dt, L
    integer, intent(out)  :: M, n

    ! simulation parameters
    dx = 0.01_dp
    dt = 0.025_dp
    L = 60._dp
    M = floor(L/dx)
    n = 5000
  end subroutine

  subroutine init_wavef(psi, x, dx, L, k, M)
    complex(dp), intent(inout) :: psi(:) 
    real(dp), intent(inout)    :: x(:)
    real(dp), intent(in)       :: dx, L, k
    integer, intent(in)        :: M
    integer :: i
    
    do i = 1,M
      x(i) = i*dx
    enddo

    ! gaussian wavepackets
    psi = exp(-0.5_dp*(x-L/2)**2)*exp(cmplx(0._dp,k*x,dp)) !+ &
    ! exp(-0.5_dp*(x-2*L/3)**2)*exp(cmplx(0._dp,-k*x,dp))

    ! normalize wavefunction
    psi = psi/sqrt(sum(abs(psi)**2*dx))
  end subroutine

  subroutine init_ops(Ax, dt, dx, M)
    complex(dp), intent(inout) :: Ax(:,:)
    real(dp), intent(in)       :: dt, dx
    integer, intent(in)        :: M

    integer :: i

    ! initialize matrix operators in band storage fmt
    do i = 1,M
      Ax(1,i) = cmplx(0._dp, -0.5_dp*dt/dx**2, dp)
      Ax(2,i) = cmplx(1._dp, dt/dx**2, dp)
      Ax(3,i) = cmplx(0._dp, -0.5_dp*dt/dx**2, dp)
    enddo
  end subroutine
end module 
