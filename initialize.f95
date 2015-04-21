module initialize
  use constants
  implicit none
  private
  public :: init_wavef, init_V, init_ops

contains
  
  subroutine init_wavef(psi,x,y,dx,L,k_x,k_y,M)
    complex(dp), intent(inout) :: psi(:,:) 
    real(dp), intent(inout)    :: x(:,:), y(:,:)
    real(dp), intent(in)       :: dx, L, k_x, k_y
    integer, intent(in)        :: M
    integer :: i, j
    
    ! create grid
    do i = 1,M
      do j = 1,M
        x(i,j) = i*dx
        y(i,j) = j*dx
      enddo
    enddo

    ! ISQW wavefunction

    psi = cmplx(sin(pi*x/L)*sin(pi*y/L),0._dp,dp) * &
      exp(cmplx(0._dp,k_x*x+k_y*y,dp))

    ! gaussian wavepackets
    !psi = exp(-0.5_dp*((x-L/2)**2 + (y-L/2)**2)) * &
    !  exp(cmplx(0._dp,k_x*x+k_y*y,dp))

    ! normalize wavefunction
    psi = psi/sqrt(sum(abs(psi)**2*dx**2))
  end subroutine

  subroutine init_V(V,x,y,L)
    real(dp), intent(in) :: x(:,:), y(:,:), L
    real(dp), intent(inout) :: V(:,:)
    
    ! block/scattering potential
    V = 0._dp
    !where(28._dp<x .and. x<32._dp) V = 1._dp
    
    ! harmonic potential
    ! V = 1._dp/4*(x-L/2)**2
  end subroutine
    
  subroutine init_ops(A_d,A_u,A_conj,V,dt,dx,M)
    complex(dp), intent(inout) :: A_d(:), A_u(:), A_conj(:,:)
    real(dp), intent(in)       :: V(:,:), dt, dx
    integer, intent(in)        :: M

    integer :: i

    ! construct matrix operators
    A_conj = (0._dp,0._dp)

    do i = 1,M
      A_conj(i,i) = cmplx(1._dp, 0.5_dp*dt*(-2._dp/(dx**2)), dp)
      A_d(i) = cmplx(1._dp, 0.5_dp*dt*2._dp/(dx**2), dp)

      if (i>1) then
        A_conj(i,i-1) = cmplx(0._dp, 0.5_dp*dt/(dx**2), dp)
      endif

      if (i<M) then
        A_conj(i,i+1) = cmplx(0._dp, dt*0.5_dp/(dx**2), dp)
        A_u(i) = cmplx(0._dp, -dt*0.5_dp/(dx**2), dp)
      endif
    enddo
  end subroutine
end module 
