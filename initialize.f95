module initialize
  use constants
  implicit none
  private
  public :: init_wavef, init_V, init_ops

contains
  
  subroutine init_wavef(psi,x,y,dx,L,k,M)
    complex(dp), intent(inout) :: psi(:,:) 
    real(dp), intent(inout)    :: x(:,:), y(:,:)
    real(dp), intent(in)       :: dx, L, k
    integer, intent(in)        :: M
    integer :: i, j
    
    ! create grid
    do i = 1,M
      do j = 1,M
        x(i,j) = i*dx
        y(i,j) = j*dx
      enddo
    enddo

    ! gaussian wavepackets
    psi = exp(-0.5_dp*((x-L/2)**2 + (y-L/2)**2))*exp(cmplx(0._dp,k_x*x+k_y*y,dp)) 

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
    
  subroutine init_ops(opp_d,opp_u,opm,V,dt,dx,M)
    complex(dp), intent(inout) :: opp_d(:), opp_u(:), opm(:,:)
    real(dp), intent(in)       :: V(:), dt, dx
    integer, intent(in)        :: M

    integer :: i

    ! construct matrix operators, x-dir
    A_c = (0._dp,0._dp)

    do i = 1,M
      A_c(i,i) = cmplx(1._dp, 0.5_dp*dt*(-2._dp/(dx**2) - 0.5_dp*V(i,j)), dp)
      A_d(i) = cmplx(1._dp, 0.5_dp*dt*(2._dp/(dx**2) + V(i)), dp)

      if (i>1) then
        A_c(i,i-1) = cmplx(0._dp, 0.5_dp*dt/(dx**2), dp)
      endif

      if (i<M) then
        A_c(i,i+1) = cmplx(0._dp, dt*0.5_dp/(dx**2), dp)
        A_u(i) = cmplx(0._dp, -dt*0.5_dp/(dx**2), dp)
      endif
    enddo
    ! y dir
    B_c = (0._dp,0._dp)

  end subroutine
end module 
