module initialize
  use constants
  use structures 
  implicit none
  private
  public :: init_param, init_wavef, init_ops

contains
  subroutine init_param(Q)
    type(modl_par), intent(inout) :: Q

    ! simulation parameters
    Q%dx = 0.01_dp
    Q%dt = 0.025_dp
    Q%L = 12._dp
    Q%M = floor(Q%L/Q%dx)
    Q%plot_interval = 100
    Q%tau = 400._dp
    Q%n = 20000
  end subroutine

  subroutine init_wavef(psi, x, Q)
    complex(dp), intent(inout) :: psi(:) 
    real(dp), intent(inout)    :: x(:)
    type(modl_par), intent(in) :: Q 

    integer :: i
    
    do i = 1,Q%M
      x(i) = i*Q%dx
    enddo

    ! gaussian wavepackets
    psi = exp(-0.5_dp*(x-Q%L/2)**2)*exp(cmplx(0._dp,Q%k*x,dp))

    ! normalize wavefunction
    psi = psi/sqrt(sum(abs(psi)**2*Q%dx))
  end subroutine

  subroutine init_ops(A, Q)
    complex(dp), intent(inout) :: A(:,:)
    type(modl_par), intent(in) :: Q

    real(dp) :: r

    r = Q%dt/Q%dx**2

    ! initialize matrix operator in band storage fmt
    A(1,:) = -0.5_dp*i_u*r
    A(2,:) = one + i_u*r
    A(3,:) = -0.5_dp*i_u*r
  end subroutine
end module 
