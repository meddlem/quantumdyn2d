module initialize
  use constants
  use structures 
  implicit none
  private
  public :: init_param, init_wavefunction, init_ops

contains
  subroutine init_param(Q, P)
    type(modl_par), intent(inout) :: Q
    type(plt_par), intent(inout)  :: P
    
    ! model parameters
    Q%dx = 0.05_dp
    Q%dt = 0.02_dp

    if (Q%V_type == 1) then
      Q%Lx = 8._dp
      Q%Ly = 8._dp
      
      P%plot_interval = 40
      P%rng = [-0.2_dp, 0.2_dp]
    elseif (Q%V_type == 2) then
      Q%Lx = 40._dp
      Q%Ly = 16._dp
      Q%Bx = Q%Lx/2
      Q%By = Q%Ly/2
      Q%Wx = Q%Lx*0.005_dp
      Q%Wy = Q%Ly*0.03_dp

      P%plot_interval = 5
      P%rng = [-1._dp, 0.2_dp]
    elseif (Q%V_type == 3) then
      Q%Lx = 12._dp
      Q%Ly = 12._dp
      
      P%plot_interval = 3
      P%rng = [-0.3_dp, 0.3_dp]
    endif
    
    Q%Mx = floor(Q%Lx/Q%dx)
    Q%My = floor(Q%Ly/Q%dx)
    Q%tau = 100._dp
    Q%N = 10000
  end subroutine
  
  subroutine init_wavefunction(psi, x, y, Q)
    complex(dp), intent(inout) :: psi(:,:) 
    real(dp), intent(inout)    :: x(:,:), y(:,:)
    type(modl_par), intent(in) :: Q
    
    real(dp), allocatable :: r(:,:), Hxy(:,:)
    real(dp)              :: A
    integer               :: i, j

    allocate(r(Q%Mx,Q%My), Hxy(Q%Mx,Q%My))
    
    ! create meshgrid
    do i = 1,Q%Mx
      do j = 1,Q%My
        x(i,j) = i*Q%dx
        y(i,j) = j*Q%dx
      enddo
    enddo
    
    if (Q%V_type == 1) then
      ! harmonic oscillator excited state
      r = sqrt((x - Q%Lx/2)**2 + (y - Q%Ly/2)**2) 
      Hxy = (x - Q%Lx/2)*(y - Q%Ly/2)
      A = 1._dp
    elseif (Q%V_type == 2) then
      ! gaussian wavepacket
      r = sqrt((x - Q%Lx/4)**2 + (y - Q%Ly/2)**2) 
      Hxy = 1._dp 
      A = 1._dp
    elseif (Q%V_type == 3) then
      ! harmonic oscillator excited state
      r = sqrt((x - Q%Lx/2)**2 + (y - Q%Ly/2)**2) 
      Hxy = (4*(x - Q%Lx/2)**2 - 2._dp)*(y - Q%Ly/2) 
      A = 1._dp
    endif

    ! calc wavefunction
    psi = Hxy*exp(-0.5_dp*A*r**2)*exp(i_u*(Q%kx*x + Q%ky*y))

    ! normalize wavefunction
    psi = psi/sqrt(sum(abs(psi)**2*Q%dx**2))

    deallocate(r, Hxy)
  end subroutine

  subroutine init_ops(Ax, Ay, Q)
    complex(dp), intent(inout) :: Ax(:,:), Ay(:,:,:)
    type(modl_par), intent(in) :: Q

    real(dp) :: r

    r = Q%dt/Q%dx**2

    ! init ADI matrix operators, x-dir, band storage fmt
    Ax(1,:) = -0.5_dp*i_u*r
    Ax(2,:) = one + i_u*r
    Ax(3,:) = -0.5_dp*i_u*r

    ! init ADI matrix operators, y-dir, band storage fmt
    Ay(1,:,:) = -0.5_dp*i_u*r
    Ay(2,:,:) = one + i_u*r
    Ay(3,:,:) = -0.5_dp*i_u*r
  end subroutine
end module 
