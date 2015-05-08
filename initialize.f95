module initialize
  use constants
  use structures 
  implicit none
  private
  public :: init_param, init_wavefunction, init_ops

contains
  subroutine init_param(Q, P)
    ! define model parameters
    type(modl_par), intent(inout) :: Q
    type(plt_par), intent(inout)  :: P
    
    ! spatial and time discretization steps
    Q%dx = 0.05_dp
    Q%dt = 0.02_dp
    
    ! time constant for adiabatic change of potential
    Q%tau = 100._dp
    
    ! number of iterations
    Q%N = 6000
    
    ! lattice length in x and y dirs 
    Q%Lx = 9._dp
    Q%Ly = 9._dp

    ! number of iterations between plot 
    P%plot_interval = 10

    ! min and max z value (colormap)
    P%rng = [-0.2_dp, 0.2_dp]

    if (any(Q%sim_type == ['hsq', 'hqa'])) then
      P%plot_interval = 20
    
    elseif (Q%sim_type == 'dsl') then
      Q%Lx = 30._dp
      Q%Ly = 13._dp

      Q%Bx = Q%Lx/2
      Q%By = Q%Ly/2
      Q%Wx = Q%Lx*0.005_dp
      Q%Wy = Q%Ly*0.03_dp
    endif
    
    ! number of points in x and y dir
    Q%Mx = floor(Q%Lx/Q%dx)
    Q%My = floor(Q%Ly/Q%dx)
  end subroutine
  
  subroutine init_wavefunction(psi, x, y, Q)
    complex(dp), intent(inout) :: psi(:,:) 
    real(dp), intent(inout)    :: x(:,:), y(:,:)
    type(modl_par), intent(in) :: Q
    
    real(dp), allocatable :: r(:,:), Hxy(:,:)
    real(dp)              :: sigma
    integer               :: i, j

    allocate(r(Q%Mx,Q%My), Hxy(Q%Mx,Q%My))
    
    ! create meshgrid
    do i = 1,Q%Mx
      do j = 1,Q%My
        x(i,j) = i*Q%dx
        y(i,j) = j*Q%dx
      enddo
    enddo
    
    sigma = 1._dp
    Hxy = 1._dp
    r = sqrt((x - Q%Lx/2)**2 + (y - Q%Ly/2)**2) 
    
    if (any(Q%sim_type == ['hsq', 'hqa'])) then
      ! harmonic oscillator excited state
      Hxy = (x - Q%Lx/2)*(y - Q%Ly/2)

    elseif (Q%sim_type == 'dsl') then
      ! gaussian wavepacket, starting off-center
      r = sqrt((x - Q%Lx/4)**2 + (y - Q%Ly/2)**2) 
      sigma = min(Q%Lx,Q%Ly)/20._dp
    
    elseif (Q%sim_type == 'har') then
      ! harmonic oscillator excited state
      Hxy = (4*(x - Q%Lx/2)**2 - 2._dp)*(y - Q%Ly/2) 
    endif

    ! define wavefunction
    psi = Hxy*exp(-0.5_dp*r**2/sigma**2)*exp(i_u*(Q%kx*x + Q%ky*y))

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
