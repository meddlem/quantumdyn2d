module initialize
  use constants
  implicit none
  private
  public :: init_param, init_wavef, init_V, init_ops

contains
  subroutine init_param(dx,dt,Lx,Ly,Mx,My,n)
    real(dp), intent(out) :: dx, dt, Lx, Ly
    integer, intent(out)  :: Mx, My, n
    
    ! model parameters
    dx = 0.05_dp
    dt = 0.01_dp
    Lx = 30._dp
    Ly = 10._dp
    Mx = floor(Lx/dx)
    My = floor(Ly/dx)
    n = 5000
  end subroutine
  
  subroutine init_wavef(psi, x, y, dx, Lx, Ly, kx, ky, Mx, My)
    complex(dp), intent(inout) :: psi(:,:) 
    real(dp), intent(inout)    :: x(:,:), y(:,:)
    real(dp), intent(in)       :: dx, Lx, Ly, kx, ky
    integer, intent(in)        :: Mx, My
    
    real(dp), allocatable :: r(:,:), Hxy(:,:)
    integer               :: i, j

    allocate(r(Mx,My), Hxy(Mx,My))
    
    ! create grid
    do i = 1,Mx
      do j = 1,My
        x(i,j) = i*dx
        y(i,j) = j*dx
      enddo
    enddo
    
    ! starting position for wavepacket
    r = sqrt((x - Lx/4)**2 + (y - Ly/2)**2) 

    ! ISQW wavefunction
    !psi = cmplx(sin(3*pi*x/L)*sin(2*pi*y/L),0._dp,dp) * &
    !  exp(cmplx(0._dp,kx*x+ky*y,dp))

    ! gaussian wavepackets
    Hxy = (x - Lx/2)*(y - Ly/2)**1
    psi = exp(-0.5_dp*r**2)*exp(cmplx(0._dp,kx*x + ky*y,dp))

    ! normalize wavefunction
    psi = psi/sqrt(sum(abs(psi)**2*dx**2))

    deallocate(r, Hxy)
  end subroutine

  subroutine init_V(V, x, y, Lx, Ly)
    real(dp), intent(in)    :: x(:,:), y(:,:), Lx, Ly
    real(dp), intent(inout) :: V(:,:)
    
    ! scattering potential
    V = 80._dp
    where(Ly*0.4_dp<y .and. y<Ly*0.6_dp) V = 0._dp
    where(x>Lx/2) V = 0._dp
    
    ! harmonic potential
    ! V = 1._dp*((x-Lx/2)**2 + (y-Ly/2)**2)
  end subroutine
    
  subroutine init_ops(Ax, Ay, V, dt, dx, Mx, My)
    complex(dp), intent(inout) :: Ax(:,:,:), Ay(:,:,:)
    real(dp), intent(in)       :: V(:,:), dt, dx
    integer, intent(in)        :: Mx, My

    integer :: i, j

    ! construct ADI matrix operators, x-dir, band storage fmt
    do i = 1,Mx
      do j = 1,My
        Ax(1,i,j) = cmplx(0._dp, -dt*0.5_dp/dx**2, dp)
        Ax(2,i,j) = cmplx(1._dp, 0.5_dp*dt*(2._dp/dx**2 + 0.5_dp*V(i,j)), dp)
        Ax(3,i,j) = cmplx(0._dp, -dt*0.5_dp/dx**2, dp)
      enddo
    enddo
    
    ! construct ADI matrix operators, y-dir, band storage fmt
    do i = 1,My
      do j = 1,Mx
        Ay(1,i,j) = cmplx(0._dp, -dt*0.5_dp/dx**2, dp)
        Ay(2,i,j) = cmplx(1._dp, 0.5_dp*dt*(2._dp/dx**2 + 0.5_dp*V(j,i)), dp)
        Ay(3,i,j) = cmplx(0._dp, -dt*0.5_dp/dx**2, dp)
      enddo
    enddo
  end subroutine
end module 
