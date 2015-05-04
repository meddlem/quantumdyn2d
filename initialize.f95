module initialize
  use constants
  implicit none
  private
  public :: init_param, init_wavef, init_ops

contains
  subroutine init_param(dx,dt,Lx,Ly,Mx,My,n)
    real(dp), intent(out) :: dx, dt, Lx, Ly
    integer, intent(out)  :: Mx, My, n
    
    ! model parameters
    dx = 0.1_dp
    dt = 0.025_dp
    Lx = 30._dp
    Ly = 10._dp
    Mx = floor(Lx/dx)
    My = floor(Ly/dx)
    n = 5000
  end subroutine
  
  subroutine init_wavef(psi, x, y, dx, Lx, Ly, kx, ky, Mx, My, V_type)
    complex(dp), intent(inout) :: psi(:,:) 
    real(dp), intent(inout)    :: x(:,:), y(:,:)
    real(dp), intent(in)       :: dx, Lx, Ly, kx, ky
    integer, intent(in)        :: Mx, My, V_type
    
    real(dp), allocatable :: r(:,:), Hxy(:,:)
    real(dp)              :: A
    integer               :: i, j

    allocate(r(Mx,My), Hxy(Mx,My))
    A = 1._dp
    
    ! create grid
    do i = 1,Mx
      do j = 1,My
        x(i,j) = i*dx
        y(i,j) = j*dx
      enddo
    enddo
    
    ! starting position for wavepacket
    if (V_type == 1) then
      r = sqrt((x - Lx/2)**2 + (y - Ly/2)**2) 
    elseif (V_type == 2) then
      r = sqrt((x - Lx/4)**2 + (y - Ly/2)**2) 
      A = 2._dp
    endif

    Hxy = (x - Lx/2)*(y - Ly/2)
    psi = exp(-0.5_dp*A*r**2)*exp(cmplx(0._dp,kx*x + ky*y,dp))

    ! normalize wavefunction
    psi = psi/sqrt(sum(abs(psi)**2*dx**2))

    deallocate(r, Hxy)
  end subroutine

  subroutine init_ops(Ax, Ay, dt, dx, Mx, My)
    complex(dp), intent(inout) :: Ax(:,:,:), Ay(:,:,:)
    real(dp), intent(in)       :: dt, dx
    integer, intent(in)        :: Mx, My

    integer :: i, j

    ! init ADI matrix operators, x-dir, band storage fmt
    do i = 1,Mx
      do j = 1,My
        Ax(1,i,j) = cmplx(0._dp, -dt*0.5_dp/dx**2, dp)
        Ax(2,i,j) = cmplx(1._dp, dt/dx**2, dp)
        Ax(3,i,j) = cmplx(0._dp, -dt*0.5_dp/dx**2, dp)
      enddo
    enddo
    
    ! init ADI matrix operators, y-dir, band storage fmt
    do i = 1,My
      do j = 1,Mx
        Ay(1,i,j) = cmplx(0._dp, -dt*0.5_dp/dx**2, dp)
        Ay(2,i,j) = cmplx(1._dp, dt/dx**2, dp)
        Ay(3,i,j) = cmplx(0._dp, -dt*0.5_dp/dx**2, dp)
      enddo
    enddo
  end subroutine
end module 
