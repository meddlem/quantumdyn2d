module initialize
  use constants
  use structures 
  implicit none
  private
  public :: init_param, init_wavef, init_ops

contains
  subroutine init_param(Q)
    type(model_parameters) :: Q
    
    ! model parameters
    Q%dx = 0.05_dp
    Q%dt = 0.02_dp
    Q%Lx = 10._dp
    Q%Ly = 10._dp
    Q%Mx = floor(Q%Lx/Q%dx)
    Q%My = floor(Q%Ly/Q%dx)
    Q%N = 5000
  end subroutine
  
  subroutine init_wavef(psi, x, y, Q)
    complex(dp), intent(inout) :: psi(:,:) 
    real(dp), intent(inout)    :: x(:,:), y(:,:)
    type(model_parameters)     :: Q
    
    real(dp), allocatable :: r(:,:), Hxy(:,:)
    real(dp)              :: A
    integer               :: i, j

    allocate(r(Q%Mx,Q%My), Hxy(Q%Mx,Q%My))
    A = 1._dp
    
    ! create grid
    do i = 1,Q%Mx
      do j = 1,Q%My
        x(i,j) = i*Q%dx
        y(i,j) = j*Q%dx
      enddo
    enddo
    
    ! starting position for wavepacket
    if (Q%V_type == 1) then
      r = sqrt((x - Q%Lx/2)**2 + (y - Q%Ly/2)**2) 
    elseif (Q%V_type == 2) then
      r = sqrt((x - Q%Lx/4)**2 + (y - Q%Ly/2)**2) 
      A = 2._dp
    endif

    Hxy = (x - Q%Lx/2)*(y - Q%Ly/2)
    psi = exp(-0.5_dp*A*r**2)*exp(i_u*(Q%kx*x + Q%ky*y))

    ! normalize wavefunction
    psi = psi/sqrt(sum(abs(psi)**2*Q%dx**2))

    deallocate(r, Hxy)
  end subroutine

  subroutine init_ops(Ax, Ay, Q)
    complex(dp), intent(inout) :: Ax(:,:,:), Ay(:,:,:)
    type(model_parameters) :: Q

    integer :: i, j

    ! init ADI matrix operators, x-dir, band storage fmt
    do i = 1,Q%Mx
      do j = 1,Q%My
        Ax(1,i,j) = -0.5_dp*i_u*Q%dt/Q%dx**2
        Ax(2,i,j) = one + i_u*Q%dt/Q%dx**2
        Ax(3,i,j) = -0.5_dp*i_u*Q%dt/Q%dx**2
      enddo
    enddo
    
    ! init ADI matrix operators, y-dir, band storage fmt
    do i = 1,Q%My
      do j = 1,Q%Mx
        Ay(1,i,j) = -0.5_dp*i_u*Q%dt/Q%dx**2 
        Ay(2,i,j) = one + i_u*Q%dt/Q%dx**2
        Ay(3,i,j) = -0.5_dp*i_u*Q%dt/Q%dx**2 
      enddo
    enddo
  end subroutine
end module 
