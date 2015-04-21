module initialize
  use constants
  implicit none
  private
  public :: init_wavef, init_V, init_ops

contains
  
  subroutine init_wavef(psi, x, y, dx, L_x, L_y, k_x, k_y, M_x, M_y)
    complex(dp), intent(inout) :: psi(:,:) 
    real(dp), intent(inout)    :: x(:,:), y(:,:)
    real(dp), intent(in)       :: dx, L_x, L_y, k_x, k_y
    integer, intent(in)        :: M_x, M_y
    
    real(dp), allocatable :: r(:,:), H_xy(:,:)
    integer               :: i, j

    allocate(r(M_x,M_y),H_xy(M_x,M_y))
    
    ! create grid
    do i = 1,M_x
      do j = 1,M_y
        x(i,j) = i*dx
        y(i,j) = j*dx
      enddo
    enddo
    
    ! distance to grid center
    r = sqrt((x-L_x/2)**2 + (y-L_y/2)**2) 

    ! ISQW wavefunction
    !psi = cmplx(sin(3*pi*x/L)*sin(2*pi*y/L),0._dp,dp) * &
    !  exp(cmplx(0._dp,k_x*x+k_y*y,dp))

    ! gaussian wavepackets
    H_xy = (x-L_x/2)*(y-L_y/2)**1
    psi = H_xy*exp(-0.5_dp*r**2)*exp(cmplx(0._dp,k_x*x+k_y*y,dp))

    ! normalize wavefunction
    psi = psi/sqrt(sum(abs(psi)**2*dx**2))

    deallocate(r)
  end subroutine

  subroutine init_V(V, x, y, L_x, L_y)
    real(dp), intent(in)    :: x(:,:), y(:,:), L_x, L_y
    real(dp), intent(inout) :: V(:,:)
    
    ! block/scattering potential
    !V = 0._dp
    !where(9._dp<x .and. x<13._dp) V = 80._dp
    
    ! harmonic potential
    V = 1._dp*((x-L_x/2)**2 + (y-L_y/2)**2)
  end subroutine
    
  subroutine init_ops(A_x_d, A_y_d, A_x_u, A_y_u, A_x_conj, A_y_conj, V, dt,& 
      dx, M_x, M_y)
    complex(dp), intent(inout) :: A_x_d(:,:), A_y_d(:,:), A_x_u(:,:), &
                                  A_y_u(:,:), A_x_conj(:,:,:), A_y_conj(:,:,:)
    real(dp), intent(in)       :: V(:,:), dt, dx
    integer, intent(in)        :: M_x, M_y

    integer :: i, j

    A_x_conj = (0._dp,0._dp)
    A_y_conj = (0._dp,0._dp)

    ! construct matrix operators, x-dir
    do i = 1,M_x
      do j = 1,M_y
        A_x_d(i,j) = cmplx(1._dp, 0.5_dp*dt*(2._dp/dx**2 + 0.5_dp*V(i,j)), dp)

        if (i<M_x) then
          A_x_u(i,j) = cmplx(0._dp, -dt*0.5_dp/dx**2, dp)
        endif
      enddo
    enddo
    
    ! construct matrix operators, y-dir
    do i = 1,M_y
      do j = 1,M_x
        A_y_d(i,j) = cmplx(1._dp, 0.5_dp*dt*(2._dp/dx**2 + 0.5_dp*V(j,i)), dp)

        if (i<M_y) then
          A_y_u(i,j) = cmplx(0._dp, -dt*0.5_dp/dx**2, dp)
        endif
      enddo
    enddo

    ! construct conjugate matrix ops, band storage fmt

    do i = 1,M_x
      do j = 1,M_y
        A_x_conj(2,i,j) = conjg(A_x_d(i,j))
        
        if (i>1) then
          A_x_conj(3,i-1,j) = conjg(A_x_u(i-1,j))
        endif
        
        if (i<M_x) then
          A_x_conj(1,i+1,j) = conjg(A_x_u(i,j))
        endif
      enddo
    enddo
    
    do i = 1,M_y
      do j = 1,M_x
        A_y_conj(2,i,j) = conjg(A_y_d(i,j))
        
        if (i>1) then
          A_y_conj(3,i-1,j) = conjg(A_y_u(i-1,j))
        endif
        
        if (i<M_y) then
          A_y_conj(1,i+1,j) = conjg(A_y_u(i,j))
        endif
      enddo
    enddo

  end subroutine
end module 
