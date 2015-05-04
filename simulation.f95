module simulation
  use constants
  use structures 
  use plotroutines
  use omp_lib
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(psi, x, y, Ax, Ay, Q)
    complex(dp), intent(inout) :: psi(:,:)
    complex(dp), intent(in)    :: Ax(:,:,:), Ay(:,:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:)
    type(modl_par), intent(in) :: Q

    real(dp), allocatable :: V(:,:)
    integer               :: i

    allocate(V(Q%Mx,Q%My))
    call animate_plot(Q)

    do i = 1,Q%N
      call potential(V, x, y, i*Q%dt, Q) 
      call solve_nxt(psi, V, Ax, Ay, Q)
      
      if (mod(i,plot_interval) == 0) then
        call plot_wavef(psi, x, y, Q)
      endif
    enddo
    
    call close_plot()
    deallocate(V)
  end subroutine

  subroutine solve_nxt(psi, V, Ax, Ay, Q)
    complex(dp), intent(inout) :: psi(:,:)
    real(dp), intent(in)       :: V(:,:)
    complex(dp), intent(in)    :: Ax(:,:,:), Ay(:,:,:)
    type(modl_par), intent(in) :: Q

    complex(dp), allocatable :: Ax_tmp(:,:,:), Ay_tmp(:,:,:)

    allocate(Ax_tmp(3,Q%Mx,Q%My), Ay_tmp(3,Q%My,Q%Mx))
    
    ! init temp arrays
    Ax_tmp = Ax
    Ay_tmp = Ay
    
    call h_sweep(psi, V, Ax_tmp, Q)
    call v_sweep(psi, V, Ay_tmp, Q)

    deallocate(Ax_tmp, Ay_tmp)
  end subroutine

  subroutine h_sweep(psi, V, Ax_tmp, Q)
    complex(dp), intent(inout) :: psi(:,:), Ax_tmp(:,:,:)
    real(dp), intent(in)       :: V(:,:)
    type(modl_par), intent(in) :: Q

    complex(dp), allocatable :: gx(:), Ax_d(:), Ax_l(:), Ax_u(:)
    integer :: i, info
    
    allocate(Ax_d(Q%Mx), Ax_l(Q%Mx-1), Ax_u(Q%Mx-1), gx(Q%Mx))

    !$omp parallel do private(Ax_d,Ax_u,Ax_l,gx)
    do i = 1,Q%My
      ! init temp arrays
      Ax_tmp(2,:,i) = Ax_tmp(2,:,i) + 0.25_dp*i_u*Q%dt*V(:,i)

      gx = Ax_tmp(2,:,i)
      Ax_d = Ax_tmp(2,:,i) 
      Ax_u = Ax_tmp(1,1:Q%Mx-1,i) 
      Ax_l = Ax_tmp(1,1:Q%Mx-1,i)

      ! explicit part of calculation, mat-vec multiplication
      call zgbmv('N', Q%Mx, Q%Mx, 1, 1, one, conjg(Ax_tmp(:,:,i)), 3, &
        psi(:,i), 1, zero, gx, 1)

      ! solve resulting tridiagonal system for psi at t=n+1/2
      call zgtsv(Q%Mx, 1, Ax_l, Ax_d, Ax_u, gx, Q%Mx, info)
      psi(:,i) = gx
    enddo
    !$omp end parallel do

    deallocate(Ax_d, Ax_l, Ax_u, gx)
  end subroutine

  subroutine v_sweep(psi, V, Ay_tmp, Q)
    complex(dp), intent(inout) :: psi(:,:), Ay_tmp(:,:,:)
    real(dp), intent(in)       :: V(:,:)
    type(modl_par), intent(in) :: Q

    complex(dp), allocatable :: gy(:), Ay_d(:), Ay_l(:), Ay_u(:)
    integer :: i, info
    
    allocate(Ay_d(Q%My), Ay_l(Q%My-1), Ay_u(Q%My-1), gy(Q%My))

    !$omp parallel do private(Ay_d,Ay_u,Ay_l,gy)
    do i = 1,Q%Mx
      ! init temp arrays
      gy = Ay_tmp(2,:,i)
      Ay_tmp(2,:,i) = Ay_tmp(2,:,i) + 0.25_dp*i_u*Q%dt*V(i,:)

      Ay_u = Ay_tmp(1,1:Q%My-1,i)
      Ay_d = Ay_tmp(2,:,i) 
      Ay_l = Ay_tmp(1,1:Q%My-1,i)

      ! explicit part of calculation, mat-vec multiplication
      call zgbmv('N', Q%My, Q%My, 1, 1, one, conjg(Ay_tmp(:,:,i)), 3, &
        psi(i,:), 1, zero, gy, 1)

      ! solve tridiagonal system for psi at t=n+1
      call zgtsv(Q%My, 1, Ay_l, Ay_d, Ay_u, gy, Q%My, info)
      psi(i,:) = gy
    enddo
    !$omp end parallel do

    deallocate(Ay_d, Ay_l, Ay_u, gy)
  end subroutine

  pure subroutine potential(V, x, y, t, Q)
    real(dp), intent(inout) :: V(:,:)
    real(dp), intent(in)    :: x(:,:), y(:,:), t
    type(modl_par), intent(in) :: Q

    real(dp) :: a
    
    ! harmonic potential
    a = 0.1_dp
    V = ((1-0.8_dp*sin(a*t))*(x-Q%Lx/2)**2 + (y-Q%Ly/2)**2)
  end subroutine
end module
