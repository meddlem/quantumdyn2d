module simulation
  use constants
  use structures 
  use plotroutines
  use omp_lib
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(psi, V, x, y, Ax, Ay, Q)
    complex(dp), intent(inout) :: psi(:,:)
    real(dp), intent(inout)    :: V(:,:)
    complex(dp), intent(in)    :: Ax(:,:,:), Ay(:,:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:)
    type(modl_par), intent(in) :: Q

    integer  :: i

    do i=1,Q%N
      call solve_nxt(psi, x, y, V, i*Q%dt, Ax, Ay, Q)
      
      if (mod(i,plot_interval) == 0) then
        call plot_wavef(psi, x, y, Q)
      endif
    enddo
  end subroutine

  subroutine solve_nxt(psi, x, y, V, t, Ax, Ay, Q)
    complex(dp), intent(inout) :: psi(:,:)
    real(dp), intent(inout)    :: V(:,:)
    complex(dp), intent(in)    :: Ax(:,:,:), Ay(:,:,:)
    real(dp), intent(in)       :: t, x(:,:), y(:,:)
    type(modl_par), intent(in) :: Q

    complex(dp), allocatable :: gx(:), gy(:), Ax_d(:), Ax_l(:), &
                                Ax_u(:), Ay_d(:), Ay_l(:), &
                                Ay_u(:), Ax_tmp(:,:,:), Ay_tmp(:,:,:)

    allocate(Ax_d(Q%Mx), Ax_l(Q%Mx-1), Ax_u(Q%Mx-1), gx(Q%Mx), Ay_d(Q%My), &
      Ay_l(Q%My-1), Ay_u(Q%My-1), Ax_tmp(3,Q%Mx,Q%My), Ay_tmp(3,Q%My,Q%Mx), &
      gy(Q%My))
    
    ! calc potential at timestep t
    call Potential(V, x, y, t, Q) 

    ! init temp arrays
    Ax_tmp = Ax
    Ay_tmp = Ay
    
    call h_sweep(psi, V, Ax_tmp, Ax_d, Ax_u, Ax_l, gx, Q)
    call v_sweep(psi, V, Ay_tmp, Ay_d, Ay_u, Ay_l, gy, Q)

    deallocate(Ax_tmp, Ay_tmp, Ax_d, Ax_l, Ax_u, Ay_d, Ay_l, Ay_u, gx, gy)
  end subroutine

  subroutine h_sweep(psi, V, Ax_tmp, Ax_d, Ax_u, Ax_l, gx, Q)
    complex(dp), intent(inout) :: psi(:,:), Ax_d(:), Ax_u(:), Ax_l(:), &
      Ax_tmp(:,:,:), gx(:)
    real(dp), intent(in)       :: V(:,:)
    type(modl_par), intent(in) :: Q

    integer :: i, info

    !$omp parallel do private(Ax_d,Ax_u,Ax_l,gx)
    do i = 1,Q%My
      ! init temp arrays
      gx = Ax_tmp(2,:,i)
      Ax_tmp(2,:,i) = Ax_tmp(2,:,i) + 0.25_dp*i_u*Q%dt*V(:,i)

      Ax_d = Ax_tmp(2,:,i) 
      Ax_u = Ax_tmp(1,1:Q%Mx-1,i) 
      Ax_l = Ax_tmp(1,1:Q%Mx-1,i)

      ! explicit part of calculation, mat-vec multiplication
      call zgbmv('N', Q%Mx, Q%Mx, 1, 1, one, conjg(Ax_tmp(:,:,i)), 3, &
        psi(:,i), 1, zero, gx, 1)

      ! solve for psi at t=n+1/2
      call zgtsv(Q%Mx, 1, Ax_l, Ax_d, Ax_u, gx, Q%Mx, info)
      psi(:,i) = gx
    enddo
    !$omp end parallel do
  end subroutine

  subroutine v_sweep(psi, V, Ay_tmp, Ay_d, Ay_u, Ay_l, gy, Q)
    complex(dp), intent(inout) :: psi(:,:), Ay_d(:), Ay_u(:), Ay_l(:), &
      Ay_tmp(:,:,:), gy(:)
    real(dp), intent(in)       :: V(:,:)
    type(modl_par), intent(in) :: Q

    integer :: i, info

    ! vertical sweep
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

      ! solve for psi at t=n+1
      call zgtsv(Q%My, 1, Ay_l, Ay_d, Ay_u, gy, Q%My, info)
      psi(i,:) = gy
    enddo
    !$omp end parallel do
  end subroutine

  pure subroutine Potential(V, x, y, t, Q)
    real(dp), intent(inout) :: V(:,:)
    real(dp), intent(in)    :: x(:,:), y(:,:), t
    type(modl_par), intent(in) :: Q

    real(dp) :: a
    
    ! harmonic potential
    a = 0.1_dp
    V = ((1-0.8_dp*sin(a*t))*(x-Q%Lx/2)**2 + (y-Q%Ly/2)**2)
  end subroutine
end module
