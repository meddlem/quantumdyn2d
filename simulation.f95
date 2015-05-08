module simulation
  use constants
  use structures 
  use plotroutines
  use omp_lib
  implicit none
  private
  public :: time_evo

contains
  subroutine time_evo(psi, x, y, Ax, Ay, Q, P)
    complex(dp), intent(inout) :: psi(:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:)
    complex(dp), intent(in)    :: Ax(:,:), Ay(:,:,:)
    type(modl_par), intent(in) :: Q
    type(plt_par), intent(in)  :: P

    real(dp), allocatable :: V(:,:), V1(:,:), V2(:,:)
    integer               :: i

    allocate(V(Q%Mx,Q%My), V1(Q%Mx,Q%My), V2(Q%Mx,Q%My))
    call animate_plot(Q, P)

    do i = 1,Q%N
      ! calculate potential using trapezoidal rule
      call potential(V1, x, y, i*Q%dt, Q) 
      call potential(V2, x, y, (i+1)*Q%dt, Q) 
      V = 0.5_dp*(V1 + V2)

      ! time integration
      call solve_nxt(psi, V, Ax, Ay, Q)
      
      if (mod(i,P%plot_interval) == 0) then
        call plot_wavef(psi, x, y, Q, P)
      endif
    enddo
    
    call close_plot()
    deallocate(V, V1, V2)
  end subroutine

  subroutine solve_nxt(psi, V, Ax, Ay, Q)
    complex(dp), intent(inout) :: psi(:,:)
    real(dp), intent(in)       :: V(:,:)
    complex(dp), intent(in)    :: Ax(:,:), Ay(:,:,:)
    type(modl_par), intent(in) :: Q

    complex(dp), allocatable :: Ax_tmp(:,:), Ay_tmp(:,:,:)

    ! init temp arrays
    allocate(Ax_tmp(3,Q%Mx), Ay_tmp(3,Q%My,Q%Mx))
    Ax_tmp = Ax 
    Ay_tmp = Ay

    ! add potential to Ay matrix  
    Ay_tmp(2,:,:) = Ay_tmp(2,:,:) + 0.5_dp*i_u*Q%dt*transpose(V)

    ! use Psi(n) to solve for intermediate Psi(n+1/2)
    call a_sweep(psi, Ax_tmp, Ay_tmp, Q)

    ! reset Ax matrix
    Ax_tmp = Ax 

    ! use Psi(n+1/2) to solve for Psi(n+1)
    call b_sweep(psi, Ax_tmp, Ay_tmp, Q)

    deallocate(Ax_tmp, Ay_tmp)
  end subroutine

  subroutine a_sweep(psi, Ax, Ay, Q)
    complex(dp), intent(inout) :: psi(:,:), Ax(:,:), Ay(:,:,:)
    type(modl_par), intent(in) :: Q

    complex(dp), allocatable :: g(:,:)
    integer                  :: i, info
    
    allocate(g(Q%Mx,Q%My))

    do i = 1,Q%Mx
      ! explicit part of calculation, mat-vec mult, using BLAS routine
      call zgbmv('N', Q%My, Q%My, 1, 1, one, conjg(Ay(:,:,i)), 3, &
        psi(i,:), 1, zero, g(i,:), 1)
    enddo

    ! solve tridiagonal system for psi at t=n+1/2, using LAPACK routine
    call zgtsv(Q%Mx, Q%My, Ax(1,1:Q%Mx-1), Ax(2,:), Ax(3,1:Q%Mx-1), &
      g, Q%Mx, info)
    
    psi = g
    deallocate(g)
  end subroutine

  subroutine b_sweep(psi, Ax, Ay, Q)
    complex(dp), intent(inout) :: psi(:,:), Ax(:,:), Ay(:,:,:)
    type(modl_par), intent(in) :: Q

    complex(dp), allocatable :: g(:,:)
    integer                  :: i, info
    
    allocate(g(Q%Mx,Q%My))

    do i = 1,Q%My
      ! explicit part of calculation, mat-vec mult, using BLAS routine
      call zgbmv('N', Q%Mx, Q%Mx, 1, 1, one, conjg(Ax), 3, &
        psi(:,i), 1, zero, g(:,i), 1)
    enddo

    !$omp parallel do
    do i = 1,Q%Mx
      ! solve tridiagonal system for psi at t=n+1, using LAPACK routine
      call zgtsv(Q%My, 1, Ay(1,1:Q%My-1,i), Ay(2,:,i), Ay(3,1:Q%My-1,i), &
        g(i,:), Q%My, info)
    enddo
    !$omp end parallel do
    
    psi = g
    deallocate(g)
  end subroutine

  pure subroutine potential(V, x, y, t, Q)
    real(dp), intent(inout)    :: V(:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:), t
    type(modl_par), intent(in) :: Q

    if (Q%sim_type == 'hsq') then
      ! adiabatic harmonic potential -> ISQW
      if (t < Q%tau) then
        V = (1._dp - t/Q%tau)**2*((x - Q%Lx/2)**2 + (y - Q%Ly/2)**2)
      else
        V = 0._dp
      endif

    elseif (Q%sim_type == 'hqa') then
      ! adiabatic harmonic potential -> quartic potential
      if (t < Q%tau) then
        V = (1._dp - t/Q%tau)**2*((x - Q%Lx/2)**2 + (y - Q%Ly/2)**2) + &
          t/Q%tau*((x - Q%Lx/2)**4 + (y - Q%Ly/2)**4)
      else
        V = (x - Q%Lx/2)**4 + (y - Q%Ly/2)**4
      endif

    elseif (Q%sim_type == 'dsl') then
      ! single slit aperture
      
      ! set barrier height
      V = 20_dp*(Q%kx**2 + Q%ky**2)
      
      ! set potential to zero outside of barrier
      where(abs(x-Q%Bx) > Q%Wx) V = 0._dp
      where(abs(y - (Q%By + 2_dp*Q%Wy)) < Q%Wy) V = 0._dp 
      where(abs(y - (Q%By - 2_dp*Q%Wy)) < Q%Wy) V = 0._dp

    elseif (Q%sim_type == 'har') then
      ! harmonic potential well
      V = (x - Q%Lx/2)**2 + (y - Q%Ly/2)**2
    endif
  end subroutine
end module
