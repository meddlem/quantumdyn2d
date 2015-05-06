module simulation
  use constants
  use structures 
  use plotroutines
  use omp_lib
  implicit none
  private
  public :: time_evo

contains
  subroutine time_evo(psi, x, y, Ax, Ay, Q)
    complex(dp), intent(inout) :: psi(:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:)
    complex(dp), intent(in)    :: Ax(:,:,:), Ay(:,:,:)
    type(modl_par), intent(in) :: Q

    real(dp), allocatable :: V(:,:), V1(:,:), V2(:,:)
    integer               :: i

    allocate(V(Q%Mx,Q%My), V1(Q%Mx,Q%My), V2(Q%Mx,Q%My))
    call animate_plot(Q)

    do i = 1,Q%N
      ! calculate potential
      call potential(V1, x, y, i*Q%dt, Q) 
      call potential(V2, x, y, (i+1)*Q%dt, Q) 
      V = 0.5_dp*(V1 + V2)

      ! time integration
      call solve_nxt(psi, V, Ax, Ay, Q)
      
      if (mod(i,Q%plot_interval) == 0) then
        call plot_wavef(psi, x, y, Q)
      endif
    enddo
    
    call close_plot()
    deallocate(V, V1, V2)
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
    Ay_tmp = conjg(Ay)
    Ay_tmp(2,:,:) = Ay_tmp(2,:,:) - 0.5_dp*i_u*Q%dt*transpose(V)

    ! sweep lattice in x and y direction    
    call x_sweep(psi, Ax_tmp, Ay_tmp, Q)
    
    ! init temp arrays
    Ax_tmp = conjg(Ax) 
    Ay_tmp = Ay
    Ay_tmp(2,:,:) = Ay(2,:,:) + 0.5_dp*i_u*Q%dt*transpose(V)

    call y_sweep(psi, Ax_tmp, Ay_tmp, Q)

    deallocate(Ax_tmp, Ay_tmp)
  end subroutine

  subroutine x_sweep(psi, Ax, Ay, Q)
    complex(dp), intent(inout) :: psi(:,:), Ax(:,:,:), Ay(:,:,:)
    type(modl_par), intent(in) :: Q

    complex(dp), allocatable :: gx(:,:)
    integer                  :: i, info
    
    allocate(gx(Q%Mx,Q%My))
    gx = psi 

    do i = 1,Q%Mx
      ! explicit part of calculation, mat-vec multiplication
      call zgbmv('N', Q%My, Q%My, 1, 1, one, Ay(:,:,i), 3, &
        psi(i,:), 1, zero, gx(i,:), 1)
    enddo

    do i = 1,Q%My
      ! solve resulting tridiagonal system for psi at t=n+1/2
      call zgtsv(Q%Mx, 1, Ax(1,1:Q%Mx-1,i), Ax(2,:,i), Ax(3,1:Q%Mx-1,i), &
        gx(:,i), Q%Mx, info)
    enddo
    
    psi = gx
    deallocate(gx)
  end subroutine

  subroutine y_sweep(psi, Ax, Ay, Q)
    complex(dp), intent(inout) :: psi(:,:), Ax(:,:,:), Ay(:,:,:)
    type(modl_par), intent(in) :: Q

    complex(dp), allocatable :: gy(:,:)
    integer                  :: i, info
    
    allocate(gy(Q%Mx,Q%My))
    gy = psi 

    do i = 1,Q%My
      ! explicit part of calculation, mat-vec multiplication
      call zgbmv('N', Q%Mx, Q%Mx, 1, 1, one, Ax(:,:,i), 3, &
        psi(:,i), 1, zero, gy(:,i), 1)
    enddo

    do i = 1,Q%Mx
      ! solve tridiagonal system for psi at t=n+1
      call zgtsv(Q%My, 1, Ay(1,1:Q%My-1,i), Ay(2,:,i), Ay(3,1:Q%My-1,i), &
        gy(i,:), Q%My, info)
    enddo
    
    psi = gy
    deallocate(gy)
  end subroutine

  pure subroutine potential(V, x, y, t, Q)
    real(dp), intent(inout)    :: V(:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:), t
    type(modl_par), intent(in) :: Q

    if (Q%V_type == 1) then
      ! adiabatic harmonic potential
      V = (x - Q%Lx/2)**2 + (y - Q%Ly/2)**2

    elseif (Q%V_type == 2) then
      ! constant scattering potential: single slit aperture
      !V = 10*(Q%kx**2+Q%ky**2)
      V = 0._dp

      where(Q%Ly*0.40_dp<y .and. y<Q%Ly*0.60_dp) V = 0._dp
      where(Q%Lx*0.49_dp>x .or. x>Q%Lx*0.51_dp) V = 0._dp
    endif
  end subroutine
end module
