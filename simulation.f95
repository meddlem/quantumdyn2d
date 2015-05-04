module simulation
  use constants
  use plotroutines
  use omp_lib
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(psi, V, x, y, n, Mx, My, Lx, Ly, Ax, Ay, dt, plot_re)
    complex(dp), intent(inout) :: psi(:,:)
    real(dp), intent(inout)    :: V(:,:)
    complex(dp), intent(in)    :: Ax(:,:,:), Ay(:,:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:), dt, Lx, Ly
    integer, intent(in)        :: n, Mx, My
    logical, intent(in)        :: plot_re

    integer  :: i

    do i=1,n
      call solve_nxt(psi, x, y, V, i*dt, dt, Mx, My, Lx, Ly, Ax, Ay)
      
      if (mod(i,plot_interval)==0) then
        call plot_wavef(psi, x, y, Mx, My, plot_re)
      endif
    enddo
  end subroutine

  subroutine solve_nxt(psi, x, y, V, t, dt, Mx, My, Lx, Ly, Ax, Ay)
    complex(dp), intent(inout) :: psi(:,:)
    real(dp), intent(inout)    :: V(:,:)
    complex(dp), intent(in)    :: Ax(:,:,:), Ay(:,:,:)
    real(dp), intent(in)       :: t, dt, Lx, Ly, x(:,:), y(:,:)
    integer, intent(in)        :: Mx, My

    complex(dp), allocatable :: gx(:), gy(:), Ax_d(:), Ax_l(:), &
                                Ax_u(:), Ay_d(:), Ay_l(:), &
                                Ay_u(:), Ax_tmp(:,:,:), Ay_tmp(:,:,:)
    integer                  :: i, info

    allocate(Ax_d(Mx), Ax_l(Mx-1), Ax_u(Mx-1), gx(Mx), Ay_d(My), &
      Ay_l(My-1), Ay_u(My-1), Ax_tmp(3,Mx,My), Ay_tmp(3,My,Mx), gy(My))
    
    ! calc potential at timestep t
    call Potential(V, x, y, t, Lx, Ly) 

    Ax_tmp = Ax
    Ay_tmp = Ay

    ! horizontal sweep
    !$omp parallel do private(Ax_d,Ax_u,Ax_l,gx)
    do i = 1,My
      ! init temp arrays
      gx = Ax_tmp(2,:,i)
      Ax_tmp(2,:,i) = Ax_tmp(2,:,i) + cmplx(0._dp, 0.25_dp*dt*V(:,i), dp)

      Ax_d = Ax_tmp(2,:,i) 
      Ax_u = Ax_tmp(1,1:Mx-1,i) 
      Ax_l = Ax_tmp(1,1:Mx-1,i)

      ! explicit part of calculation, mat-vec multiplication
      call zgbmv('N', Mx, Mx, 1, 1, one, conjg(Ax_tmp(:,:,i)), 3, psi(:,i), &
        1, zero, gx, 1)

      ! solve for psi at t=n+1/2
      call zgtsv(Mx, 1, Ax_l, Ax_d, Ax_u, gx, Mx, info)
      psi(:,i) = gx
    enddo
    !$omp end parallel do

    ! vertical sweep
    !$omp parallel do private(Ay_d,Ay_u,Ay_l,gy)
    do i = 1,Mx
      ! init temp arrays
      gy = Ay_tmp(2,:,i)
      Ay_tmp(2,:,i) = Ay_tmp(2,:,i) + cmplx(0._dp, 0.25_dp*dt*V(i,:), dp)

      Ay_u = Ay_tmp(1,1:My-1,i)
      Ay_d = Ay_tmp(2,:,i) 
      Ay_l = Ay_tmp(1,1:My-1,i)

      ! explicit part of calculation, mat-vec multiplication
      call zgbmv('N', My, My, 1, 1, one, conjg(Ay_tmp(:,:,i)), 3, psi(i,:), &
        1, zero, gy, 1)

      ! solve for psi at t=n+1
      call zgtsv(My, 1, Ay_l, Ay_d, Ay_u, gy, My, info)
      psi(i,:) = gy
    enddo
    !$omp end parallel do

    deallocate(Ax_tmp, Ay_tmp, Ax_d, Ax_l, Ax_u, Ay_d, Ay_l, Ay_u, gx, gy)
  end subroutine

  subroutine Potential(V, x, y, t, Lx, Ly)
    real(dp), intent(inout) :: V(:,:)
    real(dp), intent(in)    :: x(:,:), y(:,:), t, Lx, Ly
    
    ! scattering potential
    !V = 200._dp
    !where(Ly*0.4_dp<y .and. y<Ly*0.6_dp) V = 0._dp
    !where(Lx*0.49_dp>x .or. x>Lx*0.51_dp) V = 0._dp
    
    ! harmonic potential
    V = 1._dp*((x-Lx/2)**2 + (y-Ly/2)**2)
  end subroutine
end module
