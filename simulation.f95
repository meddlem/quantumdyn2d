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
    Ay_tmp = Ay
    Ay_tmp(2,:,:) = Ay_tmp(2,:,:) + 0.5_dp*i_u*Q%dt*transpose(V)

    ! first step, solve for intermediate Psi at t=n+1/2
    call sweep(psi, Ax_tmp, Ay_tmp, 1)

    Ax_tmp = Ax 

    ! second step, solve for Psi at t=n+1
    call sweep(psi, Ay_tmp, Ax_tmp, 2)

    deallocate(Ax_tmp, Ay_tmp)
  end subroutine

  subroutine sweep(psi, A, B, sdim)
    complex(dp), intent(inout) :: psi(:,:), A(:,:,:), B(:,:,:)
    integer, intent(in)        :: sdim

    complex(dp), allocatable :: g(:,:)
    integer                  :: i, n, m, info
    
    n = size(B,2)
    m = size(A,2)
    B = conjg(B)

    if (sdim == 1) then
      allocate(g(m,n))
      g = psi 
    else
      allocate(g(n,m))
      g = transpose(psi)
    endif

    !$omp parallel do 
    do i = 1,n
      ! mat-vec multiplication
      call zgbmv('N', n, n, 1, 1, one, B(:,:,i), 3, psi(i,:), 1, zero, &
        g(i,:), 1)
    enddo
    !$omp end parallel do
    
    B = conjg(B)

    !$omp parallel do 
    do i = 1,m
      ! solve resulting tridiagonal system
      call zgtsv(m, 1, A(1,1:m-1,i), A(2,:,i), A(3,1:m-1,i), g(:,i), m, info)
    enddo
    !$omp end parallel do
    
    if (sdim == 1) then
      psi = g
    else
      psi = transpose(g)
    endif

    deallocate(g)
  end subroutine

  pure subroutine potential(V, x, y, t, Q)
    real(dp), intent(inout)    :: V(:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:), t
    type(modl_par), intent(in) :: Q

    if (Q%V_type == 1) then
      ! adiabatic harmonic potential
      V = (1._dp - 0.5_dp*sin(Q%a*t))*(x - Q%Lx/2)**2 + &
        (1._dp - 0.7_dp*sin(Q%a*t))*(y - Q%Ly/2)**2

    elseif (Q%V_type == 2) then
      ! constant scattering potential: single slit aperture
      V = 10*(Q%kx**2+Q%ky**2)

      where(Q%Ly*0.40_dp<y .and. y<Q%Ly*0.60_dp) V = 0._dp
      where(Q%Lx*0.49_dp>x .or. x>Q%Lx*0.51_dp) V = 0._dp
    endif
  end subroutine
end module
