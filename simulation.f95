module simulation
  use constants
  use structures
  use plotroutines
  implicit none
  private
  public :: time_evo

contains
  subroutine time_evo(psi, x, A, Q)
    complex(dp), intent(inout) :: psi(:)
    complex(dp), intent(in)    :: A(:,:)
    real(dp), intent(in)       :: x(:) 
    type(modl_par), intent(in) :: Q

    real(dp), allocatable :: V(:), V1(:), V2(:)
    integer  :: i

    allocate(V(Q%M), V1(Q%M), V2(Q%M))
    call animate_plot(Q)
    
    do i = 1,Q%N
      call potential(V1, x, i*Q%dt, Q)
      call potential(V2, x, (i+1)*Q%dt, Q)
      V = 0.5_dp*(V1 + V2)

      call solve_nxt(psi, V, A, Q)

      if (mod(i,Q%plot_interval) == 0) then
        call plot_wavef(psi, x, V, Q)
      endif
    enddo

    call close_plot()
    deallocate(V, V1, V2)
  end subroutine

  subroutine solve_nxt(psi, V, A, Q)
    complex(dp), intent(inout) :: psi(:)
    real(dp), intent(inout)    :: V(:)
    complex(dp), intent(in)    :: A(:,:)
    type(modl_par), intent(in) :: Q

    complex(dp), allocatable :: g(:), A_tmp(:,:)
    integer                  :: info

    allocate(A_tmp(3,Q%M), g(Q%M))
    
    ! init temp arrays
    A_tmp = A
    A_tmp(2,:) = A_tmp(2,:) + cmplx(0._dp, 0.5_dp*Q%dt*V, dp)

    ! explicit part of calculation, mat-vec multiplication
    call zgbmv('N', Q%M, Q%M, 1, 1, one, conjg(A_tmp), 3, psi, 1, zero, g, 1)

    ! solve for wavefunction at t=n+1
    call zgtsv(Q%M, 1, A_tmp(1,1:Q%M-1), A_tmp(2,:), A_tmp(3,1:Q%M-1), g, Q%M, info)

    ! collect wavefunction at t=n+1
    psi = g

    deallocate(A_tmp, g)
  end subroutine

  pure subroutine potential(V, x, t, Q)
    real(dp), intent(inout)    :: V(:)
    real(dp), intent(in)       :: x(:), t
    type(modl_par), intent(in) :: Q

    real(dp) :: a
    
    ! changing harmonic potential
    a = 0.05_dp
    V = (1._dp - 0.8_dp*sin(a*t))*(x-Q%L/2)**2 
  end subroutine
end module
