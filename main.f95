program main
  use constants
  use initialize
  use simulation
  use plotroutines
  use io
  implicit none

  complex(dp), allocatable :: psi(:), Ax(:,:)
  real(dp), allocatable    :: x(:), V(:)
  real(dp) :: k, dx, dt, L
  integer  :: M, n

  ! init parameters
  call init_param(dx, dt, L, M, n)

  ! allocate arrays
  allocate(psi(M), x(M), V(M), Ax(3, M))
  
  ! init simulation
  call user_in(k)
  call init_wavef(psi, x, dx, L, k, M)
  call init_ops(Ax, dt, dx, M)
  call animate_plot(L)

  ! time integration
  call run_sim(psi, x, V, n, L, dt, M, Ax)

  ! close off 
  call close_plot()
  deallocate(psi, x, V, Ax)
end program
