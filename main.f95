program main
  use constants
  use initialize
  use simulation
  use plotroutines
  use io
  implicit none

  complex(dp), allocatable :: psi(:), A_x(:,:)
  real(dp), allocatable    :: x(:), V(:)
  real(dp) :: k, dx, dt, L
  integer  :: M, n

  call user_in(k,dx,dt,L,M,n)
  allocate(psi(M),x(M),V(M),A_x(3,M))
  
  call init_wavef(psi,x,dx,L,k,M)
  call init_V(V,x,L)
  call init_ops(A_x,V,dt,dx,M)
  call animate_plot(L)

  call run_sim(psi,x,V,n,M,A_x)
  
  call close_plot()
  deallocate(psi,x,V,A_x)
end program
