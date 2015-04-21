program main
  use constants
  use initialize
  use evolve
  use plotroutines
  use io
  implicit none

  complex(dp), allocatable :: psi(:), opm(:,:), opp_d(:), opp_u(:)
  real(dp), allocatable    :: x(:), V(:)
  real(dp) :: k, dx, dt, L
  integer  :: M, n

  call user_in(k,dx,dt,L,M,n)
  allocate(psi(M),x(M),V(M),opm(M,M),opp_d(M),opp_u(M-1))
  
  call init_wavef(psi,x,dx,L,k,M)
  call init_V(V,x,L)
  call init_ops(opp_d,opp_u,opm,V,dt,dx,M)
  call animate_plot(L)
  call line_plot(x,abs(psi)**2,'x','P','','',1)

  call run_sim(psi,x,V,n,M,opp_d,opp_u,opm)
  
  call close_plot()
  deallocate(psi,x,V,opm,opp_d,opp_u)
end program
