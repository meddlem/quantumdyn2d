program main
  use constants
  use initialize
  use evolve
  use plotroutines
!  use io
  implicit none

  complex(dp), allocatable :: psi(:), opm(:,:), opp_d(:), opp_du(:), opp_dl(:)
  real(dp), allocatable    :: x(:), V(:)
  real(dp) :: k, dx, dt, L
  integer  :: M, n

  ! set parameters
  k = 1._dp
  dx = 0.01_dp
  dt = 0.1_dp
  L = 3._dp
  M = floor(L/dx)
  n = 10

  !call user_in(k,dx,dt,L,M)
  allocate(psi(M),x(M),V(M),opm(M,M),opp_d(M),opp_dl(M-1),opp_du(M-1))
  call init_wavef(psi,x,dx,k,L,M)
  call init_V(V,x,L)
  call init_ops(opp_d,opp_dl,opp_du,opm,V,dt,dx,M)
  call line_plot(x,abs(psi)**2,'x','P','','',1)
!  call animate_plot()

  call run_sim(psi,x,L,n,M,opp_d,opp_dl,opp_du,opm)
  
!  call close_plot()
  deallocate(psi,x,V,opm,opp_d,opp_dl,opp_du)
end program
