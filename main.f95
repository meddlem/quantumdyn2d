program main
  use constants
  use initialize
  use evolve
  use plotroutines
  use io
  implicit none

  complex(dp), allocatable :: psi(:,:), A_x(:,:,:), A_y(:,:,:)
  real(dp), allocatable    :: x(:,:), y(:,:), V(:,:)
  real(dp) :: k_x, k_y, dx, dt, L_x, L_y
  integer  :: M_x, M_y, n

  call user_in(k_x,k_y,dx,dt,L_x,L_y,M_x,M_y,n)

  allocate(psi(M_x,M_y), x(M_x,M_y), y(M_x,M_y), V(M_x,M_y),A_x(3,M_x,M_y), &
    A_y(3,M_y,M_x))
  
  ! initialize simulation
  call init_wavef(psi,x,y,dx,L_x,L_y,k_x,k_y,M_x,M_y)
  call init_V(V,x,y,L_x,L_y)
  call init_ops(A_x,A_y,V,dt,dx,M_x,M_y)
  call animate_plot(L_x,L_y)

  call run_sim(psi,x,y,n,M_x,M_y,A_x,A_y)
  
  call close_plot()

  deallocate(psi,x,y,V,A_x,A_y)
end program
