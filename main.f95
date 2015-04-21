program main
  use constants
  use initialize
  use evolve
  use plotroutines
  use io
  implicit none

  complex(dp), allocatable :: psi(:,:), A_x_conj(:,:,:), A_y_conj(:,:,:),&
                              A_x_d(:,:), A_y_d(:,:), A_x_u(:,:), A_y_u(:,:)
  real(dp), allocatable    :: x(:,:), y(:,:), V(:,:)
  real(dp) :: k_x, k_y, dx, dt, L_x, L_y
  integer  :: M_x, M_y, n

  call user_in(k_x,k_y,dx,dt,L_x,L_y,M_x,M_y,n)

  allocate(psi(M_x,M_y), x(M_x,M_y), y(M_x,M_y), V(M_x,M_y),&
    A_x_conj(M_x,M_x,M_y), A_y_conj(M_y,M_y,M_x), A_x_d(M_x,M_y),&
    A_y_d(M_y,M_x), A_x_u(M_x-1,M_y), A_y_u(M_y-1,M_x))
  
  ! initialize simulation
  call init_wavef(psi,x,y,dx,L_x,L_y,k_x,k_y,M_x,M_y)
  call init_V(V,x,y,L_x,L_y)
  call init_ops(A_x_d,A_y_d,A_x_u,A_y_u,A_x_conj,A_y_conj,V,dt,dx,M_x,M_y)
  call animate_plot(L_x,L_y)
  
  call run_sim(psi,x,y,n,M_x,M_y,A_x_d,A_y_d,A_x_u,A_y_u,A_x_conj,A_y_conj)
  
  call close_plot()

  deallocate(psi,x,y,V,A_x_conj,A_y_conj,A_x_d,A_y_d,A_x_u,A_y_u)
end program
