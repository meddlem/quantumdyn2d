program main
  use constants
  use structures 
  use initialize
  use simulation
  use plotroutines
  use io
  implicit none

  complex(dp), allocatable :: psi(:,:), Ax(:,:,:), Ay(:,:,:)
  real(dp), allocatable    :: x(:,:), y(:,:), V(:,:)
  type(modl_par)           :: Q

  ! initialize model parameters
  call init_param(Q)
  
  ! allocate arrays
  allocate(psi(Q%Mx,Q%My), x(Q%Mx,Q%My), y(Q%Mx,Q%My), V(Q%Mx,Q%My), &
    Ax(3,Q%Mx,Q%My), Ay(3,Q%My,Q%Mx))
  
  ! initialize simulation
  call get_usr_args(Q)
  call user_in(Q)
  call init_ops(Ax, Ay, Q) 
  call init_wavef(psi, x, y, Q)
  call animate_plot(Q)

  call run_sim(psi, V, x, y, Ax, Ay, Q)
  
  call close_plot()

  deallocate(psi, x, y, V, Ax, Ay)
end program
