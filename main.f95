program main
  use constants
  use structures 
  use initialize
  use simulation
  use io
  implicit none

  complex(dp), allocatable :: psi(:,:), Ax(:,:,:), Ay(:,:,:)
  real(dp), allocatable    :: x(:,:), y(:,:)
  type(modl_par)           :: Q

  ! initialize model parameters
  call init_param(Q)
  
  ! allocate arrays
  allocate(psi(Q%Mx,Q%My), x(Q%Mx,Q%My), y(Q%Mx,Q%My), Ax(3,Q%Mx,Q%My), &
    Ay(3,Q%My,Q%Mx))
  
  ! user input
  call get_usr_args(Q)
  call user_in(Q)
  
  ! initialize model 
  call init_ops(Ax, Ay, Q) 
  call init_wavef(psi, x, y, Q)

  ! run simulation
  call run_sim(psi, x, y, Ax, Ay, Q)
  
  deallocate(psi, x, y, Ax, Ay)
end program
