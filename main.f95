program main
  use constants
  use structures 
  use initialize
  use simulation
  use io
  implicit none

  complex(dp), allocatable :: psi(:,:)
  real(dp), allocatable    :: x(:,:), y(:,:)
  type(Ops)                :: O
  type(modl_par)           :: Q

  ! initialize model parameters
  call init_param(Q)
  
  ! allocate arrays
  allocate(psi(Q%Mx,Q%My), x(Q%Mx,Q%My), y(Q%Mx,Q%My), O%Ax(3,Q%Mx,Q%My), &
    O%Ay(3,Q%My,Q%Mx))
  
  ! user input
  call get_usr_args(Q)
  call user_in(Q)
  
  ! initialize model 
  call init_ops(O, Q) 
  call init_wavef(psi, x, y, Q)

  ! run simulation
  call run_sim(psi, x, y, O, Q)
  
  deallocate(psi, x, y, O%Ax, O%Ay)
end program
