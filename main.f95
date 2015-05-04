program main
  use constants
  use initialize
  use simulation
  use plotroutines
  use io
  implicit none

  complex(dp), allocatable :: psi(:,:), Ax(:,:,:), Ay(:,:,:)
  real(dp), allocatable    :: x(:,:), y(:,:), V(:,:)
  real(dp) :: kx, ky, dx, dt, Lx, Ly
  integer  :: Mx, My, n, V_type 
  logical  :: plot_re

  ! initialize model parameters
  call init_param(dx, dt, Lx, Ly, Mx, My, n)
  
  ! allocate arrays
  allocate(psi(Mx,My), x(Mx,My), y(Mx,My), V(Mx,My), Ax(3,Mx,My), &
    Ay(3,My,Mx))
  
  ! initialize simulation
  call user_in(kx, ky)
  call get_usr_args(V_type, plot_re)
  call init_ops(Ax, Ay, dt, dx, Mx, My)
  call init_wavef(psi, x, y, dx, Lx, Ly, kx, ky, Mx, My, V_type)
  call animate_plot(Lx, Ly, plot_re)

  call run_sim(psi, V, x, y, n, Mx, My, Lx, Ly, Ax, Ay, dt, plot_re)
  
  call close_plot()

  deallocate(psi, x, y, V, Ax, Ay)
end program
