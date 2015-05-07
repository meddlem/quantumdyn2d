program main
  use constants
  use structures 
  use initialize
  use simulation
  use io
  implicit none

  type(modl_par) :: Q
  type(plt_par)  :: P

  call get_usr_args(Q, P)
  call user_in(Q)
  call init_param(Q, P)
  call run_sim(Q, P)

  contains
    subroutine run_sim(Q, P)
      type(modl_par), intent(in) :: Q
      type(plt_par), intent(in)  :: P

      complex(dp), allocatable :: psi(:,:), Ax(:,:), Ay(:,:,:)
      real(dp), allocatable    :: x(:,:), y(:,:)
      
      allocate(psi(Q%Mx,Q%My), x(Q%Mx,Q%My), y(Q%Mx,Q%My), &
        Ax(3,Q%Mx), Ay(3,Q%My,Q%Mx))
      
      call init_wavefunction(psi, x, y, Q)
      call init_ops(Ax, Ay, Q) 
      call time_evo(psi, x, y, Ax, Ay, Q, P)
      
      deallocate(psi, x, y, Ax, Ay)
    end subroutine
end program
