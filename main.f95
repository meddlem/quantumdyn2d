program main
  use constants
  use structures 
  use initialize
  use simulation
  use io
  implicit none

  type(modl_par) :: Q

  call get_usr_args(Q)
  call user_in(Q)
  call init_param(Q)
  call run_sim(Q)

  contains
    subroutine run_sim(Q)
      type(modl_par), intent(in) :: Q

      complex(dp), allocatable :: psi(:,:), Ax(:,:), Ay(:,:,:)
      real(dp), allocatable    :: x(:,:), y(:,:)
      
      allocate(psi(Q%Mx,Q%My), x(Q%Mx,Q%My), y(Q%Mx,Q%My), &
        Ax(3,Q%Mx), Ay(3,Q%My,Q%Mx))
      
      call init_wavefunction(psi, x, y, Q)
      call init_ops(Ax, Ay, Q) 
      call time_evo(psi, x, y, Ax, Ay, Q)
      
      deallocate(psi, x, y, Ax, Ay)
    end subroutine
end program
