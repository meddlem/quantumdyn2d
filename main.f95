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

      complex(dp), allocatable :: psi(:,:)
      real(dp), allocatable    :: x(:,:), y(:,:)
      type(Ops)                :: O
      
      allocate(psi(Q%Mx,Q%My), x(Q%Mx,Q%My), y(Q%Mx,Q%My), O%Ax(3,Q%Mx,Q%My), &
        O%Ay(3,Q%My,Q%Mx))
      
      call init_wavefunction(psi, x, y, Q)
      call init_ops(O, Q) 
      call time_evo(psi, x, y, O, Q)
      
      deallocate(psi, x, y, O%Ax, O%Ay)
    end subroutine
end program
