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

      complex(dp), allocatable :: psi(:), A(:,:)
      real(dp), allocatable    :: x(:) 
      
      ! allocate arrays
      allocate(psi(Q%M), x(Q%M), A(3,Q%M))
      
      ! init simulation
      call init_wavef(psi, x, Q)
      call init_ops(A, Q)

      ! time integration
      call time_evo(psi, x, A, Q)

      deallocate(psi, x, A)
    end subroutine
end program
