module evolve
  use constants
  implicit none
  private
  public :: run_sim

contains

  subroutine run_sim(psi,x,V,L,M)
    
    do i=1,n
      call solve(psi,V,)
      call plot_wavef(psi,x)
    enddo

  end subroutine
end module
