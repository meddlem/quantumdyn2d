module initialize
  use constants
  implicit none
  private
  public :: init_wavef, init_V

contains
  
subroutine init_wavef(psi_0,x,k,L,M)
  complex(dp), intent(inout) :: psi_0(:) 
  real(dp), intent(inout)    :: x(:)
  real(dp), intent(in)       :: k, L
  integer, intent(in)        :: M

  integer :: i
  
  do i = 1,M
    x = i*dx
  enddo

  psi_0 = 1/sqrt(L)*exp(cmplx(0._dp,k*x))
end subroutine

subroutine init_V(V,x,L)
  real(dp), intent(in) :: x(:), L
  real(dp), intent(inout) :: V(:)
  
  V = 0._dp
  ! block potential
  where(x>3*L/4) V = 1._dp

end subroutine

end module 
