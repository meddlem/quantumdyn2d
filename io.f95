module io
  use constants
  implicit none
  private
  public :: user_in
contains

  subroutine user_in(k_x,k_y,dx,dt,L_x,L_y,M_x,M_y,n)
    real(dp), intent(out) :: k_x, k_y, dx, dt, L_x, L_y
    integer, intent(out)  :: M_x, M_y, n
    
    ! set parameters
    dx = 0.15_dp
    dt = 0.1_dp
    L_x = 15._dp
    L_y = 15._dp
    M_x = floor(L_x/dx)
    M_y = floor(L_y/dx)
    n = 5000
  
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "k_x = " 
    read(*,*) k_x
    write(*,'(A)',advance='no') "k_y = " 
    read(*,*) k_y
    write(*,'(A)') "Running simulation..."
  end subroutine
end module
