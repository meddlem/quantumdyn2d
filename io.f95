module io
  use constants
  implicit none
  private
  public :: user_in
contains

  subroutine user_in(k_x,k_y,dx,dt,L,M,n)
    real(dp), intent(out) :: k_x, k_y, dx, dt, L
    integer, intent(out)  :: M, n
    
    ! set parameters
    dx = 0.15_dp
    dt = 0.1_dp
    L = 15._dp
    M = floor(L/dx)
    n = 5000
  
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "k_x = " 
    read(*,*) k_x
    write(*,'(A)',advance='no') "k_y = " 
    read(*,*) k_y
    write(*,'(A)') "Running simulation..."
  end subroutine

!  subroutine get_usr_args(use_T)
!    logical, intent(out) :: use_T
!    
!    character(10) :: arg
!    integer       :: i
!
!    use_T = .false. ! by default use betaJ for input
!    
!    do i=1,iargc()
!      call getarg(i,arg)
!      if (trim(arg) == '-T') then
!        use_T = .true.
!      endif
!    enddo
!  end subroutine
end module
