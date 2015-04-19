module io
  use constants
  implicit none
  private
  public :: user_in
contains

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

  subroutine user_in(k,dx,dt,L,M,n)
    real(dp), intent(out) :: k, dx, dt, L
    integer, intent(out)  :: M, n
    
    ! set parameters
    dx = 0.05_dp
    dt = 0.0005_dp
    L = 3._dp
    M = floor(L/dx)
    n = 100000
  
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "k = " 
    read(*,*) k
    write(*,'(A)') "Running simulation..."
  end subroutine
end module
