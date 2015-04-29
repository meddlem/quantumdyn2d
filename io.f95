module io
  use constants
  implicit none
  private
  public :: user_in
contains

  subroutine user_in(k)
    real(dp), intent(out) :: k
  
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "k = " 
    read(*,*) k
    write(*,'(A)') "Running simulation..."
  end subroutine
end module
