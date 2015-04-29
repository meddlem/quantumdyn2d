module io
  use constants
  implicit none
  private
  public :: user_in
contains

  subroutine user_in(kx,ky)
    real(dp), intent(out) :: kx, ky
    
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "kx = " 
    read(*,*) kx
    write(*,'(A)',advance='no') "ky = " 
    read(*,*) ky
    write(*,'(A)') "Running simulation..."
  end subroutine
end module
