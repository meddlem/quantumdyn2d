module io
  use constants
  use structures
  implicit none
  private
  public :: user_in, get_usr_args
contains

  subroutine get_usr_args(Q)
    type(modl_par), intent(inout) :: Q

    character(10) :: arg
    integer       :: i

    ! default 
    Q%V_type = 1 ! Harmonic potential
    Q%plot_re = .false. ! Plot density 

    ! check command line arguments
    do i=1,iargc()
      call getarg(i,arg)
      if ((trim(arg) == '--S') .or. (trim(arg) == '-s')) then
        Q%V_type = 2 ! scattering potential
      endif
      if ((trim(arg) == '--PlotRe') .or. (trim(arg) == '-r')) then
        Q%plot_re = .true. ! plot real part
      endif
    enddo
  end subroutine

  subroutine user_in(Q)
    type(modl_par), intent(inout) :: Q
    
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "kx = " 
    read(*,*) Q%kx
    write(*,'(A)',advance='no') "ky = " 
    read(*,*) Q%ky
    write(*,'(A)') "Running simulation..."
  end subroutine
end module
