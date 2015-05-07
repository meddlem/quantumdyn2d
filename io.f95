module io
  use constants
  use structures
  implicit none
  private
  public :: user_in, get_usr_args
contains

  subroutine get_usr_args(Q, P)
    type(modl_par), intent(inout) :: Q
    type(plt_par), intent(inout)  :: P

    character(10) :: arg
    integer       :: i

    ! default 
    Q%V_type = 3 ! Harmonic potential
    P%plot_re = .false. ! Plot density 

    ! check command line arguments
    do i=1,iargc()
      call getarg(i,arg)
      if (trim(arg) == '-a') then
        Q%V_type = 1 ! adiabatic potential change 
      endif
      if (trim(arg) == '-s') then
        Q%V_type = 2 ! scattering potential
      endif
      if ((trim(arg) == '--PlotRe') .or. (trim(arg) == '-r')) then
        P%plot_re = .true. ! plot real part
      endif
    enddo
  end subroutine

  subroutine user_in(Q)
    type(modl_par), intent(inout) :: Q
    
    if (Q%V_type == 2) then
      write(*,'(/,A,/)') '************ Input *************' 
      write(*,'(A)',advance='no') "kx = " 
      read(*,*) Q%kx
      write(*,'(A)',advance='no') "ky = " 
      read(*,*) Q%ky
    endif
    write(*,'(A)') "Running simulation..."
  end subroutine
end module
