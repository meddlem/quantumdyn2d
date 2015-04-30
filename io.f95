module io
  use constants
  implicit none
  private
  public :: user_in, get_usr_args
contains

  subroutine get_usr_args(V_type, plot_re)
    integer, intent(out) :: V_type 
    logical, intent(out) :: plot_re 

    character(10) :: arg
    integer       :: i

    ! default 
    V_type = 1 ! Harmonic potential
    plot_re = .false. ! Plot density 

    ! check command line arguments
    do i=1,iargc()
      call getarg(i,arg)
      if ((trim(arg) == '--S') .or. (trim(arg) == '-s')) then
        V_type = 2 ! scattering potential
      endif
      if ((trim(arg) == '--PlotRe') .or. (trim(arg) == '-r')) then
        plot_re = .true. ! plot real part
      endif
    enddo
  end subroutine

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
