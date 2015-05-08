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

    ! defaults 
    Q%sim_type = 'har' ! Harmonic potential
    P%plot_re = .false. ! Plot density 

    ! check command line arguments
    do i=1,iargc()
      call getarg(i,arg)
      if (trim(arg) == '--hsq') then
        Q%sim_type = 'hsq' ! adiabatic change: harmonic to isqw
      endif
      if (trim(arg) == '--dsl') then
        Q%sim_type = 'dsl' ! scattering potential
      endif
      if (trim(arg) == '--hqa') then
        Q%sim_type = 'hqa' ! adiabatic: harmonic to quartic
      endif
      if (trim(arg) == '-r') then
        P%plot_re = .true. ! plot real part
      endif
    enddo
  end subroutine

  subroutine user_in(Q)
    type(modl_par), intent(inout) :: Q
    
    if (any(Q%sim_type == ['dsl', 'har'])) then
      write(*,'(/,A,/)') '************ Input *************' 
      write(*,'(A)',advance='no') "kx = " 
      read(*,*) Q%kx
      write(*,'(A)',advance='no') "ky = " 
      read(*,*) Q%ky
    else
      Q%kx = 0._dp
      Q%ky = 0._dp
    endif

    write(*,'(A)') "Running simulation..."
  end subroutine
end module
