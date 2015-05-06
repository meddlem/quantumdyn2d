module structures
  use constants
  implicit none

  type modl_par 
    ! data type that contains all model parameters 
    
    real(dp) :: k
    real(dp) :: L ! lattice length 
    real(dp) :: dx 
    real(dp) :: dt
    real(dp) :: a  ! time constant potential

    integer :: M ! lattice points in x-dir
    integer :: N  ! number of time steps/iterations
    integer :: V_type
    integer :: plot_interval

    logical :: plot_re
  end type
end module
