module structures
  use constants
  implicit none
  
  type modl_par 
    ! data type that contains all model parameters 
    
    real(dp) :: kx 
    real(dp) :: ky
    real(dp) :: Lx ! lattice length in x-dir
    real(dp) :: Ly ! lattice length in y-dir
    real(dp) :: dx 
    real(dp) :: dt
    real(dp) :: tau ! time scale potential variation
    real(dp) :: Bx ! Barrier loc x-dir
    real(dp) :: By ! Barrier loc y-dir
    real(dp) :: Wx ! Barrier width x-dir
    real(dp) :: Wy ! Barrier width y-dir

    integer :: Mx ! number oflattice points in x-dir
    integer :: My 
    integer :: N  ! number of time steps/iterations
    integer :: V_type
    integer :: plot_interval

    logical :: plot_re
  end type
end module
