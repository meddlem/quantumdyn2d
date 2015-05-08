module structures
  use constants
  implicit none
  
  ! this module contains additional datatypes used in simulation
  ! default values are set for all parameters, later changes can be made in
  ! init_param subroutine
  type modl_par 
    ! data type that contains all model parameters 
    ! initial wavevector in x, y dir
    real(dp) :: kx = 0._dp 
    real(dp) :: ky = 0._dp 
    ! lattice length in x,y dir
    real(dp) :: Lx = 9._dp
    real(dp) :: Ly = 9._dp
    ! spatial and temporal discretization steps
    real(dp) :: dx = 0.05_dp
    real(dp) :: dt = 0.02_dp
    ! number of time steps/iterations
    integer  :: N  = 6000
    ! time scale potential variation
    real(dp) :: tau = 100._dp 
    ! Slit width x,y-dir
    real(dp) :: Wx 
    real(dp) :: Wy 
    ! Double slit spacing 
    real(dp) :: D  
    ! number of lattice points in x,y dir
    integer  :: Mx 
    integer  :: My
    ! defines potential, initial wf
    character(3) :: sim_type = 'har' 
  end type

  type plt_par
    ! data type that contains all plot parameters
    
    ! number of iterations between plots
    integer  :: plot_interval = 10 
    ! min and max z-value (color range)
    real(dp) :: rng(2) = [-0.2_dp, 0.2_dp] 
    ! plot real part or density of wavefunction
    logical  :: plot_re = .false.
  end type
end module
