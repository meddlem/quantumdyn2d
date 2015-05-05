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
    real(dp) :: a  ! time constant potential

    integer :: Mx ! lattice size in x-dir
    integer :: My
    integer :: N  ! number of time steps
    integer :: V_type

    logical :: plot_re
  end type
  
  type Ops
    ! data type that contains the ADI matrix operators

    complex(dp), allocatable :: Ax(:,:,:) 
    complex(dp), allocatable :: Ay(:,:,:) 
  end type
end module
