module structures
  use constants
  implicit none
  
  type modl_par 
    ! data type that contains all model parameters
    
    real(dp) :: kx 
    real(dp) :: ky
    real(dp) :: Lx
    real(dp) :: Ly
    real(dp) :: dx
    real(dp) :: dt

    integer :: Mx
    integer :: My
    integer :: N
    integer :: V_type

    logical :: plot_re
  end type
  
  type Ops
    ! data type that contains the ADI matrix operators

    complex(dp), allocatable :: Ax(:,:,:) 
    complex(dp), allocatable :: Ay(:,:,:) 
  end type
end module
