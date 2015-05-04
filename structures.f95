module structures
  use constants
  implicit none
  
  type modl_par 
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
   complex(dp), allocatable :: Ax(:,:,:) 
   complex(dp), allocatable :: Ay(:,:,:) 
  end type
end module
