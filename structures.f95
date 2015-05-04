module structures
  use constants
  implicit none
  
  type model_parameters 
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
end module
