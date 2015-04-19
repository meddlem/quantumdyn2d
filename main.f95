program main
  use constants
  use initialize
  use io
  implicit none

  complex(dp), allocatable :: psi(:)
  real(dp), allocatable    :: x(:), V(:)
  real(dp) :: k, dx, L
  integer  :: M

  call user_in(k,dx,L,M)
  allocate(psi(M),x(M),V(M))
  call init_wavef(psi,x,k,L,M)
  call init_V(V,x,L)
  call run_sim(psi,x,V,L,M)
  deallocate(psi,x)
end program
