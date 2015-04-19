module initialize
  use constants
  implicit none
  private
  public :: init_wavef, init_V, init_ops

contains
  
  subroutine init_wavef(psi_0,x,dx,k,L,M)
    complex(dp), intent(inout) :: psi_0(:) 
    real(dp), intent(inout)    :: x(:)
    real(dp), intent(in)       :: dx, k, L
    integer, intent(in)        :: M

    integer :: i
    
    do i = 1,M
      x(i) = i*dx
    enddo

    psi_0 = 1/sqrt(L)*exp(cmplx(0._dp,k*x,dp))
    psi_0(1) = 0._dp
    psi_0(M) = 0._dp
  end subroutine

  subroutine init_V(V,x,L)
    real(dp), intent(in) :: x(:), L
    real(dp), intent(inout) :: V(:)
    
    V = 0._dp
    ! block potential
!    where(x>3*L/4) V = 1._dp
  end subroutine
    
  subroutine init_ops(opp_d,opp_du,opm,V,dt,dx,M)
    complex(dp), intent(inout) :: opp_d(:), opp_du(:), opm(:,:)
    real(dp), intent(in)       :: V(:), dt, dx
    integer, intent(in)        :: M

    integer :: i, j
    character(30) :: rowfmt

    opm = cmplx(0._dp,0._dp,dp)
    ! construct operators 
    do i = 1,M
      opm(i,i) = cmplx(1._dp , (-dt*2._dp/(dx**2) - V(i))/2._dp , dp)
      opp_d(i) = cmplx(1._dp , (dt*2._dp/(dx**2) + V(i))/2._dp , dp)

      if (i>1) then
        opm(i,i-1) = cmplx(0._dp , dt*0.5_dp/(dx**2) , dp)
      endif

      if (i<M) then
        opm(i,i+1) = cmplx(0._dp , dt*0.5_dp/(dx**2) , dp)
        opp_du(i) = cmplx(0._dp , -dt*0.5_dp/(dx**2) , dp)
      endif
    enddo

    write(rowfmt,'(A,I4,A)') '(',M,'(2X,F4.1,1X,F4.1))'

    do i=1,M
      write(*,rowfmt) (real(opm(i,j)),aimag(opm(i,j)),j=1,M)
    enddo

  end subroutine
end module 
