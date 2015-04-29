module simulation
  use constants
  use plotroutines
  use omp_lib
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(psi, x, y, n, Mx, My, Ax, Ay)
    complex(dp), intent(inout) :: psi(:,:)
    complex(dp), intent(in)    :: Ax(:,:,:), Ay(:,:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:)
    integer, intent(in)        :: n, Mx, My

    integer :: i

    do i=1,n
      call solve_nxt(psi, Mx, My, Ax, Ay)
      
      if (mod(i,5)==0) call plot_wavef(psi, x, y, Mx, My)
    enddo
  end subroutine

  subroutine solve_nxt(psi, Mx, My, Ax, Ay)
    complex(dp), intent(inout) :: psi(:,:)
    complex(dp), intent(in)    :: Ax(:,:,:), Ay(:,:,:)
    integer, intent(in)        :: Mx, My

    complex(dp), allocatable :: gx(:), gy(:), Ax_d(:), Ax_l(:), &
                                Ax_u(:), Ay_d(:), Ay_l(:), &
                                Ay_u(:)
    integer                  :: i, info

    allocate(Ax_d(Mx), Ax_l(Mx-1), Ax_u(Mx-1), gx(Mx), Ay_d(My), &
      Ay_l(My-1), Ay_u(My-1), gy(My))

    ! horizontal sweep
    !$omp parallel do private(Ax_d,Ax_u,Ax_l,gx)
    do i=1,My
      ! init temp arrays
      gx = Ax(2,:,i)
      Ax_d = Ax(2,:,i)
      Ax_u = Ax(1,1:Mx-1,i)
      Ax_l = Ax(1,1:Mx-1,i)

      ! explicit part of calculation, mat-vec multiplication
      call zgbmv('N', Mx, Mx, 1, 1, one, conjg(Ax(:,:,i)), 3, psi(:,i), 1, &
        zero, gx, 1)

      ! solve for psi at t=n+1/2
      call zgtsv(Mx, 1, Ax_l, Ax_d, Ax_u, gx, Mx, info)
      psi(:,i) = gx
    enddo
    !$omp end parallel do

    ! vertical sweep
    !$omp parallel do private(Ay_d,Ay_u,Ay_l,gy)
    do i=1,Mx
      ! init temp arrays
      gy = Ay(2,:,i)
      Ay_d = Ay(2,:,i)
      Ay_u = Ay(1,1:My-1,i)
      Ay_l = Ay(1,1:My-1,i)

      ! explicit part of calculation, mat-vec multiplication
      call zgbmv('N', My, My, 1, 1, one, conjg(Ay(:,:,i)), 3, psi(i,:), 1, &
        zero, gy, 1)

      ! solve for psi at t=n+1
      call zgtsv(My, 1, Ay_l, Ay_d, Ay_u, gy, My, info)
      psi(i,:) = gy
    enddo
    !$omp end parallel do

    deallocate(Ax_d, Ax_l, Ax_u, Ay_d, Ay_l, Ay_u, gx, gy)
  end subroutine
end module
