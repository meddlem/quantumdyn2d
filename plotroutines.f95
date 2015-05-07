module plotroutines
  use constants
  use structures 
  implicit none
  private
  public :: plot_wavef, close_plot, animate_plot

contains
  subroutine animate_plot(Q)
    type(modl_par), intent(in) :: Q

    integer :: ret
    
    ! creates fifo pipe: plotfifo.dat
    call system("rm -f plotfifo.dat; mkfifo plotfifo.dat", ret)
    
    ! create a gnuplot command file
    open(10,access = 'sequential',file = 'matplot.plt')
      !write(10,*) 'set term x11' 
      write(10,*) 'set view map'
      write(10,*) 'set palette defined ( 0 "#000090", 1 "#000fff",\' 
      write(10,*) '2 "#0090ff", 3 "#0fffee", 4 "#90ff70", 5 "#ffee00",\' 
      write(10,*) '6 "#ff7000", 7 "#ee0000", 8 "#7f0000")'
      write(10,*) 'set pm3d'
      write(10,*) 'set size ratio',Q%Ly/Q%Lx
      write(10,*) 'set hidden3d'
      write(10,*) 'set xrange [',Q%Lx/2,':',Q%Lx,']'
      write(10,*) 'set yrange [0:',Q%Ly,']'
      write(10,*) 'set xlabel "x"'
      write(10,*) 'set ylabel "y" norotate'
      if (.not. Q%plot_re) then
        !write(10,*) 'set cbrange [0:0.3]'
        write(10,*) 'set cblabel "density"'
      elseif (Q%plot_re) then
        write(10,*) 'set cbrange [-0.2:0.2]'
        write(10,*) 'set cblabel "Re(Psi)"'
      endif
      write(10,*) 'load "loop.plt"'
    close(10)
    
    ! create plot/animate instruction
    open(10,access = 'sequential', file = 'loop.plt')
      write(10,*) 'splot "< cat plotfifo.dat" using 1:2:3 with pm3d notitle'
      ! write(10,*) '"" using 1:2:4 with lines' !with lines
      write(10,*) 'pause 0.2'
      write(10,*) 'reread'
    close(10)
    
    ! now fork instance of gnuplot to plot/animate the lattice
    call system("gnuplot matplot.plt &",ret)
  end subroutine
  
  subroutine plot_wavef(psi, x, y, Q)
    complex(dp), intent(in)    :: psi(:,:)
    real(dp), intent(in)       :: x(:,:), y(:,:)
    type(modl_par), intent(in) :: Q

    integer :: i, j
    character(50) :: rfmt

    rfmt = '(F10.5,1X,F10.5,1X,F10.5)' 
    
    if (Q%plot_re) then
      open(11,access = 'sequential',status = 'replace',file = 'plotfifo.dat')
        do i = 1,Q%Mx
          do j = 1,Q%My
            write(11,rfmt) x(i,j), y(i,j), real(psi(i,j)) ! write plot data
          enddo
          write(11,*) '' ! add space between different xvals, for gnuplot
        enddo
      close(11)
    else 
      open(11,access = 'sequential',status = 'replace',file = 'plotfifo.dat')
        do i = 1,Q%Mx
          do j = 1,Q%My
            write(11,rfmt) x(i,j), y(i,j), abs(psi(i,j))**2 ! write plot data
          enddo
          write(11,*) '' ! add space between different xvals, for gnuplot
        enddo
      close(11)
    endif
  end subroutine

  subroutine close_plot()
    call system('pkill gnuplot')
    call system('rm -f plotfifo.dat')
  end subroutine
end module 
