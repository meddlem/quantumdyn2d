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
    call system("rm -f plotfifo.dat; mkfifo plotfifo.dat",ret)
    
    ! create a gnuplot command file
    open(10,access = 'sequential',file = 'matplot.plt')
      !write(10,*) 'set term x11' 
      write(10,*) 'set style line 1 lt 1 lc rgb "blue" lw 2 pt 2 ps 0.6'
      write(10,*) 'set style line 2 lt 1 lc rgb "red" lw 2 pt 2 ps 0.6'
      write(10,*) 'set style line 3 lt 1 lc rgb "purple" lw 2 pt 2 ps 0.6'
      write(10,*) 'set style line 4 lt 1 lc rgb "black" lw 1 pt 2 ps 0.6'
      write(10,*) 'set grid'
      write(10,*) 'set xlabel "x"'
      !write(10,*) 'set ylabel "Psi^2"'
      write(10,*) 'set xrange [0:',Q%L,']'
      write(10,*) 'set yrange [-1.1:1.1]'
      write(10,*) 'load "loop.plt"'
    close(10)
    
    ! create plot/animate instruction
    open(10,access = 'sequential', file = 'loop.plt')
      !write(10,*) 'plot "< cat plotfifo.dat" using 1:2 with lines ls 1 title "Re(Psi)",\'
      ! write(10,*) '"" using 1:3 with lines ls 2 title "Im(Psi)",\'
      write(10,*) 'plot "< cat plotfifo.dat" using 1:4 with lines ls 1 title "P",\'
      write(10,*) '"" using 1:5 with lines ls 4 title "V"'
      write(10,*) 'pause 0.05'
      write(10,*) 'reread'
    close(10)
    
    ! now fork instance of gnuplot to plot/animate the lattice
    call system("gnuplot matplot.plt &",ret)
  end subroutine
  
  subroutine plot_wavef(psi, x, V, Q)
    complex(dp), intent(in)    :: psi(:)
    real(dp), intent(in)       :: x(:), V(:)
    type(modl_par), intent(in) :: Q

    integer :: i
    character(50) :: rfmt

    rfmt = '(F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5)' 
    
    open(11,access = 'sequential',status = 'replace',file = 'plotfifo.dat')
      do i = 1,Q%M
        write(11,rfmt) x(i), real(psi(i)), aimag(psi(i)), abs(psi(i))**2, V(i) ! write plot data to pipe 
      enddo
    close(11)
  end subroutine

  subroutine close_plot()
    call system('pkill gnuplot')
    call system('rm -f plotfifo.dat')
  end subroutine
end module 
