module plotroutines
  use constants
  implicit none
  private
  public :: line_plot, plot_wavef, close_plot, animate_plot

contains
  subroutine animate_plot(Lx,Ly)
    real(dp), intent(in) :: Lx, Ly
    integer :: ret
    
    ! creates fifo pipe: plotfifo.dat
    call system("rm -f plotfifo.dat; mkfifo plotfifo.dat",ret)
    
    ! create a gnuplot command file
    open(10,access = 'sequential',file = 'matplot.plt')
      write(10,*) 'set term x11' 
      !write(10,*) 'set style line 1 lt 1 lc rgb "blue" lw 2 pt 2 ps 0.6'
      write(10,*) 'set grid back'
      !write(10,*) 'set palette grey'
      write(10,*) 'set palette defined ( 0 "#000090", 1 "#000fff",\' 
      write(10,*) '2 "#0090ff", 3 "#0fffee", 4 "#90ff70", 5 "#ffee00",\' 
      write(10,*) '6 "#ff7000", 7 "#ee0000", 8 "#7f0000")'
      write(10,*) 'set view map'
      write(10,*) 'set pm3d'
      write(10,*) 'set hidden3d'
      write(10,*) 'set xlabel "x"'
      write(10,*) 'set ylabel "y"'
      write(10,*) 'set xrange [0:',Lx,']'
      write(10,*) 'set yrange [0:',Ly,']'
      !write(10,*) 'set zrange [0:0.2]'
      write(10,*) 'load "loop.plt"'
    close(10)
    
    ! create plot/animate instruction
    open(10,access = 'sequential', file = 'loop.plt')
      write(10,*) 'splot "< cat plotfifo.dat" using 1:2:3 with pm3d'! lines linecolor rgb "blue"' !with pm3d
!      write(10,*) '"" using 1:2:4 with lines' !with lines
      write(10,*) 'pause 0.2'
      write(10,*) 'reread'
    close(10)
    
    ! now fork instance of gnuplot to plot/animate the lattice
    call system("gnuplot matplot.plt &",ret)
  end subroutine
  
  subroutine plot_wavef(psi, x, y, Mx, My)
    complex(dp), intent(in) :: psi(:,:)
    real(dp), intent(in)    :: x(:,:), y(:,:)
    integer, intent(in)     :: Mx, My

    integer :: i, j
    character(50) :: rfmt

    rfmt = '(F10.5,1X,F10.5,1X,F10.5)' 
    
    open(11,access = 'sequential',status = 'replace',file = 'plotfifo.dat')
      do i = 1,Mx
        do j = 1,My
          write(11,rfmt) x(i,j), y(i,j), abs(psi(i,j))**2 ! write plot data to pipe 
        enddo
        write(11,*) '' ! add space between different xvals, needed for gnuplot
      enddo
    close(11)
  end subroutine

  subroutine close_plot()
    call system('pkill gnuplot')
    call system('rm -f plotfifo.dat')
  end subroutine

  subroutine line_plot(x,y1,xlabel,ylabel,label1,title,plot_no,y2,label2)
    real(dp), intent(in) :: x(:), y1(:)
    real(dp), intent(in), optional :: y2(:)
    character(*), intent(in) :: xlabel, ylabel, label1, title
    character(*), intent(in), optional :: label2
    integer, intent(in) :: plot_no
    character(1024) :: filename
    integer :: i, ret, m
    real(dp) :: xmin, xmax, xrange(2), ymin, ymax, yrange(2)
    
    if (size(y1)/=size(x)) print *, "error, arguments must be same size"
    
    if (present(y2)) then
      if (size(y1)/=size(y2)) then 
        print *, "error, arguments must be same size"
        return
      endif
    endif
    
    m = size(x)
    xmin = minval(x)
    xmax = maxval(x)
    ymin = minval(y1)
    ymax = maxval(y1)

    xrange = [xmin-(xmax-xmin)*0.1_dp, xmax+(xmax-xmin)*0.1_dp]
    yrange = [ymin-(ymax-ymin)*0.1_dp, ymax+(ymax-ymin)*0.1_dp]

    open(10,access = 'sequential',file = 'xydata.dat')
      do i=1,m
        if (present(y2)) then
          write(10,*) x(i),y1(i),y2(i) ! write datapoints to file
        else
          write(10,*) x(i),y1(i) ! write datapoints to file
        endif
      enddo
    close(10)
    
    ! create gnuplot command file
    write(filename,'(A,I0,A)') 'set output "plot',plot_no,'.png"'
    open(10,access = 'sequential',file = 'gplot.txt')
      
      ! set output terminal  
      write(10,*) 'set term pngcairo size 640,480 enhanced font "Verdana,10"'
      !write(10,*) 'set term epscairo size 15cm,9cm font "Verdana,15"'

      write(10,*) filename
      ! set line color definitions
      write(10,*) &
        'set style line 1 lt 1 lc rgb "#ff0000" lw 2 #red'
      write(10,*) &
        'set style line 2 lt 1 lc rgb "#0000ff" lw 2 #blue'
      ! axes 
      write(10,*) 'set style line 11 lc rgb "#808080" lt 1'
      write(10,*) 'set border 31 back ls 11'
      write(10,*) 'set tics nomirror scale 0.75'
      write(10,*) 'set key right center'
      write(10,*) 'set mxtics 2'
      write(10,*) 'set mytics 2'
      ! grid 
      write(10,*) 'set style line 12 lc rgb "#808080" lt 0 lw 1'
      write(10,*) 'set grid back ls 12'
      ! plotrange
      write(10,*) 'set xrange [',xrange(1),':',xrange(2),']'
      write(10,*) 'set yrange [',yrange(1),':',yrange(2),']'
      ! plot labels
      write(10,*) 'set title "'//TRIM(title)//'"'
      write(10,*) 'set xlabel '//'"'//TRIM(xlabel)//'"'
      write(10,*) 'set ylabel '//'"'//TRIM(ylabel)//'"'
      
      if (m > 0) then
        write(10,*) 'plot "xydata.dat" using 1:2 with line ls 10 t "", \'
        write(10,*) &
        ' "xydata.dat" using 1:2 with line ls 1 t "'//TRIM(label1)//'", \'
        if (present(y2)) then
          write(10,*) &
          ' "xydata.dat" using 1:3 with line ls 2 t "'//TRIM(label2)//'"'
        endif
      endif
    close(10)

    ! now call gnuplot and plot the curves
    call system('gnuplot gplot.txt',ret)
    call system('rm gplot.txt; rm xydata.dat',ret)
  end subroutine 
end module 
