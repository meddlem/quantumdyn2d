module plotroutines
  use constants
  use structures 
  implicit none
  private
  public :: plot_wavef, close_plot, animate_plot, snapshot

contains
  subroutine animate_plot(Q, P)
    ! creates animate plot instruction for gnuplot
    type(modl_par), intent(in) :: Q
    type(plt_par), intent(in)  :: P

    integer :: ret
    
    ! creates fifo pipe: plotfifo.dat
    call system("rm -f plotfifo.dat; mkfifo plotfifo.dat", ret)
    
    ! create a gnuplot command file
    open(10,access = 'sequential',file = 'matplot.plt')
      ! write(10,*) 'set term x11' ! uncomment line in case of errors
      write(10,*) 'set palette defined ( 0 "#000090", 1 "#000fff",\' 
      write(10,*) '2 "#0090ff", 3 "#0fffee", 4 "#90ff70", 5 "#ffee00",\' 
      write(10,*) '6 "#ff7000", 7 "#ee0000", 8 "#7f0000")'
      write(10,*) 'set size ratio', Q%Ly/Q%Lx
      write(10,*) 'dx =', Q%dx
      write(10,*) 'set xrange [0:',Q%Lx,']'
      write(10,*) 'set yrange [0:',Q%Ly,']'
      write(10,*) 'set xlabel "x"'
      write(10,*) 'set ylabel "y" norotate'
      if (.not. P%plot_re) then
        write(10,*) 'set cbrange [0:',P%rng(2),']'
        write(10,*) 'set cblabel "density"'
      elseif (P%plot_re) then
        write(10,*) 'set cbrange [',P%rng(1),':',P%rng(2),']'
        write(10,*) 'set cblabel "Re(Psi)"'
      endif
      write(10,*) 'load "loop.plt"'
    close(10)
    
    ! create plot/animate instruction
    open(10,access = 'sequential', file = 'loop.plt')
      write(10,*) 'plot "< cat plotfifo.dat" u ($1*dx):($2*dx):3 matrix \'
      write(10,*) 'with image notitle'
      write(10,*) 'pause 0.2'
      write(10,*) 'reread'
    close(10)
    
    ! now fork instance of gnuplot to plot/animate the lattice
    call system("gnuplot matplot.plt &",ret)
  end subroutine
  
  subroutine plot_wavef(psi, Q, P, auto)
    complex(dp), intent(in)    :: psi(:,:)
    type(modl_par), intent(in) :: Q
    type(plt_par), intent(in)  :: P
    logical, intent(in)        :: auto

    integer :: i, j
    character(50) :: rfmt, filename

    write(rfmt, '(A,I0,A)') '(',Q%Mx,'(10X,F10.5))' 
    
    if (auto) then
      filename = 'plotfifo.dat'
    else
      filename = 'wf.dat'
    endif

    ! write plotdata to fifo pipe
    if (.not. P%plot_re) then
      open(11,status = 'replace',file = filename)
        do i = 1,Q%My
          write(11,rfmt) (abs(psi(j,i))**2, j=1,Q%Mx)
        enddo
      close(11)
    else
      open(11,status = 'replace',file = 'plotfifo.dat')
        do i = 1,Q%My
          write(11,rfmt) (real(psi(j,i)), j=1,Q%Mx)
        enddo
      close(11)
    endif
  end subroutine

  subroutine close_plot()
    call system('pkill gnuplot')
    call system('rm -f plotfifo.dat')
  end subroutine

  subroutine snapshot(psi, Q, P, pl_no)
    ! creates single plot off wavefunction
    complex(dp), intent(in)    :: psi(:,:)
    type(modl_par), intent(in) :: Q
    type(plt_par), intent(in)  :: P
    integer, intent(in)        :: pl_no

    integer       :: ret
    character(40) :: filename

    write(filename,'(A,I0,A)') 'set output "plot',pl_no,'.png"'
    
    ! create a gnuplot command file
    open(10,access = 'sequential',file = 'wf.plt')
      write(10,*) 'set term pngcairo size 640,480'
      !write(10,*) 'set term epscairo size 12cm,9cm enhanced font "Verdana,15"'
      write(10,*) filename
      write(10,*) 'set palette defined ( 0 "#000090", 1 "#000fff",\' 
      write(10,*) '2 "#0090ff", 3 "#0fffee", 4 "#90ff70", 5 "#ffee00",\' 
      write(10,*) '6 "#ff7000", 7 "#ee0000", 8 "#7f0000")'
      write(10,*) 'set size ratio', Q%Ly/Q%Lx
      write(10,*) 'dx =', Q%dx
      write(10,*) 'set xrange [0:',Q%Lx,']'
      write(10,*) 'set yrange [0:',Q%Ly,']'
      write(10,*) 'set xlabel "x"'
      write(10,*) 'set ylabel "y" norotate'
      if (.not. P%plot_re) then
        write(10,*) 'set cbrange [0:',P%rng(2),']'
        write(10,*) 'set cblabel "density"'
      elseif (P%plot_re) then
        write(10,*) 'set cbrange [',P%rng(1),':',P%rng(2),']'
        write(10,*) 'set cblabel "Re(Psi)"'
      endif
      write(10,*) 'plot "wf.dat" u ($1*dx):($2*dx):3 matrix \'
      write(10,*) 'with image notitle'
    close(10)

    call plot_wavef(psi, Q, P, .false.)
    
    ! now fork instance of gnuplot to plot/animate the lattice
    call system("gnuplot wf.plt &",ret)
  end subroutine
end module 
