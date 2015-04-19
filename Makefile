FC = gfortran
FFLAGS = -ffast-math -Wall -march=native -O3 -fopenmp #-Warray-temporaries #compiler flags
LDFLAGS = -fopenmp #link flags

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

PROG = main #program name

#required objects: 
OBJS =
OBJS += main.o

all: $(PROG)

main: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f95
	$(COMPILE) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) $(PROG) $(OBJS) *.mod
	$(RM) *.png *.txt *.plt *.dat
