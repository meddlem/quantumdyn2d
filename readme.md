# 2d quantum dynamics simulation

Time integration of the 2d Schrodinger equation, using a second order Alternating Direction Implicit (ADI) scheme. 
The resulting systems are solved using LAPACK and BLAS routines for tridiagonal matrices. The program relies on gnuplot for visualization of the wavefunction. 
The result of the simulation is an animation time evolution of the density (default), or the real part of the wavefunction (use flag: -r).
##3 experiments may be selected by the user: 

1. Default: excited state of harmonic potential
2. Use --dsl flag for double slit diffraction, with gaussian wavepacket
3. Use --hsq flag for adiabatic potential change: Harmonic well to ISQW
4. Use --hqa flag for adiabatic potential change: Harmonic well to quartic potential

Example: $./main --dsl 
