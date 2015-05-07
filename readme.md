# 2d quantum dynamics simulation

Time integration of the 2d Schrodinger equation, using a second order Alternating Direction Implicit (ADI) scheme. 
The resulting systems are solved using LAPACK routines for tridiagonal matrices. 
The result of the simulation is an animation time evolution of the density (default), or the real part of the wavefunction (use flag: -r).
##3 experiments may be selected by the user: 

1. Default: Adiabatic potential variation: Harmonic well to ISQW
2. Use -s flag for single slit diffraction
3. Use -d flag for double slit diffraction

Example: $./main -d 
