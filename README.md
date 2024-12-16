ERmod (Energy Representation Module) is a program to calculate the solvation free energy based on the energy representation method. The program allows users to calculate the solvation free energy in arbitrary solvent, including inhomogeneous systems, and runs in cooperation with state-of-art molecular simulation softwares, such as NAMD, GROMACS, and AMBER.

## Typical Installation
This package is built with autotools. To compile the program,
1. configure the package with "configure",
2. then compile the package with "make".

If you have Intel MKL, configure program with:
    $ ./configure --with-mkl

If your computer has BLAS & LAPACK library, and FFTW (http://www.fftw.org/) version 3, try:
    $ ./configure --with-fft=fftw

If configuration finishes successfully, type
    $ make
to start compilation.

## Non-typical Installation
If you want to use a specific version of Intel MKL, try
    ./configure --with-mkl=(version number) --with-fft=mkl

If you want to use a specific BLAS library, try
    ./configure --with-blas=(BLAS library specification)

Current configuration script only supports gfortran and ifort as a compiler.
For other compilers, please specify F77 / FC / FFLAGS / FCFLAGS environment variables to appropriate ones.
