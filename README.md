This is an EXPERIMENTAL version. For information about the original code, see https://sourceforge.net/projects/ermod/.

ERmod (Energy Representation Module) is a program to calculate the solvation free energy based on the energy representation method. The program allows users to calculate the solvation free energy in arbitrary solvent, including inhomogeneous systems, and runs in cooperation with state-of-art molecular simulation softwares, such as NAMD, GROMACS, and AMBER.

## Requirements
- NVIDIA HPC SDK (https://developer.nvidia.com/hpc-sdk)
- BLAS & LAPACK library (BLAS and LAPACK, OpenBLAS, MKL etc) for erdst

## Typical Installation
This package is built with autotools. To compile the program,
1. cd erdst (or cd slvfe)
2. aclocal
3. autoconf
4. automake (or automake --add-missing)
5. configure the package with "configure" (see below)
6. then compile the package with "make"

## Configuration
### erdst
For multi-GPU,
```
$ ./configure FC=nvfortran CC=nvc MPIFC=mpif90
```
or for single-GPU,
```
$ ./configure --disable-mpi FC=nvfortran CC=nvc
```

If your computer has OpenBLAS library, try:
```
$ ./configure --with-openblas [other options]
```
If you have Intel MKL, configure program with:
```
$ ./configure --with-mkl [other options]
```

It is set to single-precision calculations by default.
If you want to use double-precision calculations,
```
$ ./configure --enable-double [other options]
```

If configuration finishes successfully, type
```
$ make
```
to start compilation.

### slvfe
slvfe command only supports single-GPU,
```
$ ./configure FC=nvfortran
```
If configuration finishes successfully, type
```
$ make
```
to start compilation.

Current configuration script only supports nvfortran as a compiler.
