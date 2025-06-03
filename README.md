This package is the OpenACC version of ERmod. For details on how to use it and information about the original code, see https://sourceforge.net/projects/ermod/.

ERmod (Energy Representation Module) is a program to calculate the solvation free energy based on the energy representation method. The program allows users to calculate the solvation free energy in arbitrary solvent, including inhomogeneous systems, and runs in cooperation with state-of-art molecular simulation softwares, such as NAMD, GROMACS, and AMBER.

## Requirements
- NVIDIA HPC SDK (https://developer.nvidia.com/hpc-sdk)

## Typical Installation
This package is built with autotools. To compile the program,
1. cd erdst (or cd slvfe)
2. autoreconf -fiv
3. configure the package with "configure" (see below)
4. then compile the package with "make"

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

## Reference
ERmod-OpenACC: GPU Acceleration of Solvation Free Energy Calculation with Energy-Representation Theory
Yutaka Maruyama, Hidekazu Kojima, Nobuyuki Matubayasi
J. Comput. Chem. (in press)
https://doi.org/10.1002/jcc.70152