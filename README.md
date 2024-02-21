# Simulate percolating network

A fast incremental network solver, using low rank matrix updates, as part of a deterministic simulation of electrical interactions between tin particles to investigate if they act as a hardware-based neural network.

[Download the report](Low_rank_updates.pdf).

A system of groups (groups of tin particles) and the conductances between those groups is provided via .mtx format or using a triplet object. This data must be in lower triangular form. The program `make_symmetric.cpp` can be used to convert .mtx files into lower triangular form.

As demonstrated in `example.c`, this triplet is provided to the solver and conductance and voltage updates can be applied over a number ofg iterations.

The solver works by first computing the Cholesky decomposition of the system then solving for the voltage at every group. For each iteration, conductances can be updated. The decomposition is then updated without recomputing from the system matrix.

The CHOLMOD library is used to perform the matrix computations. CHOLMOD is a subset of [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html). CHOLMOD requires BLAS, used to perform elementary linear algebra computations optimally on your CPU. If your computer does not already have BLAS, consider installing OpenBLAS.

The solver is written in C++ but can also interface with C so that it may be used within the existing simulator.

## Installation

* Download, compile, and install [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html). Only the CHOLMOD module is used (and its dependances).
* Set the environment variable `export SSPATH=/path/to/SuiteSparse-X.X.X`. You may also need to change the compiler command you are using at the `CC=..` and `CCC=` lines in the MakeFile. I used GNU gcc and g++ 9.
* Run `make` to compile the example.

## Usage
Assuming all the files compiled successfully, run the example code with `make && ./example.o `. You should see the following output after the compilation info.

```
Input:
CHOLMOD triplet:  4-by-4, nz 5, lower.  OK
3 2 0.500000
3 0 0.200000
2 1 0.300000
2 0 0.100000
1 0 0.400000

Iter: 0

Group Index: Voltage
0: 0.622222
1: 1
2: 0
3: 0.177778

Iter: 1

Group Index: Voltage
0: 0.52766
1: 1
2: 0
3: 0.27234
```
This solves the example system given in the accompaning documentation, but switches between adding and subtracting `1.0/3.0` from the connection between group 0 and group 3 every iteration.

## Using .mtx files
You should be able to find .mtx exporters online for whatever language/software you are using. CHOLMOD also has inbuild functions for dealing with .mtx files. Make sure you export in coordinate format, if it gives you the option. You may wish to use `make_symmetric.cpp` to ensure your .mtx file is in the correct format.

Run `./make_symmetric.o < matrices/simple.mtx > simple.out.mtx`.

This will take (from stdin) the `simple.mtx` file from the `matrices` folder and checks if it is positive definate and symmetric. If it is, it writes to (stdout) `simple.out.mtx`.
