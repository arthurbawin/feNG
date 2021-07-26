# feNG

feNG is a finite element solver for simple incompressible flow applications. It is very much a work in progress.
Its current dependencies are PETSc or PARDISO. Current supported finite element spaces are Lagrange elements of degree 1 to 4 on triangles, 
and implemented weak forms allow to solve 
<<<<<<< HEAD
    * Poisson equation
    * Stokes equations
    * Navier-Stokes equation
=======
- Poisson equation
- Stokes equations
- Navier-Stokes equations
>>>>>>> 2251dfe91238f3e62b7d9129bdb299b0bbdfe57e

Supported mesh file format is GMSH's native .msh format (versions 2.2 and 4+).

To create a Makefile with CMake and build :

<<<<<<< HEAD
 	mkdir build
 	cd build
 	cmake -DPETSC_DIR=/path/to/petsc -DPETSC_ARCH=architecture ..
    make
=======
    mkdir build
    cd build
    cmake -DPETSC_DIR=/path/to/petsc -DPETSC_ARCH=architecture ..
    make
>>>>>>> 2251dfe91238f3e62b7d9129bdb299b0bbdfe57e
