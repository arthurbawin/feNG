# feNG

feNG is a 1-2D finite element solver for simple incompressible flow applications. It is very much a work in progress.
Its current dependencies are PETSc or MKL PARDISO.

Supported mesh file format is GMSH's native .msh format (versions 2.2 and 4+).

To create a Makefile with CMake and build :

 	mkdir build
 	cd build
 	cmake -DPETSC_DIR=/path/to/petsc -DPETSC_ARCH=architecture ..
  make
