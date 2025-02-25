# feNG

[![Debug](https://github.com/arthurbawin/feNG/blob/master/.github/workflows/debug.yml/badge.svg)](https://github.com/arthurbawin/feNG/blob/master/.github/workflows/debug.yml)

feNG (short for finite element NextGeneration, a somewhat ironic name) is a C++ finite element library for the resolution of heat transfer and incompressible flow applications on simplicial meshes and in a continuous Galerkin setting. The library was loosely inspired by [MFEM](https://mfem.org/) and takes the form of a toolbox, where continuous operators are assembled to form a partial differential equation.

The library features:
- Lagrange finite elements of order 0 to 4 on triangles, and of arbitrary order on tets. Finite elements for the resolution of saddle-point problems (e.g. the Stokes equations) are also available, such as the mixed Crouzeix-Raviart $(\mathcal{P}^{2+},\mathcal{P}^{-1})$ element for velocity and pressure discretization (quadratic Lagrange enriched with cubic bubble combined with discontinuous linear Lagrange).
- Support for high-order (quadratic, a.k.a. $\mathcal{P}^2$) meshes
- Error estimation based on the [Zhang & Naga recovery operator](https://epubs.siam.org/doi/abs/10.1137/S1064827503402837)
- Riemannian metric computation for anisotropic adaptive remeshing with the [Mmg](https://www.mmgtools.org/) library

It is very much a work in progress.

## Dependencies
- feNG interfaces both [PETSc](https://petsc.org/) and Intel MKL Pardiso as sparse linear solvers, at least one of these libraries is required.
- The supported mesh file format is [Gmsh](https://gmsh.info/)'s native .msh format (version 4+).
- Visualization files are written in legacy VTK format for visualization with e.g. [ParaView](https://www.paraview.org/)
- Both [Mmg](https://www.mmgtools.org/) and Gmsh's development kit (SDK) are required for metric computations and anisotropic mesh adaptation with linear finite elements. For higher-order finite elements, the optimization library (simplex method) [SoPlex](https://soplex.zib.de/) is additionally required.

## Installation

The easiest approach is to use CMake and specify the location of either PETSc or Intel MKL (or both).
For instance, with PETSc and starting from the feNG/ directory :
```
mkdir build && cd build
cmake -DENABLE_PETSC=1 -DPETSC_DIR=/path/to/petsc -DPETSC_ARCH=architecture ..
make
```
This will compile the library, the executables and the regression tests.
To check everything works as intended, simply run :
```
ctest
```

## Authors
- Arthur Bawin (main developer)
- Baptiste Berlioux
- Simon Berthelin
- André Garon
