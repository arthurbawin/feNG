#include <iostream>
#include <fstream>

#include "feNG.h"
#include "feFunction.h"
#include "feMesh.h"
#include "feQuadrature.h"
#include "feSpace.h"
#include "feSpaceTriangle.h"
#include "feCncGeo.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feSysElm.h"
#include "feBilinearForm.h"
#include "feSolver.h"
#include "feLinearSystemPETSc.h"
#include "feExporter.h"
#ifdef HAVE_MKL
#include "feLinearSystemMklPardiso.h"
#endif

double fSol(const double t, const std::vector<double> x, const std::vector<double> par) {
  return pow(x[0], 6);
  // return pow(x[0],3);
  // return pow(x[0],2);
}

double fSol2(const double t, const std::vector<double> x, const std::vector<double> par) {
  return pow(x[1], 6);
  // return pow(x[0],3);
  // return pow(x[0],2);
}

double fSource(const double t, const std::vector<double> x, const std::vector<double> par) {
  double kd = par[0];
  return kd * 30. * pow(x[0], 4);
  // return kd*6*x[0];
  // return kd*2;
}

double fSource2(const double t, const std::vector<double> x, const std::vector<double> par) {
  double kd = par[0];
  return kd * 30. * pow(x[1], 4);
  // return kd*6*x[0];
  // return kd*2;
}

int main(int argc, char **argv) {
#ifdef HAVE_PETSC
  petscInitialize(argc, argv);
#endif

  double kd = 1.0;
  feFunction *funSol = new feFunction(fSol, {kd});
  feFunction *funSol2 = new feFunction(fSol2, {kd});
  feFunction *funSource = new feFunction(fSource, {kd});
  feFunction *funSource2 = new feFunction(fSource2, {kd});

  int nIter = 1;
  std::vector<double> normL2(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);

  for(int iter = 0; iter < nIter; ++iter) {
    std::string meshName = "../../data/Square/square_1edge_1.msh";

    feMesh2DP1 *mesh = new feMesh2DP1(meshName, false);
    nElm[iter] = mesh->getNbInteriorElems();

    feSpace1DP3 U_bord = feSpace1DP3(mesh, "U", "Bord", funSol);
    feSpaceTriP3 U_surface = feSpaceTriP3(mesh, "U", "Surface", funSol);

    std::vector<feSpace *> fespace = {&U_bord, &U_surface};
    std::vector<feSpace *> feEssBC = {&U_bord};

    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);

    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);

    int nQuad = 8;
    std::vector<feSpace *> spaceDiffusion2D_U = {&U_surface};
    std::vector<feSpace *> spaceSource2D_U = {&U_surface};
    feBilinearForm *diffU =
      new feBilinearForm(spaceDiffusion2D_U, mesh, nQuad, new feSysElm_2D_Diffusion(kd, nullptr));
    feBilinearForm *sourceU =
      new feBilinearForm(spaceSource2D_U, mesh, nQuad, new feSysElm_2D_Source(1.0, funSource));

    std::vector<feBilinearForm *> formMatrices = {diffU};
    std::vector<feBilinearForm *> formResiduals = {diffU, sourceU};

    feNorm *norm = new feNorm(&U_surface, mesh, nQuad, funSol);
    std::vector<feNorm *> norms = {norm};

    feTolerances tol{1e-9, 1e-8, 3};

    // long int cnt = 0;
    // #pragma omp parallel for private(cnt)
    // for(int i = 0; i < 100; ++i){
    //   // printf("Printing %3d from thread %d\n", i, omp_get_thread_num());
    //   printf("Thread %d has printed %2d times\n", omp_get_thread_num(), cnt++);
    // }

    // printf("Measuring time from here\n");
    // tic();
#ifdef HAVE_MKL    
    feLinearSystemMklPardiso *linearSystem;
    linearSystem = new feLinearSystemMklPardiso(formMatrices, formResiduals, metaNumber, mesh);
    // toc();

    // return 1;

    StationarySolver *solver =
      new StationarySolver(tol, metaNumber, linearSystem, sol, norms, mesh);
    solver->makeSteps(0, fespace);

    std::string vtkFile = "sol.vtk";
    feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);

    delete norm;
    delete linearSystem;
    delete sourceU;
    delete diffU;
    delete sol;
    delete metaNumber;
    delete mesh;
#endif    
  }
  delete funSource;
  delete funSol;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2[2 * i + 1] =
      -log(normL2[2 * i] / normL2[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  }
  printf("%12s \t %12s \t %12s\n", "nElm", "||E||", "p");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e\n", nElm[i], normL2[2 * i], normL2[2 * i + 1]);

#ifdef HAVE_PETSC
  petscFinalize();
#endif

  return 0;
}
