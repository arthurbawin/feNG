#include <iostream>
#include <functional>
#include <vector>
#include <cmath>

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
#include "feLinearSystemMklPardiso.h"
#include "feExporter.h"
#define USING_PETSC
#ifdef HAVE_PETSC
#include "feLinearSystemPETSc.h"
#endif
#ifdef HAVE_MKL
#include "feLinearSystemMklPardiso.h"
#endif

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   return sin(2 * M_PI * t) * cos(2 * M_PI * x1) + 1;
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double c1 = par[0];
//   return -2 * M_PI * cos(2 * M_PI * t) * cos(2 * M_PI * x1) -
//          c1 * 2 * M_PI * 2 * M_PI * sin(2 * M_PI * t) * cos(2 * M_PI * x1);
// }

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double c2 = par[1];
//   double x1 = x[0];
//   return x1 * (x1 * x1 - 25.) * pow(t, 6) + x1 * x1 * x1 * pow(t, c2) +1.;
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double c1 = par[0];
//   double c2 = par[1];
//   double beta = (c2 == 0.) ? 0. : c2 * pow(t, c2 - 1.);
//   double x1 = x[0];
//   return -x1 * (x1 * x1 - 25.) * 6. * pow(t, 5.) + c1 * 6. * x1 * pow(t, 6.) +
//          c1 * 6. * x1 * pow(t, c2) - x1 * x1 * x1 * beta;
// }

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  return x1 * x1 * x1 * pow(t, 6.);
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  return -6. * pow(t, 5.) * x1 * x1 * x1 + 6 * x1 * c1 * pow(t, 6.);
}

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double a = par[1];
//   return (x1*x1*+1)*pow(t, a);
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double c1 = par[0];
//   double a = par[1];
//   double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
//   return -beta*(x1*x1*+1) + c1*2*pow(t, a);
// }

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double a = par[1];
//   return pow(t, a);
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double c1 = par[0];
//   double a = par[1];
//   double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
//   return -beta;
// }

int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double xa = 0.;
  double xb = 5.;
  double exposant = 6.;
  double kd = 0.1;
  std::vector<double> par = {kd, exposant};
  feFunction *funSol = new feFunction(fSol, par);
  feFunction *funSource = new feFunction(fSource, par);

  int nIter = 3;
  std::vector<double> normL2_BDF2(2 * nIter, 0.0);
  std::vector<double> normL2_DC3(2 * nIter, 0.0);

  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);

  for(int iter = 0; iter < nIter; ++iter) {
    std::string meshName = "../../data/Square/square_2.0_" + std::to_string(iter + 1) + ".msh";
    // Maillage
    feMesh2DP1 *mesh = new feMesh2DP1(meshName, false);
    nElm[iter] = mesh->getNbInteriorElems();
    // Espaces d'interpolation

    feSpace1DP3 U_angle = feSpace1DP3(mesh, "U", "Angle", funSol);
    feSpace1DP3 U_haut = feSpace1DP3(mesh, "U", "Haut", funSol);
    feSpaceTriP3 U_surface = feSpaceTriP3(mesh, "U", "Surface", funSol);
    feSpace1DP3 U_gauche = feSpace1DP3(mesh, "U", "Gauche", funSol);
    std::vector<feSpace *> fespace = {&U_angle, &U_surface, &U_haut, &U_gauche};
    std::vector<feSpace *> feEssBC = {&U_angle, &U_haut, &U_gauche};
    // Numerotations
    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // Solution
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 5 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feSolution *solDC3 = new feSolution(mesh, fespace, feEssBC, metaNumber);
    // solDC3->initializeTemporalSolution(t0, t1, nTimeSteps);
    // solDC3->initializeUnknowns(mesh, metaNumber);
    // solDC3->initializeEssentialBC(mesh, metaNumber);
    // Formes (bi)lineaires
    int degQuad = 11;
    std::vector<feSpace *> spaceDiffusion2D_U = {&U_surface};
    std::vector<feSpace *> spaceSource2D_U = {&U_surface};
    std::vector<feSpace *> spaceMasse2D_U = {&U_surface};

    feBilinearForm *diff_U_M2D =
      new feBilinearForm(spaceDiffusion2D_U, mesh, degQuad, new feSysElm_2D_Diffusion(kd, nullptr));
    feBilinearForm *source_U_M2D =
      new feBilinearForm(spaceSource2D_U, mesh, degQuad, new feSysElm_2D_Source(1.0, funSource));
    feBilinearForm *masse_U_M2D =
      new feBilinearForm(spaceMasse2D_U, mesh, degQuad, new feSysElm_2D_Masse(1.0, nullptr));

    std::vector<feBilinearForm *> formMatrices = {diff_U_M2D, masse_U_M2D};
    std::vector<feBilinearForm *> formResiduals = {diff_U_M2D, masse_U_M2D, source_U_M2D};
    // Norme de la solution
    feNorm *normBDF2 = new feNorm(&U_surface, mesh, degQuad, funSol);
    feNorm *normDC3 = new feNorm(&U_surface, mesh, degQuad, funSol);
    std::vector<feNorm *> norms = {normBDF2, normDC3};
    // Systeme lineaire
    // feLinearSystemPETSc *linearSystem =
    // new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_MKL
    feLinearSystemMklPardiso *linearSystem =
      new feLinearSystemMklPardiso(formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
    // linearSystem->initialize();
    // Resolution
    feTolerances tol{1e-9, 1e-9, 20};
    std::string CodeIni = "BDF1/DCF"; // Define the way of initialization |"SolEx"->for exact
                                      // solution|  |"BDF1/DCF"->using only initial conditions|
    DC3Solver solver(tol, metaNumber, linearSystem, solDC3, norms, mesh, t0, t1, nTimeSteps,
                     CodeIni);
    solver.makeSteps(nTimeSteps, fespace);
    std::vector<double> &normL2BDF2 = solver.getNorm(0);
    std::vector<double> &normL2DC3 = solver.getNorm(1);
    normL2_BDF2[2 * iter] = *std::max_element(normL2BDF2.begin(), normL2BDF2.end());
    normL2_DC3[2 * iter] = *std::max_element(normL2DC3.begin(), normL2DC3.end());
#endif

    delete normBDF2;
    delete normDC3;
    delete linearSystem;
    delete masse_U_M2D;
    delete source_U_M2D;
    delete diff_U_M2D;
    delete solDC3;
    delete metaNumber;
    delete mesh;
#endif
  }
  delete funSource;
  delete funSol;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_BDF2[2 * i + 1] = log(normL2_BDF2[2 * (i - 1)] / normL2_BDF2[2 * i]) / log(2.);
    normL2_DC3[2 * i + 1] = log(normL2_DC3[2 * (i - 1)] / normL2_DC3[2 * i]) / log(2.);
  }
  printf("%12s \t %12s \t %12s \t %12s\t %12s \t %12s\n", "nSteps", "nElm", "||E-BDF||",
         "Taux BDF2", "||E-DC3||", "Taux DC3");
  for(int i = 0; i < nIter; ++i)

    printf("%12d \t %12d \t %12.6e \t %12.6e\t %12.6e \t %12.6e\n", TT[i], nElm[i],
           normL2_BDF2[2 * i], normL2_BDF2[2 * i + 1], normL2_DC3[2 * i], normL2_DC3[2 * i + 1]);
#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}
