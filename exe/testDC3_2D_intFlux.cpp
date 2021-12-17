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
#include "feExporter.h"
#include "feLinearSystemMklPardiso.h"
#define USING_PETSC
#ifdef HAVE_PETSC
#include "feLinearSystemPETSc.h"
#endif
#ifdef HAVE_MKL
#include "feLinearSystemMklPardiso.h"
#endif

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double y1 = x[1];
  int a = par[1];
  return pow(t, a) * (x1 * x1 + 1) * y1;
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double y1 = x[1];
  double c1 = par[0];
  double a = par[1];
  double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
  return -beta * (x1 * x1 + 1) * y1 + c1 * 2 * pow(t, a) * y1;
}

double dfdy(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double y1 = x[1];
  double c1 = par[0];
  double a = par[1];
  return pow(t, a) * (x1 * x1 + 1);
}

double dfdx(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double y1 = x[1];
  double c1 = par[0];
  double a = par[1];
  return pow(t, a) * y1 * (2 * x1);
}

double fSolDot(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double c1 = par[0];
  double a = par[1];
  double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
  double x1 = x[0];
  double y1 = x[1];
  return beta * (x1 * x1 + 1) * y1;
}

double f0(const double t, const std::vector<double> &x, const std::vector<double> par) { return 0; }

void fgrad(const double t, const std::vector<double> &x, const std::vector<double> &par,
           std::vector<double> &res) {
  double x1 = x[0];
  double y1 = x[1];
  double c1 = par[0];
  double a = par[1];
  res[0] = pow(t, a) * y1 * (2 * x1);
  res[1] = pow(t, a) * (x1 * x1 + 1);
}

int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double exposant = 4.;
  double kd = 0.1;
  std::vector<double> par = {kd, exposant};
  feFunction *funSol = new feFunction(fSol, par);
  feFunction *funSolDot = new feFunction(fSolDot, par);
  feFunction *funSource = new feFunction(fSource, par);
  feFunction *dfunSoldx = new feFunction(dfdx, par);
  feFunction *dfunSoldy = new feFunction(dfdy, par);
  feFunction *fun0 = new feFunction(f0, par);
  feVectorFunction *fungrad = new feVectorFunction(fgrad, par);

  int nIter = 4;
  std::vector<double> normL2_BDF2(2 * nIter, 0.0);
  std::vector<double> normL2_DC3(2 * nIter, 0.0);
  std::vector<double> normBDF2_bord(2 * nIter, 0.0);
  std::vector<double> normDC3_bord(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);

  for(int iter = 0; iter < nIter; ++iter) {
    // std::string meshName = "../../data/circle.msh";
    // std::string meshName = "../../data/squareTest2.msh";
    std::string meshName = "../../data/Square/square_1edge_" + std::to_string(iter + 1) + ".msh";
    // Maillage
    feMesh2DP1 *mesh = new feMesh2DP1(meshName, false);
    nElm[iter] = mesh->getNbInteriorElems();
    // Espaces d'interpolation

    feSpace1DP3 U_bord = feSpace1DP3(mesh, "U", "Bord", fun0);
    feSpaceTriP3 U_surface = feSpaceTriP3(mesh, "U", "Surface", funSol);

    feSpace1DP3 L_bord = feSpace1DP3(mesh, "L", "Bord", funSol);
    feSpace1DP3 V_bord = feSpace1DP3(mesh, "V", "Bord", funSol);

    std::vector<feSpace *> fespace = {&U_bord, &U_surface, &L_bord, &V_bord};
    std::vector<feSpace *> feEssBC = {};
    // Numerotations
    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // Solution
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 10 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feSolution *solDC3 = new feSolution(mesh, fespace, feEssBC, metaNumber);
    // solDC3->initializeTemporalSolution(t0, t1, nTimeSteps);
    // solDC3->initializeUnknowns(mesh, metaNumber);
    // solDC3->initializeEssentialBC(mesh, metaNumber);
    // std::string vtkFile = "../../data/TestDC3avantIni2.vtk";
    // feExporterVTK writer1(vtkFile, mesh, solDC3, metaNumber, fespace);
    // Formes (bi)lineaires
    int degQuad = 11;
    std::vector<feSpace *> spaceDiffusion2D_U = {&U_surface};
    std::vector<feSpace *> spaceSource2D_U = {&U_surface};
    std::vector<feSpace *> spaceMasse2D_U = {&U_surface};
    std::vector<feSpace *> spaceWeak_bord = {&U_bord, &L_bord, &V_bord};

    feBilinearForm *diff_U_M2D =
      new feBilinearForm(spaceDiffusion2D_U, mesh, degQuad, new feSysElm_2D_Diffusion(kd, nullptr));
    feBilinearForm *source_U_M2D =
      new feBilinearForm(spaceSource2D_U, mesh, degQuad, new feSysElm_2D_Source(1.0, funSource));
    feBilinearForm *masse_U_M2D =
      new feBilinearForm(spaceMasse2D_U, mesh, degQuad, new feSysElm_2D_Masse(1.0, nullptr));
    feBilinearForm *weakBC_bord = new feBilinearForm(spaceWeak_bord, mesh, degQuad,
                                                     new feSysElm_1D_weakBC_edo1(1.0, funSolDot));

    std::vector<feBilinearForm *> formMatrices = {diff_U_M2D, masse_U_M2D, weakBC_bord};
    std::vector<feBilinearForm *> formResiduals = {diff_U_M2D, masse_U_M2D, source_U_M2D,
                                                   weakBC_bord};

    // Norme de la solution
    feNorm *normBDF2 = new feNorm(&U_surface, mesh, degQuad, funSol);
    feNorm *normDC3 = new feNorm(&U_surface, mesh, degQuad, funSol);
    // Choose "NormeL2" to calculate the L2 Norm and "IntFlux" to calculate the integrale of the
    // flux
    feNorm *int_Lambda_BDF2 = new feNorm(&L_bord, mesh, degQuad, funSol, fungrad, "IntFlux");
    feNorm *int_Lambda_DC3 = new feNorm(&L_bord, mesh, degQuad, funSol, fungrad, "IntFlux");
    std::vector<feNorm *> norms = {normBDF2, normDC3, int_Lambda_BDF2, int_Lambda_DC3};

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
    std::vector<double> &normL2_BDF2_Lambda_bord = solver.getInt(2);
    std::vector<double> &normL2_DC3_Lambda_bord = solver.getInt(3);
    for(int i = 0; i < normL2_BDF2_Lambda_bord.size(); i++)
      normL2_BDF2_Lambda_bord[i] = fabs(normL2_BDF2_Lambda_bord[i]);
    for(int i = 0; i < normL2_DC3_Lambda_bord.size(); i++)
      normL2_DC3_Lambda_bord[i] = fabs(normL2_DC3_Lambda_bord[i]);
    normL2_BDF2[2 * iter] = *std::max_element(normL2BDF2.begin(), normL2BDF2.end());
    normL2_DC3[2 * iter] = *std::max_element(normL2DC3.begin(), normL2DC3.end());
    normBDF2_bord[2 * iter] =
      *std::max_element(normL2_BDF2_Lambda_bord.begin(), normL2_BDF2_Lambda_bord.end());
    normBDF2_bord[2 * iter] =
      *std::max_element(normL2_BDF2_Lambda_bord.begin(), normL2_BDF2_Lambda_bord.end());
    normDC3_bord[2 * iter] =
      *std::max_element(normL2_DC3_Lambda_bord.begin(), normL2_DC3_Lambda_bord.end());
    normDC3_bord[2 * iter] =
      *std::max_element(normL2_DC3_Lambda_bord.begin(), normL2_DC3_Lambda_bord.end());
#endif
    std::string vtkFile = "../../data/TestDC3.vtk";
    feExporterVTK writer(vtkFile, mesh, solDC3, metaNumber, fespace);
    delete normBDF2;
    delete normDC3;
    delete int_Lambda_BDF2;
    delete int_Lambda_DC3;
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
    normBDF2_bord[2 * i + 1] = log(normBDF2_bord[2 * (i - 1)] / normBDF2_bord[2 * i]) / log(2.);
    normBDF2_bord[2 * i + 1] = log(normBDF2_bord[2 * (i - 1)] / normBDF2_bord[2 * i]) / log(2.);
    normDC3_bord[2 * i + 1] = log(normDC3_bord[2 * (i - 1)] / normDC3_bord[2 * i]) / log(2.);
    normDC3_bord[2 * i + 1] = log(normDC3_bord[2 * (i - 1)] / normDC3_bord[2 * i]) / log(2.);
  }
  printf("%12s \t %12s \t %12s \t %12s\t %12s \t %12s\n", "nSteps", "nElm", "||E-BDF||",
         "Taux BDF2", "||E-DC3||", "Taux DC3");
  for(int i = 0; i < nIter; ++i)

    printf("%12d \t %12d \t %12.6e \t %12.6e\t %12.6e \t %12.6e\n", TT[i], nElm[i],
           normL2_BDF2[2 * i], normL2_BDF2[2 * i + 1], normL2_DC3[2 * i], normL2_DC3[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("Lambda based on BDF2");
  printf("\n");
  printf("%12s \t %12s \t %16s \t %16s \n", "nTimeSteps", "nElm", "||E_lambda||", "Taux lambda");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %16.6e \t %16.6e \n", TT[i], nElm[i], normBDF2_bord[2 * i],
           normBDF2_bord[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("Lambda based on DC3");
  printf("\n");
  printf("%12s \t %12s \t %16s \t %16s \n", "nTimeSteps", "nElm", "||E_lambda||", "Taux lambda");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %16.6e \t %16.6e \n", TT[i], nElm[i], normDC3_bord[2 * i],
           normDC3_bord[2 * i + 1]);
#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}