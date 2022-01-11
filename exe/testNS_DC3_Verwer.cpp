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
#include "feExporter.h"

#ifdef HAVE_PETSC
#include "feLinearSystemPETSc.h"
#endif
#ifdef HAVE_MKL
#include "feLinearSystemMklPardiso.h"
#endif

double fSolU(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q - 1.0;
  double a = par[3];
  return -pow(y, p) * (1. - pow(x, q)) * pow(t, a);
}

double fSolV(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q - 1.0;
  double a = par[3];
  return pow(x, p) * (1. - pow(y, q)) * pow(t, a);
}
void fVecVeloc(const double t, const std::vector<double> &pos, const std::vector<double> &par,
               std::vector<double> &res) {
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q - 1.0;
  double a = par[3];
  res[0] = -pow(y, p) * (1. - pow(x, q)) * pow(t, a); // u
  res[1] = pow(x, p) * (1. - pow(y, q)) * pow(t, a); // v
  res[2] = pow(y, p) * q * pow(x, p) * pow(t, a); // dudx
  res[3] = (q == 1.) ? 0. : -p * pow(y, p - 1.) * (1 - pow(x, q)) * pow(t, a); // dudy
  res[4] = (q == 1.) ? 0. : p * pow(x, p - 1.) * (1 - pow(y, q)) * pow(t, a); // dvdx
  res[5] = -pow(x, p) * q * pow(y, p) * pow(t, a); // dvdy
}
double fSolP(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q - 1.0;
  double a = par[3];
  return (2. + x + y);
  // return 0;
}
double fUDot(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q - 1.0;
  double a = par[3];
  double uHat = -pow(y, p) * (1 - pow(x, q));

  return a * pow(t, a - 1.) * uHat;
}

double fVDot(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q - 1.0;
  double a = par[3];

  double vHat = pow(x, p) * (1 - pow(y, q));

  return a * pow(t, a - 1.) * vHat;
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 0.0;
}

void fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par,
             std::vector<double> &res) {
  double mu = par[1];
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q - 1.0;
  double a = par[3];
  double ft = pow(t, a);
  double u = -pow(y, p) * (1 - pow(x, q)) * ft;
  double v = pow(x, p) * (1 - pow(y, q)) * ft;
  double uHat = -pow(y, p) * (1 - pow(x, q));
  double vHat = pow(x, p) * (1 - pow(y, q));

  // dx
  double dudx = pow(y, p) * q * pow(x, p) * ft;
  double d2udx2 = (q == 1.) ? 0. : pow(y, p) * q * p * pow(x, p - 1.) * ft;
  double dvdx = (q == 1.) ? 0. : p * pow(x, p - 1.) * (1 - pow(y, q)) * ft;
  double d2vdx2 = (q == 1.) ? 0. : p * (p - 1.) * pow(x, p - 2.) * (1 - pow(y, q)) * ft;
  // dy
  double dudy = (q == 1.) ? 0. : -p * pow(y, p - 1.) * (1 - pow(x, q)) * ft;
  double d2udy2 = (q == 1.) ? 0. : -p * (p - 1.) * pow(y, p - 2.) * (1 - pow(x, q)) * ft;
  double dvdy = -pow(x, p) * q * pow(y, p) * ft;
  double d2vdy2 = (q == 1.) ? 0. : -pow(x, p) * q * p * pow(y, p - 1.) * ft;

  double dpdx = 1;
  double dpdy = 1;
  double rho = par[0];

  res[0] =
    rho * a * pow(t, a - 1.) * uHat + rho * (u * dudx + v * dudy) - mu * (d2udx2 + d2udy2) + dpdx;
  res[1] =
    rho * a * pow(t, a - 1.) * vHat + rho * (u * dvdx + v * dvdy) - mu * (d2vdx2 + d2vdy2) + dpdy;
}
void fSolDot(const double t, const std::vector<double> &pos, const std::vector<double> &par,
             std::vector<double> &res) {
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q - 1.0;
  double a = par[3];

  double uHat = -pow(y, p) * (1 - pow(x, q));
  double vHat = pow(x, p) * (1 - pow(y, q));

  res[0] = a * pow(t, a - 1.) * uHat;
  res[1] = a * pow(t, a - 1.) * vHat;
}

void fgrad(const double t, const std::vector<double> &pos, const std::vector<double> &par,
           std::vector<double> &res) {
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q - 1.0;
  double a = par[3];
  double ft = pow(t, a);

  res[0] = pow(y, p) * q * pow(x, p) * ft;
  res[1] = -p * pow(y, p - 1.) * (1 - pow(x, q)) * ft;
}

int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double rho = 1.0;
  double mu = 1.0;
  double q = 2.;
  double exposant = 1.;
  std::vector<double> stokesParam = {rho, mu, q, exposant};

  feFunction *funSolU = new feFunction(fSolU, stokesParam);
  feFunction *funSolV = new feFunction(fSolV, stokesParam);
  feFunction *funSolP = new feFunction(fSolP, stokesParam);
  feFunction *funZero = new feFunction(fZero, {});
  // feVectorFunction *funSolDot = new feVectorFunction(fSolDot, stokesParam);
  feFunction *funUDot = new feFunction(fUDot, stokesParam);
  feFunction *funVDot = new feFunction(fVDot, stokesParam);
  feVectorFunction *funSource = new feVectorFunction(fSource, stokesParam);
  feVectorFunction *funVecVeloc = new feVectorFunction(fVecVeloc, stokesParam);
  feVectorFunction *funGrad = new feVectorFunction(fgrad, stokesParam);

  int nIter = 3;
  std::vector<double> normU_BDF2(2 * nIter, 0.0);
  std::vector<double> normU_DC3(2 * nIter, 0.0);
  std::vector<double> normP_BDF2(2 * nIter, 0.0);
  std::vector<double> normP_DC3(2 * nIter, 0.0);
  std::vector<double> normL_BDF2(2 * nIter, 0.0);
  std::vector<double> normL_DC3(2 * nIter, 0.0);
  std::vector<double> intL_BDF2(2 * nIter, 0.0);
  std::vector<double> intL_DC3(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);

  for(int iter = 0; iter < nIter; ++iter) {
    std::string meshName =
      "../../data/Square/SquareNS_1Edge/squareNS" + std::to_string(iter + 1) + "_1Edge.msh";
    // std::string meshName = "../../data/Square/squareTest" + std::to_string(iter + 1) + ".msh";

    feMesh2DP1 *mesh = new feMesh2DP1(meshName, false);
    nElm[iter] = mesh->getNbInteriorElems();

    feSpace1DP3 U_Bord = feSpace1DP3(mesh, "U", "Bord", funSolU);
    feSpaceTriP3 U_surface = feSpaceTriP3(mesh, "U", "Surface", funSolU);

    feSpace1DP3 V_Bord = feSpace1DP3(mesh, "V", "Bord", funSolV);
    feSpaceTriP3 V_surface = feSpaceTriP3(mesh, "V", "Surface", funSolV);

    feSpace1DP0 P_point = feSpace1DP0(mesh, "P", "PointPression", funSolP);
    feSpaceTriP2 P_surface = feSpaceTriP2(mesh, "P", "Surface", funSolP);

    feSpace1DP3 Lu_Bord = feSpace1DP3(mesh, "Lu", "Bord", funSolU);
    feSpace1DP3 Lv_Bord = feSpace1DP3(mesh, "Lv", "Bord", funSolV);

    feSpace1DP3 Bu_Bord = feSpace1DP3(mesh, "Bu", "Bord", funSolU);
    feSpace1DP3 Bv_Bord = feSpace1DP3(mesh, "Bv", "Bord", funSolV);

    std::vector<feSpace *> fespace = {&U_Bord,  &U_surface, &V_Bord,  &V_surface, &P_surface,
                                      &Lu_Bord, &Lv_Bord,   &Bu_Bord, &Bv_Bord};

    std::vector<feSpace *> feEssBC = {&P_point};

    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);

    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 5 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);

    // Formes (bi)lineaires
    int nQuad = 16;
    std::vector<feSpace *> spacesNS2D = {&U_surface, &V_surface, &P_surface};
    std::vector<feSpace *> spaceWeak_bordU = {&U_Bord, &Lu_Bord, &Bu_Bord};
    std::vector<feSpace *> spaceWeak_bordV = {&V_Bord, &Lv_Bord, &Bv_Bord};

    // feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, nQuad, new
    // feSysElm_2D_NavierStokes(stokesParam, funSource));
    feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, nQuad,
                                              new feSysElm_2D_NavierStokes(stokesParam, funSource));
    feBilinearForm *weakBC_bordU =
      new feBilinearForm(spaceWeak_bordU, mesh, nQuad, new feSysElm_1D_weakBC_edo1(1.0, funUDot));
    feBilinearForm *weakBC_bordV =
      new feBilinearForm(spaceWeak_bordV, mesh, nQuad, new feSysElm_1D_weakBC_edo1(1.0, funVDot));

    std::vector<feBilinearForm *> formMatrices = {NS2D, weakBC_bordU, weakBC_bordV};
    std::vector<feBilinearForm *> formResiduals = {NS2D, weakBC_bordU, weakBC_bordV};

    std::vector<feSpace *> fespaceNorm = {&U_surface, &V_surface};
    feNorm *normVecVeloc_BDF2 =
      new feNorm(fespaceNorm, mesh, nQuad, funZero, funVecVeloc, "NormL2");
    feNorm *normVecVeloc_DC3 = new feNorm(fespaceNorm, mesh, nQuad, funZero, funVecVeloc, "NormL2");
    feNorm *normPression_BDF2 = new feNorm(&P_surface, mesh, nQuad, funSolP);
    feNorm *normPression_DC3 = new feNorm(&P_surface, mesh, nQuad, funSolP);
    feNorm *normLambda_BDF2 = new feNorm(&Lu_Bord, mesh, nQuad, funSolP, funGrad, "IntForce");
    feNorm *normLambda_DC3 = new feNorm(&Lu_Bord, mesh, nQuad, funSolP, funGrad, "IntForce");
    feNorm *IntLambda_BDF2 = new feNorm(&Lu_Bord, mesh, nQuad, funSolP, funGrad, "IntLambda");
    feNorm *IntLambda_DC3 = new feNorm(&Lu_Bord, mesh, nQuad, funSolP, funGrad, "IntLambda");
    // Group by type of fields, not type of solver
    std::vector<feNorm *> norms = {normVecVeloc_BDF2, normVecVeloc_DC3, normPression_BDF2,
                                   normPression_DC3,  normLambda_BDF2,  normLambda_DC3,
                                   IntLambda_BDF2,    IntLambda_DC3};
    // Systeme lineaire

    // feLinearSystemPETSc *linearSystem =
    // new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_MKL
    feLinearSystemMklPardiso *linearSystem =
      new feLinearSystemMklPardiso(formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
    // Resolution
    feTolerances tol{1e-9, 1e-9, 20};
    std::string CodeIni = "BDF1/DCF"; // Define the way of initialization |"SolEx"->for exact
                                      // solution|  |"BDF1/DCF"->using only initial conditions|
    DC3Solver solver(tol, metaNumber, linearSystem, sol, norms, mesh, t0, t1, nTimeSteps, CodeIni);
    solver.makeSteps(nTimeSteps, fespace);

    // std::string vtkFile = "../../data/solutionManufacturee" + std::to_string(iter+1) + ".vtk";
    // std::string vtkFile = "../../data/CasTest_NS_steady" + std::to_string(iter+1) + ".vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);

    std::vector<double> &normUBDF2 = solver.getNorm(0);
    std::vector<double> &normUDC3 = solver.getNorm(1);
    std::vector<double> &normPBDF2 = solver.getNorm(2);
    std::vector<double> &normPDC3 = solver.getNorm(3);
    std::vector<double> &normLBDF2 = solver.getNorm(4);
    std::vector<double> &normLDC3 = solver.getNorm(5);
    std::vector<double> &IntLBDF2 = solver.getNorm(6);
    std::vector<double> &IntLDC3 = solver.getNorm(7);

    normLDC3[1] = 0.;
    normLDC3[2] = 0.;
    IntLDC3[1] = 0.;
    IntLDC3[2] = 0.;

    for(int i = 0; i < normLBDF2.size(); i++) {
      // std::cout<<normLBDF2[i]<<std::endl;
      normLBDF2[i] = fabs(normLBDF2[i]);
    }
    for(int i = 0; i < normLDC3.size(); i++) normLDC3[i] = fabs(normLDC3[i]);
    for(int i = 0; i < IntLBDF2.size(); i++) {
      // std::cout<<IntLBDF2[i]<<std::endl;
      IntLBDF2[i] = fabs(IntLBDF2[i]);
    }
    for(int i = 0; i < IntLDC3.size(); i++) IntLDC3[i] = fabs(IntLDC3[i]);
    for(int i = 0; i < normUBDF2.size(); i++) {
      std::cout << "U_BDF2  " << normUBDF2[i] << std::endl;
    }
    for(int i = 0; i < normUDC3.size(); i++) { std::cout << "U_DC3  " << normUDC3[i] << std::endl; }
    for(int i = 0; i < normLBDF2.size(); i++) {
      std::cout << "L_BDF2  " << normLBDF2[i] << std::endl;
    }
    for(int i = 0; i < normLDC3.size(); i++) { std::cout << "L_DC3  " << normLDC3[i] << std::endl; }

    normU_BDF2[2 * iter] = *std::max_element(normUBDF2.begin(), normUBDF2.end());
    normU_DC3[2 * iter] = *std::max_element(normUDC3.begin(), normUDC3.end());
    normP_BDF2[2 * iter] = *std::max_element(normPBDF2.begin(), normPBDF2.end());
    normP_DC3[2 * iter] = *std::max_element(normPDC3.begin(), normPDC3.end());
    normL_BDF2[2 * iter] = *std::max_element(normLBDF2.begin(), normLBDF2.end());
    normL_DC3[2 * iter] = *std::max_element(normLDC3.begin(), normLDC3.end());
    intL_BDF2[2 * iter] = *std::max_element(IntLBDF2.begin(), IntLBDF2.end());
    intL_DC3[2 * iter] = *std::max_element(IntLDC3.begin(), IntLDC3.end());

#endif
    // std::string vtkFile = "../../data/TestNS_DC3.vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);
    delete normVecVeloc_BDF2;
    delete normVecVeloc_DC3;
    delete normPression_BDF2;
    delete normPression_DC3;
    delete linearSystem;
    delete NS2D;
    delete sol;
    delete metaNumber;
    delete mesh;
#endif
  }
  delete funSource;
  delete funSolU;
  delete funSolV;
  delete funSolP;
  delete funZero;

  // Calcul du taux de convergence
  // for(int i = 0; i < normL_BDF2.size(); i++)
  // normL_BDF2[i] = fabs(normL_BDF2[i]);
  //   for(int i = 0; i < normL_DC3.size(); i++)
  //     normL_DC3[i] = fabs(normL_DC3[i]);
  for(int i = 1; i < nIter; ++i) {
    normU_BDF2[2 * i + 1] = -log(normU_BDF2[2 * i] / normU_BDF2[2 * (i - 1)]) / log(2.);
    normU_DC3[2 * i + 1] = -log(normU_DC3[2 * i] / normU_DC3[2 * (i - 1)]) / log(2.);
    normP_BDF2[2 * i + 1] = -log(normP_BDF2[2 * i] / normP_BDF2[2 * (i - 1)]) / log(2.);
    normP_DC3[2 * i + 1] = -log(normP_DC3[2 * i] / normP_DC3[2 * (i - 1)]) / log(2.);
    normL_BDF2[2 * i + 1] = -log(normL_BDF2[2 * i] / normL_BDF2[2 * (i - 1)]) / log(2.);
    normL_DC3[2 * i + 1] = -log(normL_DC3[2 * i] / normL_DC3[2 * (i - 1)]) / log(2.);
  }
  printf("\n");
  printf("\n");
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t \n", "nSteps", "nElm", "||u-uh||_BDF2",
         "tauxU_BDF2", "||u-uh||_DC3", "tauxU_DC3");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \n", TT[i], nElm[i],
           normU_BDF2[2 * i], normU_BDF2[2 * i + 1], normU_DC3[2 * i], normU_DC3[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \n", "nSteps", "nElm", "||p-ph||_BDF2",
         "tauxP_BDF2", "||p-ph||_DC3", "tauxP_DC3");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \n", TT[i], nElm[i],
           normP_BDF2[2 * i], normP_BDF2[2 * i + 1], normP_DC3[2 * i], normP_DC3[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \n", "nSteps", "nElm", "||L-Lh||_BDF2",
         "tauxL_BDF2", "||L-Lh||_DC3", "tauxL_DC3");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \n", TT[i], nElm[i],
           normL_BDF2[2 * i], normL_BDF2[2 * i + 1], normL_DC3[2 * i], normL_DC3[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("%12s \t %12s \t %12s \t %12s \n", "nSteps", "nElm", "||L||_BDF2", "||L||_DC3");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e \n", TT[i], nElm[i], intL_BDF2[2 * i],
           intL_DC3[2 * i]);

#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}
