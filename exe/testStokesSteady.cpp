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

double fSolU(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  return 2. * exp(x) * (x - 1.) * (x - 1.) * x * x * (y * y - y) * (2. * y - 1.);
  // return x * x * x * x * y * y * y * y;
}

double fSolV(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  return exp(x) * (x - 1.) * x * (x * x + 3. * x - 2.) * (y * y - 1.) * y * y;
  // return -4. / 5. * x * x * x * y * y * y * y * y;
}

double fSolP(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  return -424. + 156. * exp(1) +
         (y * y - y) *
           (-456. +
            exp(x) * (456. + x * x * (228. - 5. * (y * y - y)) + 2. * x * (y * y - y - 228.) +
                      2. * x * x * x * (y * y - y - 36.) + x * x * x * x * (y * y - y + 12.)));
  // return x * x * y * y;
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 0.0;
}

void fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par,
             std::vector<double> &res) {
  // double rho = par[0];
  double mu = par[1];
  // double nu = mu/rho;
  double x = pos[0];
  double y = pos[1];
  // double dudx = 2*exp(x)*x*y*(2*y*y-3*y+1)*(x*x*x + 2*x*x - 5*x + 2);
  // double dudy = 2*exp(x)*x*x*(x-1)*(x-1)*(6*y*y-6*y+1);
  // double dvdx = exp(x)*y*y*(y-1)*(y-1)*(x*x*x*x + 6*x*x*x + x*x -8*x+2); // Signe ? Pour avoir
  // div = 0 double dvdy = 2*exp(x)*x*y*(2*y*y-3*y+1)*(x*x*x + 2*x*x - 5*x + 2); double d2udx2 =
  // 2*exp(x)*y*(x*x*x*x + 6*x*x*x + x*x -8*x + 2) * (2*y*y - 3*y + 1); double d2udxy =
  // 2*exp(x)*x*(x*x*x + 2*x*x - 5*x + 2) * (6*y*y-6*y+1); double d2udy2 =
  // 12*exp(x)*x*x*(x-1)*(x-1)*(2*y-1); double d2vdx2 = -exp(x)*y*y*(x*x*x*x + 10*x*x*x + 19*x*x -
  // 6*x - 6)*(y-1)*(y-1); double d2vdxy = -2*exp(x)*y*(x*x*x*x + 6*x*x*x + x*x -8*x + 2) * (2*y*y -
  // 3*y + 1); double d2vdy2 = -2*exp(x)*x*(x-1)*(x*x+3*x-2)*(6*y*y-6*y+1); double ddxxdx = 2 *
  // d2udx2; double ddyxdy = d2udy2 + d2vdxy; double ddxydx = d2udxy + d2vdx2; double ddyydy = 2 *
  // d2vdy2;

  // double dpdx = (y*y - y)*(exp(x)*(4*x*x*x*(y*y - y + 12) + 6*x*x*(y*y - y - 36) + 2*(228 -
  // 5*(y*y - y)) + 2*(y*y - y - 228))
  // + exp(x)*(x*x*x*x*(y*y - y + 12) + 2*x*x*x*(y*y - y - 36) + x*x*(228 - 5*(y*y - y)) + 2*x*(y*y
  // - y - 228) + 456));

  // double dpdy = exp(x)*(y*y - y)*(x*x*x*x*(2*y - 1) + 2*x*x*x* (2*y - 1) - 5*x*x*(2*y - 1) +
  // 2*x*(2*y - 1))
  // + (2*y - 1)*(exp(x)*(x*x*x*x*(y*y - y + 12) + 2*x*x*x*(y*y - y - 36) + x*x*(228 - 5*(y*y - y))
  // + 2*x*(y*y - y - 228) + 456) - 456);

  // -div (nu*d) + grad(p)

  // res[0] = - (ddxxdx + ddyxdy) + dpdx;
  // res[1] = - (ddxydx + ddyydy) + dpdy;

  res[0] = -mu * (12 * (x * x * y * y * y * y + x * x * x * x * y * y)) + 2 * x * y * y;
  res[1] = -mu * (-8. / 5. * x * y * y * y * (10 * x * x + 3 * y * y)) + 2 * y * x * x;
}

int main(int argc, char **argv) {
#ifdef HAVE_PETSC
  petscInitialize(argc, argv);
#endif

  double rho = 1.0;
  double mu = 1.0;
  std::vector<double> stokesParam = {rho, mu};

  feFunction *funSolU = new feFunction(fSolU, {});
  feFunction *funSolV = new feFunction(fSolV, {});
  feFunction *funSolP = new feFunction(fSolP, {});
  feFunction *funZero = new feFunction(fZero, {});
  feVectorFunction *funSource = new feVectorFunction(fSource, stokesParam);

  int nIter = 4;
  std::vector<double> normL2_U(2 * nIter, 0.0);
  std::vector<double> normL2_V(2 * nIter, 0.0);
  std::vector<double> normL2_P(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);

  for(int iter = 0; iter < nIter; ++iter) {
    // std::string meshName = "../../data/squarePression" + std::to_string(iter+1) + ".msh";
    std::string meshName = "../../data/square" + std::to_string(iter + 1) + ".msh";
    // std::string meshName = "../../data/squareCoarse" + std::to_string(iter+1) + ".msh";

    feMesh2DP1 *mesh = new feMesh2DP1(meshName, false);
    nElm[iter] = mesh->getNbInteriorElems();

    feSpace1DP2 U_angle(mesh, "U", "Angle", funSolU);
    feSpace1DP2 U_haut(mesh, "U", "Haut", funSolU);
    feSpace1DP2 U_gauche(mesh, "U", "Gauche", funSolU);
    feSpaceTriP2 U_surface(mesh, "U", "Surface", funZero);

    feSpace1DP2 V_angle(mesh, "V", "Angle", funSolV);
    feSpace1DP2 V_haut(mesh, "V", "Haut", funSolV);
    feSpace1DP2 V_gauche(mesh, "V", "Gauche", funSolV);
    feSpaceTriP2 V_surface(mesh, "V", "Surface", funZero);

    // feSpace1DP2    P_angle(mesh, "P", "Angle",  funSolP);
    feSpace1DP1 P_haut(mesh, "P", "Haut", funSolP);
    // feSpace1DP2   P_gauche(mesh, "P", "Gauche", funSolP);
    feSpaceTriP1 P_surface(mesh, "P", "Surface", funZero);
    // feSpace1DP0 P_point(mesh, "P", "PointPression", funSolP);

    std::vector<feSpace *> fespace = {&U_angle, &U_haut, &U_gauche, &U_surface, &V_angle, &V_haut,
                                      &V_gauche, &V_surface,
                                      // &P_surface, &P_point};
                                      &P_surface, &P_haut};
    // &P_angle, &P_haut, &P_gauche, &P_surface};
    std::vector<feSpace *> feEssBC = {&U_angle, &U_haut,   &U_gauche, &V_angle,
                                      &V_haut,  &V_gauche, &P_haut};
    // &P_point};

    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // metaNumber->printNumberings();

    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);
    sol->initializeUnknowns(mesh, metaNumber);
    sol->initializeEssentialBC(mesh, metaNumber);

    // Formes (bi)lineaires
    int dQuad = 5;
    std::vector<feSpace *> spacesStokes2D = {&U_surface, &V_surface, &P_surface};

    feBilinearForm *stokes2D = new feBilinearForm(spacesStokes2D, mesh, dQuad,
                                                  new feSysElm_2D_Stokes(stokesParam, funSource));

    std::vector<feBilinearForm *> formMatrices = {stokes2D};
    std::vector<feBilinearForm *> formResiduals = {stokes2D};

    feNorm *normU = new feNorm(&U_surface, mesh, dQuad, funSolU);
    feNorm *normV = new feNorm(&V_surface, mesh, dQuad, funSolV);
    feNorm *normP = new feNorm(&P_surface, mesh, dQuad, funSolP);
    std::vector<feNorm *> norms = {normU, normV, normP};
    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
    linearSystem->initialize();
    linearSystem->assembleResiduals(sol);
    linearSystem->assembleMatrices(sol);
    feTolerances tol{1e-9, 1e-8, 10};
    solveStationary(&normL2_U[2 * iter], tol, metaNumber, linearSystem, formMatrices, formResiduals,
                    sol, norms, mesh);
    linearSystem->finalize();

    std::string vtkFile = "../../data/stokes1.vtk";
    feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);

    normU->computeL2Norm(metaNumber, sol, mesh);
    normV->computeL2Norm(metaNumber, sol, mesh);
    normP->computeL2Norm(metaNumber, sol, mesh);
    normL2_U[2 * iter] = normU->getNorm();
    normL2_V[2 * iter] = normV->getNorm();
    normL2_P[2 * iter] = normP->getNorm();

    delete normU;
    delete normV;
    delete normP;
    delete linearSystem;
    delete stokes2D;
    delete sol;
    delete metaNumber;
    delete mesh;
  }
  delete funSource;
  delete funSolU;
  delete funSolV;
  delete funSolP;
  delete funZero;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_U[2 * i + 1] =
      -log(normL2_U[2 * i] / normL2_U[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normL2_V[2 * i + 1] =
      -log(normL2_V[2 * i] / normL2_V[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normL2_P[2 * i + 1] =
      -log(normL2_P[2 * i] / normL2_P[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  }
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s\n", "nElm", "||u-uh||", "tauxU",
         "||v-vh||", "tauxV", "||p-ph||", "tauxP");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n", nElm[i],
           normL2_U[2 * i], normL2_U[2 * i + 1], normL2_V[2 * i], normL2_V[2 * i + 1],
           normL2_P[2 * i], normL2_P[2 * i + 1]);

#ifdef HAVE_PETSC
  petscFinalize();
#endif
  return 0;
}
