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

// double fSolU(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
//   double x = pos[0];
//   double y = pos[1];
//   double q = par[2];
//   double p = q -1.0;
//   return -pow(y,p) * (1.-pow(x,q))*pow(t,q);
// }

// double fSolV(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
//   double x = pos[0];
//   double y = pos[1];
//   double q = par[2];
//   double p = q -1.0;
//   return pow(x,p) * (1.-pow(y,q))*pow(t,q);
// }
// void fVecVeloc(const double t, const std::vector<double> &pos, const std::vector<double> &par,
//              std::vector<double> &res) {
//   double x = pos[0];
//   double y = pos[1];
//   double q = par[2];
//   double p = q -1.0;
//   res[0] = -pow(y,p) * (1.-pow(x,q))*pow(t,q);       //u
//   res[1] = pow(x,p) * (1.-pow(y,q))*pow(t,q);        //v
//   res[2] = pow(y,p)*q*pow(x,p)*pow(t,q);             //dudx
//   res[3] = -p*pow(y,p-1.) * (1-pow(x,q))*pow(t,q);   //dudy
//   res[4] = p*pow(x,p-1.) * (1-pow(y,q))*pow(t,q);    //dvdx
//   res[5] = -pow(x,p)*q*pow(y,p)*pow(t,q);            //dvdy

// }
// double fSolP(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
//   double x = pos[0];
//   double y = pos[1];
//   double q = par[2];
//   double p = q -1.0;
//   // return (2. +x +y );
//   return 0;
// }

// double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
//   return 0.0;
// }

// void fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par,
//              std::vector<double> &res) {
//   double mu = par[1];
//   double x = pos[0];
//   double y = pos[1];
//   double q = par[2];
//   double p = q -1.0;
//   double ft = pow(t,q) ;
//   double u = -pow(y,p) * (1-pow(x,q))*ft;
//   double v = pow(x,p) * (1-pow(y,q))*ft;
//   double uHat = -pow(y,p) * (1-pow(x,q));
//   double vHat = pow(x,p) * (1-pow(y,q)) ;

//   //dx
//   double dudx   = pow(y,p)*q*pow(x,p)*ft;
//   double d2udx2 = pow(y,p)*q*p*pow(x,p-1.)*ft;
//   double dvdx   = p*pow(x,p-1.) * (1-pow(y,q))*ft;
//   double d2vdx2 = p*(p-1.)*pow(x,p-2.) * (1-pow(y,q))*ft;
//   //dy
//   double dudy   = -p*pow(y,p-1.) * (1-pow(x,q))*ft;
//   double d2udy2 = -p*(p-1.)*pow(y,p-2.) * (1-pow(x,q))*ft;
//   double dvdy   = -pow(x,p)*q*pow(y,p)*ft;
//   double d2vdy2 = -pow(x,p)*q*p*pow(y,p-1.)*ft;

//   double dpdx = 0;
//   double dpdy = 0;
//   double rho = par[0];

//   res[0] = rho*q*pow(t,q-1.)*uHat + rho*(u*dudx + v* dudy)  - mu* (d2udx2 + d2udy2)+ dpdx;
//   res[1] = rho*q*pow(t,q-1.)*vHat + rho*(u*dvdx + v* dvdy) - mu*(d2vdx2 + d2vdy2 ) + dpdy;

// }
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
  res[3] = -p * pow(y, p - 1.) * (1 - pow(x, q)) * pow(t, a); // dudy
  res[4] = p * pow(x, p - 1.) * (1 - pow(y, q)) * pow(t, a); // dvdx
  res[5] = -pow(x, p) * q * pow(y, p) * pow(t, a); // dvdy
}
double fSolP(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q - 1.0;
  double a = par[3];
  // return (2. +x +y );
  return 0;
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
  double d2udx2 = pow(y, p) * q * p * pow(x, p - 1.) * ft;
  double dvdx = p * pow(x, p - 1.) * (1 - pow(y, q)) * ft;
  double d2vdx2 = p * (p - 1.) * pow(x, p - 2.) * (1 - pow(y, q)) * ft;
  // dy
  double dudy = -p * pow(y, p - 1.) * (1 - pow(x, q)) * ft;
  double d2udy2 = -p * (p - 1.) * pow(y, p - 2.) * (1 - pow(x, q)) * ft;
  double dvdy = -pow(x, p) * q * pow(y, p) * ft;
  double d2vdy2 = -pow(x, p) * q * p * pow(y, p - 1.) * ft;

  double dpdx = 0;
  double dpdy = 0;
  double rho = par[0];

  res[0] =
    rho * a * pow(t, a - 1.) * uHat + rho * (u * dudx + v * dudy) - mu * (d2udx2 + d2udy2) + dpdx;
  res[1] =
    rho * a * pow(t, a - 1.) * vHat + rho * (u * dvdx + v * dvdy) - mu * (d2vdx2 + d2vdy2) + dpdy;
}

int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double rho = 1.0;
  double mu = 1.0;
  double q = 2.;
  double exposant = 5.;
  std::vector<double> stokesParam = {rho, mu, q, exposant};

  feFunction *funSolU = new feFunction(fSolU, stokesParam);
  feFunction *funSolV = new feFunction(fSolV, stokesParam);
  feFunction *funSolP = new feFunction(fSolP, stokesParam);
  feFunction *funZero = new feFunction(fZero, {});
  feVectorFunction *funSource = new feVectorFunction(fSource, stokesParam);
  feVectorFunction *funVecVeloc = new feVectorFunction(fVecVeloc, stokesParam);

  int nIter = 3;
  std::vector<double> normU_BDF1(2 * nIter, 0.0);
  std::vector<double> normU_DC2F(2 * nIter, 0.0);
  std::vector<double> normU_DC3F(2 * nIter, 0.0);
  std::vector<double> normP_BDF1(2 * nIter, 0.0);
  std::vector<double> normP_DC2F(2 * nIter, 0.0);
  std::vector<double> normP_DC3F(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);

  for(int iter = 0; iter < nIter; ++iter) {
    std::string meshName = "../../data/Square/squareNS" + std::to_string(iter + 1) + ".msh";

    feMesh2DP1 *mesh = new feMesh2DP1(meshName, false);
    nElm[iter] = mesh->getNbInteriorElems();

    feSpace1DP3 U_angle = feSpace1DP3(mesh, "U", "Angle", funSolU);
    feSpace1DP3 U_haut = feSpace1DP3(mesh, "U", "Haut", funSolU);
    feSpace1DP3 U_gauche = feSpace1DP3(mesh, "U", "Gauche", funSolU);
    feSpaceTriP3 U_surface = feSpaceTriP3(mesh, "U", "Surface", funSolU);

    feSpace1DP3 V_angle = feSpace1DP3(mesh, "V", "Angle", funSolV);
    feSpace1DP3 V_haut = feSpace1DP3(mesh, "V", "Haut", funSolV);
    feSpace1DP3 V_gauche = feSpace1DP3(mesh, "V", "Gauche", funSolV);
    feSpaceTriP3 V_surface = feSpaceTriP3(mesh, "V", "Surface", funSolV);

    feSpace1DP0 P_point = feSpace1DP0(mesh, "P", "PointPression", funSolP);
    feSpaceTriP2 P_surface = feSpaceTriP2(mesh, "P", "Surface", funSolP);

    std::vector<feSpace *> fespace = {&U_angle, &U_haut,   &U_gauche,  &U_surface, &V_angle,
                                      &V_haut,  &V_gauche, &V_surface, &P_point,   &P_surface};

    std::vector<feSpace *> feEssBC = {&U_angle, &U_haut,   &U_gauche, &V_angle,
                                      &V_haut,  &V_gauche, &P_point};

    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);

    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 5 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);

    // Formes (bi)lineaires
    int nQuad = 16;
    std::vector<feSpace *> spacesNS2D = {&U_surface, &V_surface, &P_surface};

    // feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, nQuad, new
    // feSysElm_2D_NavierStokes(stokesParam, funSource));
    feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, nQuad,
                                              new feSysElm_2D_NavierStokes(stokesParam, funSource));

    std::vector<feBilinearForm *> formMatrices = {NS2D};
    std::vector<feBilinearForm *> formResiduals = {NS2D};

    std::vector<feSpace *> fespaceNorm = {&U_surface, &V_surface};
    feNorm *normVecVeloc_BDF1 =
      new feNorm(fespaceNorm, mesh, nQuad, funZero, funVecVeloc, "NormL2");
    feNorm *normVecVeloc_DC2F =
      new feNorm(fespaceNorm, mesh, nQuad, funZero, funVecVeloc, "NormL2");
    feNorm *normVecVeloc_DC3F =
      new feNorm(fespaceNorm, mesh, nQuad, funZero, funVecVeloc, "NormL2");
    feNorm *normPression_BDF1 = new feNorm(&P_surface, mesh, nQuad, funSolP);
    feNorm *normPression_DC2F = new feNorm(&P_surface, mesh, nQuad, funSolP);
    feNorm *normPression_DC3F = new feNorm(&P_surface, mesh, nQuad, funSolP);
    std::vector<feNorm *> norms = {normVecVeloc_BDF1, normVecVeloc_DC2F, normVecVeloc_DC3F,
                                   normPression_BDF1, normPression_DC2F, normPression_DC3F};
    // Systeme lineaire

    // feLinearSystemPETSc *linearSystem =
    // new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_MKL
    feLinearSystemMklPardiso *linearSystem =
      new feLinearSystemMklPardiso(formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
    // Resolution
    feTolerances tol{1e-6, 1e-9, 20};

    DC3FSolver solver(tol, metaNumber, linearSystem, sol, norms, mesh, t0, t1, nTimeSteps);
    solver.makeSteps(nTimeSteps, fespace);

    // maxNormL2DC3_U[2 * iter] = fmax(maxNormL2DC3_U[2 * iter], normL2DC3[3 * i]);
    // maxNormL2DC3_V[2 * iter] = fmax(maxNormL2DC3_V[2 * iter], normL2DC3[3 * i + 1]);
    // maxNormL2DC3_P[2 * iter] = fmax(maxNormL2DC3_P[2 * iter], normL2DC3[3 * i + 2]);
    // std::string vtkFile = "../../data/solutionManufacturee" + std::to_string(iter+1) + ".vtk";
    // std::string vtkFile = "../../data/CasTest_NS_steady" + std::to_string(iter+1) + ".vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);

    // normU->computeL2Norm(metaNumber, sol, mesh);
    // normV->computeL2Norm(metaNumber, sol, mesh);
    // normP->computeL2Norm(metaNumber, sol, mesh);
    normU_BDF1[2 * iter] = normVecVeloc_BDF1->getNorm();
    normU_DC2F[2 * iter] = normVecVeloc_DC2F->getNorm();
    normU_DC3F[2 * iter] = normVecVeloc_DC3F->getNorm();
    normP_BDF1[2 * iter] = normPression_BDF1->getNorm();
    normP_DC2F[2 * iter] = normPression_DC2F->getNorm();
    normP_DC3F[2 * iter] = normPression_DC3F->getNorm();
#endif
    // std::string vtkFile = "../../data/TestNS_DC3.vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);
    delete normVecVeloc_BDF1;
    delete normVecVeloc_DC2F;
    delete normVecVeloc_DC3F;
    delete normPression_BDF1;
    delete normPression_DC2F;
    delete normPression_DC3F;
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
  for(int i = 1; i < nIter; ++i) {
    normU_BDF1[2 * i + 1] =
      -log(normU_BDF1[2 * i] / normU_BDF1[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normU_DC2F[2 * i + 1] =
      -log(normU_DC2F[2 * i] / normU_DC2F[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normU_DC3F[2 * i + 1] =
      -log(normU_DC3F[2 * i] / normU_DC3F[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normP_BDF1[2 * i + 1] =
      -log(normP_BDF1[2 * i] / normP_BDF1[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normP_DC2F[2 * i + 1] =
      -log(normP_DC2F[2 * i] / normP_DC2F[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normP_DC3F[2 * i + 1] =
      -log(normP_DC3F[2 * i] / normP_DC3F[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  }
  printf("\n");
  printf("\n");
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t \n", "nSteps", "nElm",
         "||u-uh||_BDF1", "tauxU_BDF1", "||u-uh||_DC2F", "tauxU_DC2F", "||u-uh||_DC3F",
         "tauxU_DC3F");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \n", TT[i],
           nElm[i], normU_BDF1[2 * i], normU_BDF1[2 * i + 1], normU_DC2F[2 * i],
           normU_DC2F[2 * i + 1], normU_DC3F[2 * i], normU_DC3F[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t \n", "nSteps", "nElm",
         "||p-ph||_BDF1", "tauxP_BDF1", "||p-ph||_DC2F", "tauxP_DC2F", "||p-ph||_DC3F",
         "tauxP_DC3F");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \n", TT[i],
           nElm[i], normP_BDF1[2 * i], normP_BDF1[2 * i + 1], normP_DC2F[2 * i],
           normP_DC2F[2 * i + 1], normP_DC3F[2 * i], normP_DC3F[2 * i + 1]);

#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}
