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
#include "feCurvedAdapt.h"

double fSolU(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double F = exp(-2.0 * t);
  return F * sin(x) * cos(y); // Taylor-Green
  // return pow(x,4) + pow(y,4);
}

double fSolV(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double F = exp(-2.0 * t);
  return -F * cos(x) * sin(y); // Taylor-Green
  // return -4.0 * pow(x,3) * y;
}

double fSolP(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  // double rho = par[0];
  double F = exp(-2.0 * t);
  return F * F * 1.0 / 4.0 * (cos(2. * x) + cos(2. * y)); // Taylor-Green
  // return pow(x,4) * pow(y,4);
}

double fOne(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 1.0;
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 0.0;
}

void fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par,
             std::vector<double> &res) {
  // double rho = par[0];
  // double mu = par[1];
  // double nu = mu/rho;
  double x = pos[0];
  double y = pos[1];

  // res[0] = -mu * (12*(x*x*y*y*y*y + x*x*x*x*y*y)) + 2*x*y*y;
  // res[1] = -mu * (-8./5. *x*y*y*y*(10*x*x + 3*y*y)) + 2*y*x*x;

  // res[0] = 4.0*pow(x,3)*(pow(x,4) + pow(y,4)) - 12.0*pow(x,3)*pow(y,4) - 12.0*x*x - 12.0*y*y;
  // res[1] = 4.0*pow(x,4)*pow(y,3) + 24.0*x*y + 16.0*pow(x,6)*y - 12.0*x*x*y*(pow(x,4) + pow(y,4));
  if(fabs(t) <= 1e-10) {
    res[0] = 2.0 * cos(y) * sin(x);
    res[1] = -2.0 * cos(x) * sin(y);
  } else {
    res[0] = 0.0;
    res[1] = 0.0;
  }
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
  feFunction *funSolP = new feFunction(fSolP, stokesParam);
  feFunction *funZero = new feFunction(fZero, {});
  feFunction *funOne = new feFunction(fOne, {});
  feVectorFunction *funSource = new feVectorFunction(fSource, stokesParam);

  int nIter = 2;
  std::vector<double> normL2_U(2 * nIter, 0.0);
  std::vector<double> normL2_V(2 * nIter, 0.0);
  std::vector<double> normL2_P(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);

  std::string meshName = "../../data/squareTaylorGreen2.msh";
  std::string root = "../../data/squareTaylorGreen";
  std::string metricMeshName, nextMeshName;

  feMesh2DP1 *mesh = new feMesh2DP1(meshName, false);
  feMesh2DP1 *otherMesh;
  feMesh2DP1::mapType physicalEntities = mesh->getPhysicalEntitiesDescription();

  feSpace *U_angle, *U_haut, *U_gauche, *U_surface;
  feSpace *V_angle, *V_haut, *V_gauche, *V_surface;
  feSpace *P_point, *P_surface;

  feSpace *U_angle2, *U_haut2, *U_gauche2, *U_surface2;
  feSpace *V_angle2, *V_haut2, *V_gauche2, *V_surface2;
  feSpace *P_point2, *P_surface2;

  feMetaNumber *metaNumber, *otherMetaNumber;
  feSolution *sol, *otherSol;

  // BDF2Solver *solver;

  std::string elemType = "P2P1";

  for(int iter = 0; iter < nIter; ++iter) {
    // std::string meshName = "../../data/squareTaylorGreen2.msh";
    std::string meshName = "../../data/Square/squareNS" + std::to_string(iter + 1) + ".msh";
    // std::string meshName = "../../data/squarePression" + std::to_string(iter+1) + ".msh";

    // mesh = new feMesh2DP1(meshName, false);
    nElm[iter] = mesh->getNbInteriorElems();

    if(elemType == "P2P1") {
      U_angle = new feSpace1DP2(mesh, "U", "Angle", funSolU);
      U_haut = new feSpace1DP2(mesh, "U", "Haut", funSolU);
      U_gauche = new feSpace1DP2(mesh, "U", "Gauche", funSolU);
      U_surface = new feSpaceTriP2(mesh, "U", "Surface", funSolU);

      V_angle = new feSpace1DP2(mesh, "V", "Angle", funSolV);
      V_haut = new feSpace1DP2(mesh, "V", "Haut", funSolV);
      V_gauche = new feSpace1DP2(mesh, "V", "Gauche", funSolV);
      V_surface = new feSpaceTriP2(mesh, "V", "Surface", funSolV);

      P_point = new feSpace1DP0(mesh, "P", "PointPression", funSolP);
      P_surface = new feSpaceTriP1(mesh, "P", "Surface", funZero);
    } else if(elemType == "P3P2") {
      U_angle = new feSpace1DP3(mesh, "U", "Angle", funSolU);
      U_haut = new feSpace1DP3(mesh, "U", "Haut", funSolU);
      U_gauche = new feSpace1DP3(mesh, "U", "Gauche", funSolU);
      U_surface = new feSpaceTriP3(mesh, "U", "Surface", funSolU);

      V_angle = new feSpace1DP3(mesh, "V", "Angle", funSolV);
      V_haut = new feSpace1DP3(mesh, "V", "Haut", funSolV);
      V_gauche = new feSpace1DP3(mesh, "V", "Gauche", funSolV);
      V_surface = new feSpaceTriP3(mesh, "V", "Surface", funSolV);

      P_point = new feSpace1DP0(mesh, "P", "PointPression", funSolP);
      P_surface = new feSpaceTriP2(mesh, "P", "Surface", funZero);
    } else {
      printf("Error - Unknown finite element type.\n");
    }

    std::vector<feSpace *> fespace = {U_angle, U_haut,   U_gauche,  U_surface, V_angle,
                                      V_haut,  V_gauche, V_surface, P_point,   P_surface};

    std::vector<feSpace *> feEssBC = {
      U_angle, U_haut, U_gauche, V_angle, V_haut, V_gauche, P_point,
      // P_surface
    };

    metaNumber = new feMetaNumber(mesh, fespace, feEssBC);

    sol = new feSolution(mesh, fespace, feEssBC, metaNumber);

    // Formes (bi)lineaires
    int dQuad = 15;
    std::vector<feSpace *> spacesNS2D = {U_surface, V_surface, P_surface};

    feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, dQuad,
                                              new feSysElm_2D_NavierStokes(stokesParam, funSource));
    // feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, dQuad, new
    // feSysElm_2D_NavierStokes(stokesParam, nullptr));

    std::vector<feBilinearForm *> formMatrices = {NS2D};
    std::vector<feBilinearForm *> formResiduals = {NS2D};

    feNorm *normU = new feNorm(U_surface, mesh, dQuad, funSolU);
    feNorm *normV = new feNorm(V_surface, mesh, dQuad, funSolV);
    feNorm *normP = new feNorm(P_surface, mesh, dQuad, funSolP);
    std::vector<feNorm *> norms = {normU, normV, normP};

    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);

    feTolerances tol{1e-9, 1e-8, 10};

    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 10 * pow(2, iter);

    // BDF2Solver *solver = new BDF2Solver(tol, metaNumber, linearSystem, sol, norms, mesh, t0, t1,
    // nTimeSteps);
    StationarySolver *solver =
      new StationarySolver(tol, metaNumber, linearSystem, sol, norms, mesh);
    solver->makeSteps(nTimeSteps, fespace);

    // std::vector<double> normL2BDF(3 * nTimeSteps, 0.0);
    // std::vector<double> normL2DC3(3 * nTimeSteps, 0.0);
    // solveBDF2(normL2BDF, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, norms,
    // mesh, fespace);
    // sol->initializeTemporalSolution(t0, t1, nTimeSteps);
    // solveBDF2(normL2BDF, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol,
    //          norms, mesh, fespace);

    // normL2_U[2*iter] = *std::max_element(solver->getNorm(0).begin(), solver->getNorm(0).end());
    // normL2_V[2*iter] = *std::max_element(solver->getNorm(1).begin(), solver->getNorm(1).end());
    // normL2_P[2*iter] = *std::max_element(solver->getNorm(2).begin(), solver->getNorm(2).end());
    // normL2_U[2*iter] = solver->getNorm(0)[nTimeSteps-1];
    // normL2_V[2*iter] = solver->getNorm(1)[nTimeSteps-1];
    // normL2_P[2*iter] = solver->getNorm(2)[nTimeSteps-1];
    normL2_U[2 * iter] = solver->getNorm(0);
    normL2_V[2 * iter] = solver->getNorm(1);
    normL2_P[2 * iter] = solver->getNorm(2);

    // normU->computeL2Norm(metaNumber, sol, mesh);
    // normV->computeL2Norm(metaNumber, sol, mesh);
    // normP->computeL2Norm(metaNumber, sol, mesh);
    // normL2_U[2 * iter] = normU->getNorm();
    // normL2_V[2 * iter] = normV->getNorm();
    // normL2_P[2 * iter] = normP->getNorm();

    // std::string vtkFile = "../../data/solutionManufacturee" + std::to_string(iter+1) + ".vtk";
    // std::string vtkFile = "../../data/taylorGreenUnsteady" + std::to_string(iter+1) + ".vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);

    metricMeshName = root + "sizeField" + std::to_string(iter + 1) + ".msh";
    nextMeshName = root + "adapted" + std::to_string(iter + 1) + ".msh";

    std::vector<double> estErreur(2, 0.);
    feRecovery *recU = new feRecovery(metaNumber, U_surface, mesh, sol, estErreur, funZero,
                                      meshName, metricMeshName);
    // feRecovery *recV = new feRecovery(metaNumber, V_surface, mesh, sol, estErreur, funZero,
    // meshName, metricMeshName); feRecovery *recP = new feRecovery(metaNumber, P_surface, mesh,
    // sol, estErreur, funZero);

    double modelSize = 2. * M_PI;

    feMetricOptions metricOptions;
    metricOptions.computationMethod = 3;
    metricOptions.polynomialDegree = 2;
    metricOptions.nTargetVertices = 2500;
    metricOptions.eTargetError = 5e-4;
    metricOptions.hMin = modelSize / 1000000.;
    metricOptions.hMax = modelSize / 6.;
    metricOptions.LpNorm = 2.0;
    metricOptions.nPhi = 801;

    std::vector<feRecovery *> rec = {recU};
    int useAnalytical = 0;
    feCurvedAdapt foo(mesh, rec, metricOptions, meshName, metricMeshName, nextMeshName,
                      useAnalytical);
    std::string cmdGMSH = "gmsh " + nextMeshName + " &";
    system(cmdGMSH.c_str());

    std::cout << "Reading " << nextMeshName << std::endl;
    otherMesh = new feMesh2DP1(nextMeshName, false, physicalEntities);

    if(elemType == "P2P1") {
      U_angle2 = new feSpace1DP2(otherMesh, "U", "Angle", funSolU);
      U_haut2 = new feSpace1DP2(otherMesh, "U", "Haut", funSolU);
      U_gauche2 = new feSpace1DP2(otherMesh, "U", "Gauche", funSolU);
      U_surface2 = new feSpaceTriP2(otherMesh, "U", "Surface", funSolU);

      V_angle2 = new feSpace1DP2(otherMesh, "V", "Angle", funSolV);
      V_haut2 = new feSpace1DP2(otherMesh, "V", "Haut", funSolV);
      V_gauche2 = new feSpace1DP2(otherMesh, "V", "Gauche", funSolV);
      V_surface2 = new feSpaceTriP2(otherMesh, "V", "Surface", funSolV);

      P_point2 = new feSpace1DP0(otherMesh, "P", "PointPression", funSolP);
      P_surface2 = new feSpaceTriP1(otherMesh, "P", "Surface", funZero);
    } else if(elemType == "P3P2") {
      U_angle2 = new feSpace1DP3(otherMesh, "U", "Angle", funSolU);
      U_haut2 = new feSpace1DP3(otherMesh, "U", "Haut", funSolU);
      U_gauche2 = new feSpace1DP3(otherMesh, "U", "Gauche", funSolU);
      U_surface2 = new feSpaceTriP3(otherMesh, "U", "Surface", funSolU);

      V_angle2 = new feSpace1DP3(otherMesh, "V", "Angle", funSolV);
      V_haut2 = new feSpace1DP3(otherMesh, "V", "Haut", funSolV);
      V_gauche2 = new feSpace1DP3(otherMesh, "V", "Gauche", funSolV);
      V_surface2 = new feSpaceTriP3(otherMesh, "V", "Surface", funSolV);

      P_point2 = new feSpace1DP0(otherMesh, "P", "PointPression", funSolP);
      P_surface2 = new feSpaceTriP2(otherMesh, "P", "Surface", funZero);
    } else {
      printf("Error - Unknown finite element type.\n");
    }

    std::vector<feSpace *> otherSpaces = {U_angle2, U_haut2,   U_gauche2,  U_surface2, V_angle2,
                                          V_haut2,  V_gauche2, V_surface2, P_point2,   P_surface2};

    std::vector<feSpace *> otherSpacesEss = {
      U_angle2, U_haut2, U_gauche2, V_angle2, V_haut2, V_gauche2, P_point2,
      // P_surface
    };

    otherMetaNumber = new feMetaNumber(otherMesh, otherSpaces, otherSpacesEss);
    otherSol = new feSolution(otherMesh, otherSpaces, otherSpacesEss, otherMetaNumber);

    mesh->transfer(otherMesh, metaNumber, otherMetaNumber, solver->getSolutionContainer(), fespace,
                   feEssBC, otherSpaces);

    otherSol->setSolFromContainer(solver->getSolutionContainer(), 0);

    mesh = otherMesh;
    meshName = nextMeshName;

    delete normU;
    delete normV;
    delete normP;
    delete linearSystem;
    delete NS2D;
  }
  delete sol;
  delete metaNumber;
  delete mesh;
  delete funSource;
  delete funSolU;
  delete funSolV;
  delete funSolP;
  delete funZero;
  delete funOne;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_U[2 * i + 1] =
      -log(normL2_U[2 * i] / normL2_U[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normL2_V[2 * i + 1] =
      -log(normL2_V[2 * i] / normL2_V[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normL2_P[2 * i + 1] =
      -log(normL2_P[2 * i] / normL2_P[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  }
  // for(int i = 1; i < nIter; ++i) {
  //   maxNormL2DC3_U[2 * i + 1] = log(maxNormL2DC3_U[2 * (i - 1)] / maxNormL2DC3_U[2 * i]) /
  //   log(2.); maxNormL2DC3_V[2 * i + 1] = log(maxNormL2DC3_V[2 * (i - 1)] / maxNormL2DC3_V[2 * i])
  //   / log(2.); maxNormL2DC3_P[2 * i + 1] = log(maxNormL2DC3_P[2 * (i - 1)] / maxNormL2DC3_P[2 *
  //   i]) / log(2.);
  // }
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s\n", "nElm", "||u-uh||", "tauxU",
         "||v-vh||", "tauxV", "||p-ph||", "tauxP");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n", nElm[i],
           normL2_U[2 * i], normL2_U[2 * i + 1], normL2_V[2 * i], normL2_V[2 * i + 1],
           normL2_P[2 * i], normL2_P[2 * i + 1]);

    // printf("%12s \t %12s \t %12s\t %12s \t %12s\t %12s \t %12s\n", "nElm", "||E-DC3_U||",
    //        "Taux DC3_U", "||E-DC3_V||", "Taux DC3_V", "||E-DC3_P||", "Taux DC3_P");
    // for(int i = 0; i < nIter; ++i)
    //   printf("%12d \t %12.6e \t %12.6e\t %12.6e \t %12.6e\t %12.6e \t %12.6e\n", nElm[i],
    //          maxNormL2DC3_U[2 * i], maxNormL2DC3_U[2 * i + 1], maxNormL2DC3_V[2 * i],
    //          maxNormL2DC3_V[2 * i + 1], maxNormL2DC3_P[2 * i], maxNormL2DC3_P[2 * i + 1]);
#ifdef HAVE_PETSC
  petscFinalize();
#endif
  return 0;
}
