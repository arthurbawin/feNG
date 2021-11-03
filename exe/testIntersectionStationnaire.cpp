/*  A test case where we solve a Poisson equation over a 2D domain with features.
    The goal is to create a metric field from the computed solution,
    then intersect it with the metric field computed by Gmsh and obtained from the
    features of the geometry.

    Author : Arthur Bawin
*/

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
#include "feLinearSystemMklPardiso.h"
#include "feExporter.h"
#include "feRecovery.h"
#include "feMetric.h"
#include "feCurvedAdapt.h"

static double a = 1e-4;
static double b = 0.3;
static double c = 3.0;

static double xmin = -27.;
static double xmax = 573.;
static double ymin = -20.;
static double ymax = 420.;

static double xc = (xmin + xmax) / 2;
static double yc = ymin;
static double r = (ymax - ymin) / 2;

static double x_0 = 0.5;
static double y_0 = 0.5;

/* The dimensions of the 2D domain are (-27.3,573) x (-20,420).
   The analytical solution is an ondulating smoothed Heaviside function :
   u(x,y) = (1 + tanh(a*((x-x0) - b*sin(c*pi*y))))/2     */
double fSol(const double t, const std::vector<double> x, const std::vector<double> par) {
  double xi = x[0];
  double yi = x[1];
  // return (1. + tanh(a*((yi-y_0) - b*sin(c*M_PI*xi))))/2.;
  return pow((xi - x_0), 4) * pow((yi - y_0), 4);
}

// Demi-cercle
double fSol2(const double t, const std::vector<double> x, const std::vector<double> par) {
  double xi = x[0];
  double yi = x[1];
  double arg = r * r - ((xi - xc) * (xi - xc) + (yi - yc) * (yi - yc));
  return (1. + tanh(a * arg)) / 2.;
}

double fSource(const double t, const std::vector<double> x, const std::vector<double> par) {
  double kd = par[0];
  double xi = x[0];
  double yi = x[1];
  double d = tanh(a * (b * sin(M_PI * c * xi) - (yi - y_0)));
  // return - kd * a*a * d*(d*d-1.) - (a*b*c*c*M_PI*M_PI*sin(M_PI*c*xi)*(d*d-1.))/2. -
  // a*a*b*b*c*c*M_PI*M_PI*d*cos(M_PI*c*xi)*cos(M_PI*c*xi)*(d*d-1.);
  return kd * 12. * (xi - x_0) * (xi - x_0) * (yi - y_0) * (yi - y_0) *
         ((xi - x_0) * (xi - x_0) + (yi - y_0) * (yi - y_0));
}

double fSource2(const double t, const std::vector<double> x, const std::vector<double> par) {
  double kd = par[0];
  double xi = x[0];
  double yi = x[1];
  double arg1 = -r * r + (xi - xc) * (xi - xc) + (yi - yc) * (yi - yc);
  double arg2 = tanh(a * arg1);
  double arg3 = pow(arg2, 2);
  double lap = 2. * a * (arg3 - 1.) - 4. * a * a * arg2 * (arg3 - 1) * (xi - xc) * (xi - xc) -
               4. * a * a * arg2 * (arg3 - 1) * (yi - yc) * (yi - yc);
  return kd * lap;
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 0.0;
}

int main(int argc, char **argv) {
#ifdef HAVE_PETSC
  petscInitialize(argc, argv);
#endif

  double kd = 1.0;
  feFunction *funSol = new feFunction(fSol2, {kd});
  feFunction *funSource = new feFunction(fSource2, {kd});
  feFunction *funZero = new feFunction(fZero, {});

  int nIter = 1;
  std::vector<double> normL2(2 * nIter, 0.0);
  std::vector<double> normL2_U(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);

  std::string meshName;
  feMesh2DP1::mapType physicalEntities;

  bool computeOnlyErrorOnNextMesh = true;

  std::string root = "../../data/intersectionStationnaire/demiCercle/";

  for(int iter = 0; iter < nIter; ++iter) {
    if(iter == 0) {
      // meshName = root + "veryCoarseMsh2.msh";
      // meshName = root + "veryCoarseWithEntities.msh";
      // meshName = root + "veryCoarseWithEntitiesWithoutPoints.msh";
      meshName = root + "veryCoarseWithEntitiesWithoutPointsRenumberedFine.msh";
      // meshName = root + "veryCoarseWithoutPointsMsh2.msh";
    }

    feMesh2DP1 *mesh;
    if(iter == 0) {
      mesh = new feMesh2DP1(meshName, false);
      physicalEntities = mesh->getPhysicalEntitiesDescription();
    } else
      mesh = new feMesh2DP1(meshName, false, physicalEntities);

    nElm[iter] = mesh->getNbInteriorElems();

    feSpace *U_bord, *U_surface;

    int deg = 1;
    switch(deg) {
    case 1:
      U_bord = new feSpace1DP1(mesh, "U", "Bord", funSol);
      U_surface = new feSpaceTriP1(mesh, "U", "Domaine", funSol);
      break;
    case 2:
      U_bord = new feSpace1DP2(mesh, "U", "Bord", funSol);
      U_surface = new feSpaceTriP2(mesh, "U", "Domaine", funSol);
      break;
    case 3:
      U_bord = new feSpace1DP3(mesh, "U", "Bord", funSol);
      U_surface = new feSpaceTriP3(mesh, "U", "Domaine", funSol);
      break;
    case 4:
      U_bord = new feSpace1DP4(mesh, "U", "Bord", funSol);
      U_surface = new feSpaceTriP4(mesh, "U", "Domaine", funSol);
      break;
    default:
      U_bord = new feSpace1DP1(mesh, "U", "Bord", funSol);
      U_surface = new feSpaceTriP1(mesh, "U", "Domaine", funSol);
    }

    std::vector<feSpace *> fespace = {U_bord, U_surface};
    std::vector<feSpace *> feEssBC = {U_bord};

    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);

    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);

    int dQuad = 15;
    std::vector<feSpace *> spaceDiffusion2D = {U_surface};
    std::vector<feSpace *> spaceSource2D = {U_surface};

    feBilinearForm *diffU =
      new feBilinearForm(spaceDiffusion2D, mesh, dQuad, new feSysElm_2D_Diffusion(kd, nullptr));
    feBilinearForm *sourceU =
      new feBilinearForm(spaceSource2D, mesh, dQuad, new feSysElm_2D_Source(1.0, funSource));

    std::vector<feBilinearForm *> formMatrices = {diffU};
    std::vector<feBilinearForm *> formResiduals = {diffU, sourceU};

    feNorm *normU = new feNorm(U_surface, mesh, dQuad, funSol);
    std::vector<feNorm *> norms = {normU};

    // feLinearSystemPETSc *linearSystem = new feLinearSystemPETSc(argc, argv, formMatrices,
    // formResiduals, metaNumber, mesh);
    feLinearSystemMklPardiso *linearSystem =
      new feLinearSystemMklPardiso(formMatrices, formResiduals, metaNumber, mesh);

    feTolerances tol{1e-9, 1e-8, 5};

    StationarySolver solver(tol, metaNumber, linearSystem, sol, norms, mesh);
    solver.makeSteps(0, fespace);
    normL2[2 * iter] = solver.getNorm(0);
    normL2_U[2 * iter] = normU->getNorm();

    std::string vtkFile = root + "sol.vtk";
    feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);

    normU->computeL2Norm(metaNumber, sol, mesh);
    double errorL2 = normU->getNorm();
    printf("========== ITER %d ========== Computed L2 norm on mesh %s = %10.10e\n", iter,
           meshName.c_str(), errorL2);

    // if(iter < 3 || !computeOnlyErrorOnNextMesh){

    std::string metricMeshName = root + "sizeField" + std::to_string(iter + 1) + ".msh";
    std::string nextMeshName = root + "adapted" + std::to_string(iter + 1) + ".msh";

    std::vector<double> estErreur(2, 0.);
    feRecovery *recU = new feRecovery(metaNumber, U_surface, mesh, sol, estErreur, funZero,
                                      meshName, metricMeshName);

    feMetricOptions metricOptions;
    metricOptions.computationMethod = 3;
    metricOptions.polynomialDegree = deg;
    metricOptions.eTargetError = 1e-3;
    metricOptions.nTargetVertices = 1000;

    double modelSize = fmax(ymax - ymin, xmax - xmin);
    metricOptions.hMin = modelSize / 10000;
    metricOptions.hMax = modelSize / 10;
    metricOptions.LpNorm = 2.0;
    metricOptions.nPhi = 501;

    switch(metricOptions.computationMethod) {
    case 1: {
      feMetric *metric = new feMetric(recU, metricOptions);
      metric->writeSizeFieldGmsh(meshName, metricMeshName);
      metric->writeSizeFieldSol2D("toIntersect2D_sol2D.sol");
      metric->writeSizeFieldSol3D("toIntersect2D.sol");
      metric->drawEllipsoids("ellipsesFromSimulation.pos");

      std::string cmdMMGS = "mmg2d " + metricMeshName + " -o " + nextMeshName + " -hgrad 3";
      std::string cmdGMSH = "gmsh " + nextMeshName + " &";
      system(cmdMMGS.c_str());
      system(cmdGMSH.c_str());

      delete metric;
    } break;

    case 2: break;

    case 3: {
      std::vector<feRecovery *> rec = {recU};
      int useAnalytical = 0;
      feCurvedAdapt foo(mesh, rec, metricOptions, meshName, metricMeshName, nextMeshName,
                        useAnalytical, funSol);
      // system("gmsh output.msh &");
      // nextMeshName = "output.msh";
    } break;

    default:;
    }

    meshName = nextMeshName;

    delete recU;
    delete normU;
    delete linearSystem;
    // }
    delete sourceU;
    delete diffU;
    delete sol;
    delete metaNumber;
    delete U_bord;
    delete U_surface;
    delete mesh;
  }
  delete funSource;
  delete funSol;

  // // Calcul du taux de convergence
  // for(int i = 1; i < nIter; ++i) {
  //   normL2_U[2 * i + 1] =
  //     -log(normL2_U[2 * i] / normL2_U[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  //   normL2_V[2 * i + 1] =
  //     -log(normL2_V[2 * i] / normL2_V[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  // }
  // printf("%12s \t %12s \t %12s \t %12s \t %12s\n", "nElm", "||E_U||", "p_U", "||E_V||", "p_V");
  // for(int i = 0; i < nIter; ++i)
  //   printf("%12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n", nElm[i], normL2_U[2 * i],
  //          normL2_U[2 * i + 1], normL2_V[2 * i], normL2_V[2 * i + 1]);

#ifdef HAVE_PETSC
  petscFinalize();
#endif

  return 0;
}