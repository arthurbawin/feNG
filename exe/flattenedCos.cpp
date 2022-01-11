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

static double x_0 = 0.5;
static double y_0 = 0.5;

static double b = 15.;

double fSolU(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  // return sqrt((1+b*b)/(1+b*b*sin(x)*cos(y)*sin(x)*cos(y))) * sin(x)*cos(y); // flattenedCos
  return pow((x - x_0), 4) * pow((y - y_0), 4); // x4y4
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 0.0;
}

double fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double kd = par[0];
  double x = pos[0];
  double y = pos[1];
  double cx = cos(x);
  double cy = cos(y);
  double sx = sin(x);
  double sy = sin(y);

  double f1 = b * b * cy * cy * sx * sx + 1.;

  double f2 = sqrt(1. + b * b);

  double lap =
    2 * b * b * cy * cy * cy * sx * sx * sx * f2 / pow(f1, 3. / 2.) - 2. * cy * sx * f2 / sqrt(f1) -
    3. * b * b * cx * cx * cy * cy * cy * sx * f2 / pow(f1, 3. / 2.) -
    3 * b * b * cy * sx * sx * sx * sy * sy * f2 / pow(f1, 3. / 2.) +
    3 * b * b * b * b * cx * cx * cy * cy * cy * cy * cy * sx * sx * sx * f2 / pow(f1, 5. / 2.) +
    3 * b * b * b * b * cy * cy * cy * sx * sx * sx * sx * sx * sy * sy * f2 / pow(f1, 5. / 2.);

  // return kd * lap; // flattenedCos
  return kd * 12. * (x - x_0) * (x - x_0) * (y - y_0) * (y - y_0) *
         ((x - x_0) * (x - x_0) + (y - y_0) * (y - y_0)); // x4y4
}

int main(int argc, char **argv) {
#ifdef HAVE_PETSC
  petscInitialize(argc, argv);
#endif

  double kd = 1.0;

  feFunction *funSolU = new feFunction(fSolU, {});
  feFunction *funZero = new feFunction(fZero, {});
  feFunction *funSource = new feFunction(fSource, {kd});

  int nIter = 5;
  std::vector<double> normL2_U(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);

  std::string dir = "../../data/convergenceFlattenedCos/x4y4/";
  std::string meshName = dir + "mesh.msh";
  std::string metricMeshName, nextMeshName;
  std::string infoFileName = dir + "initialL2Error.txt";

  feMesh2DP1 *mesh = new feMesh2DP1(meshName, false);
  feMesh2DP1 *otherMesh;
  feMesh2DP1::mapType physicalEntities = mesh->getPhysicalEntitiesDescription();

  feSpace *U_angle, *U_haut, *U_gauche, *U_surface;
  feSpace *U_angle2, *U_haut2, *U_gauche2, *U_surface2;

  feMetaNumber *metaNumber, *otherMetaNumber;
  feSolution *sol, *otherSol;

  for(int iter = 0; iter < nIter; ++iter) {
    nElm[iter] = mesh->getNbInteriorElems();

    int deg = 1;
    U_angle = new feSpace1DP1(mesh, "U", "Angle", funSolU);
    U_haut = new feSpace1DP1(mesh, "U", "Haut", funSolU);
    U_gauche = new feSpace1DP1(mesh, "U", "Gauche", funSolU);
    U_surface = new feSpaceTriP1(mesh, "U", "Surface", funSolU);

    std::vector<feSpace *> fespace = {U_angle, U_haut, U_gauche, U_surface};

    std::vector<feSpace *> feEssBC = {U_angle, U_haut, U_gauche};

    metaNumber = new feMetaNumber(mesh, fespace, feEssBC);

    sol = new feSolution(mesh, fespace, feEssBC, metaNumber);

    // Formes (bi)lineaires
    int dQuad = 15;
    std::vector<feSpace *> spaceDiffusion2D = {U_surface};
    std::vector<feSpace *> spaceSource2D = {U_surface};

    feBilinearForm *diffU =
      new feBilinearForm(spaceDiffusion2D, mesh, dQuad, new feSysElm_2D_Diffusion(kd, nullptr));
    feBilinearForm *sourceU =
      new feBilinearForm(spaceSource2D, mesh, dQuad, new feSysElm_2D_Source(1.0, funSource));

    std::vector<feBilinearForm *> formMatrices = {diffU};
    std::vector<feBilinearForm *> formResiduals = {diffU, sourceU};

    feNorm *normU = new feNorm(U_surface, mesh, dQuad, funSolU);
    std::vector<feNorm *> norms = {normU};

    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);

    feTolerances tol{1e-9, 1e-8, 10};

    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 10 * pow(2, iter);

    StationarySolver *solver =
      new StationarySolver(tol, metaNumber, linearSystem, sol, norms, mesh);
    solver->makeSteps(nTimeSteps, fespace);

    normL2_U[2 * iter] = solver->getNorm(0);
    std::cout << "Computed error = " << normL2_U[2 * iter] << std::endl;

    FILE *infoFile = fopen(infoFileName.c_str(), "a");
    fprintf(infoFile, "erreur de la norme = %12.5e\n", normL2_U[2 * iter]);

    std::string vtkFile = "../../data/flattenedCos" + std::to_string(iter + 1) + ".vtk";
    feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);

    double modelSize = 1.0;

    feMetricOptions metricOptions;
    metricOptions.computationMethod = 3;
    metricOptions.polynomialDegree = 1;
    metricOptions.nTargetVertices = 2500;
    metricOptions.hMin = modelSize / 1000000.;
    metricOptions.hMax = modelSize;
    metricOptions.LpNorm = 2.0;
    metricOptions.nPhi = 801;

    int eCibleExp = 6;
    double eCible = pow(10, -eCibleExp);
    metricOptions.eTargetError = eCible;

    int useAnalytical = 0;

    metricMeshName = dir + "sizeFieldJF" + std::to_string(iter + 1) + ".msh";
    nextMeshName = dir + "adaptedJF" + std::to_string(iter + 1) + ".msh";

    if(metricOptions.computationMethod == 1) {
      metricMeshName = dir + "Kuate_P" + std::to_string(deg) + "_" +
                       std::to_string(metricOptions.nTargetVertices) + "vert_iter" +
                       std::to_string(iter + 1) + "SizeField.msh";
      nextMeshName = dir + "Kuate_P" + std::to_string(deg) + "_" +
                     std::to_string(metricOptions.nTargetVertices) + "vert_iter" +
                     std::to_string(iter + 1) + "Adapted.msh";
      infoFileName = dir + "errorL2_Kuate_P" + std::to_string(deg) + "_" +
                     std::to_string(metricOptions.nTargetVertices) + "vert_iter" +
                     std::to_string(iter + 1) + ".txt";
    } else {
      if(useAnalytical > 0) {
        metricMeshName = dir + "JF_P" + std::to_string(deg) + "_1em" + std::to_string(eCibleExp) +
                         "_iter" + std::to_string(iter + 1) + "SizeField.msh";
        nextMeshName = dir + "JF_P" + std::to_string(deg) + "_1em" + std::to_string(eCibleExp) +
                       "_iter" + std::to_string(iter + 1) + "Adapted.msh";
        infoFileName = dir + "errorL2_JF_P" + std::to_string(deg) + "_1em" +
                       std::to_string(eCibleExp) + "_iter" + std::to_string(iter + 1) + ".txt";
      } else {
        metricMeshName = dir + "recovery_JF_P" + std::to_string(deg) + "_1em" +
                         std::to_string(eCibleExp) + "_iter" + std::to_string(iter + 1) +
                         "SizeField.msh";
        nextMeshName = dir + "recovery_JF_P" + std::to_string(deg) + "_1em" +
                       std::to_string(eCibleExp) + "_iter" + std::to_string(iter + 1) +
                       "Adapted.msh";
        infoFileName = dir + "recovery_errorL2_JF_P" + std::to_string(deg) + "_1em" +
                       std::to_string(eCibleExp) + "_iter" + std::to_string(iter + 1) + ".txt";
      }
    }

    std::vector<double> estErreur(2, 0.);
    feRecovery *recU = new feRecovery(metaNumber, U_surface, mesh, sol, estErreur, funZero,
                                      meshName, metricMeshName);

    // feMetric *metric = new feMetric(recU, metricOptions);
    // metric->writeSizeFieldGmsh(meshName, metricMeshName);
    // metric->writeSizeFieldSol2D("toIntersect2D_sol2D.sol");
    // metric->writeSizeFieldSol3D("toIntersect2D.sol");
    // metric->drawEllipsoids("ellipsesFromSimulation.pos");

    // std::string cmdMMGS = "mmg2d " + metricMeshName + " -o " + nextMeshName + " -hgrad 100";
    // std::string cmdGMSH = "gmsh " + nextMeshName + " &";
    // system(cmdMMGS.c_str());
    // system(cmdGMSH.c_str());

    // delete metric;

    std::vector<feRecovery *> rec = {recU};

    feCurvedAdapt foo(mesh, rec, metricOptions, metricMeshName, nextMeshName, useAnalytical,
                      funSolU);
    std::string cmdGMSH = "gmsh " + nextMeshName + " &";
    system(cmdGMSH.c_str());

    std::cout << "Reading " << nextMeshName << std::endl;
    otherMesh = new feMesh2DP1(nextMeshName, false, physicalEntities);

    U_angle2 = new feSpace1DP1(otherMesh, "U", "Angle", funSolU);
    U_haut2 = new feSpace1DP1(otherMesh, "U", "Haut", funSolU);
    U_gauche2 = new feSpace1DP1(otherMesh, "U", "Gauche", funSolU);
    U_surface2 = new feSpaceTriP1(otherMesh, "U", "Surface", funSolU);

    std::vector<feSpace *> otherSpaces = {
      U_angle2,
      U_haut2,
      U_gauche2,
      U_surface2,
    };

    std::vector<feSpace *> otherSpacesEss = {U_angle2, U_haut2, U_gauche2};

    otherMetaNumber = new feMetaNumber(otherMesh, otherSpaces, otherSpacesEss);
    otherSol = new feSolution(otherMesh, otherSpaces, otherSpacesEss, otherMetaNumber);

    mesh->transfer(otherMesh, metaNumber, otherMetaNumber, solver->getSolutionContainer(), fespace,
                   feEssBC, otherSpaces);

    otherSol->setSolFromContainer(solver->getSolutionContainer(), 0);

    mesh = otherMesh;
    meshName = nextMeshName;

    delete normU;
    delete linearSystem;
  }
  delete sol;
  delete metaNumber;
  delete mesh;
  delete funSource;
  delete funSolU;
  delete funZero;

  // // Calcul du taux de convergence
  // for(int i = 1; i < nIter; ++i) {
  //   normL2_U[2 * i + 1] = -log(normL2_U[2 * i] / normL2_U[2 * (i - 1)]) / log(sqrt(nElm[i]) /
  //   sqrt(nElm[i - 1]));
  // }
  // printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s\n", "nElm", "||u-uh||", "tauxU",
  //        "||v-vh||", "tauxV", "||p-ph||", "tauxP");
  // for(int i = 0; i < nIter; ++i)
  //   printf("%12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n", nElm[i],
  //          normL2_U[2 * i], normL2_U[2 * i + 1], normL2_V[2 * i], normL2_V[2 * i + 1],
  //          normL2_P[2 * i], normL2_P[2 * i + 1]);
#ifdef HAVE_PETSC
  petscFinalize();
#endif
  return 0;
}
