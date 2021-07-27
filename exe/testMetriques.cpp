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
#include "feRecovery.h"
#include "feMetric.h"
#include "fullMatrix.h"

double fSol(const double t, const std::vector<double> x, const std::vector<double> par) {
  // return pow(x[0],6);
  // return pow(x[0],3);
  // return pow(x[0],2);
  // return pow(x[0],4);
  return pow(x[0], 4) * pow(x[1], 4);
}

double fSource(const double t, const std::vector<double> x, const std::vector<double> par) {
  double kd = par[0];
  // return kd*30.*pow(x[0],4);
  // return kd*6*x[0];
  // return kd*2;
  // return kd * 12.0 * x[0]*x[0];
  return kd * 12.0 * x[0] * x[0] * x[1] * x[1] * (x[0] * x[0] + x[1] * x[1]);
}

int main(int argc, char **argv) {
  double kd = 1.0;
  feFunction *funSol = new feFunction(fSol, {kd});
  feFunction *funSource = new feFunction(fSource, {kd});

  int nIter = 1;
  std::vector<double> normL2_exacte(2 * nIter, 0.0);
  std::vector<double> normL2_recons(2 * nIter, 0.0);
  std::vector<double> normL2_recons_exacte(2 * nIter, 0.0);
  std::vector<int> nNodes(nIter, 0);

  std::string root, meshName, metricMeshName, nextMeshName;
  feMesh2DP1::mapType physicalEntities;

  for(int iter = 0; iter < nIter; ++iter) {
    if(iter == 0) {
      meshName = "../../data/square" + std::to_string(iter + 1) + "Msh2.msh";
      // meshName = "../../data/squareCoarse1.msh";
      size_t lastindex = meshName.find_last_of(".");
      root = meshName.substr(0, lastindex);
    }
    metricMeshName = root + "SizeField" + std::to_string(iter + 1) + ".msh";
    nextMeshName = root + "Adapted" + std::to_string(iter + 1) + ".msh";

    feMesh2DP1 *mesh;
    if(iter == 0) {
      mesh = new feMesh2DP1(meshName, false);
      physicalEntities = mesh->getPhysicalEntitiesDescription();
    } else {
      mesh = new feMesh2DP1(meshName, false, physicalEntities);
    }
    nNodes[iter] = mesh->getNbNodes();

    feSpace1DP2 U_angle = feSpace1DP2(mesh, "U", "Angle", funSol);
    feSpace1DP2 U_haut = feSpace1DP2(mesh, "U", "Haut", funSol);
    feSpace1DP2 U_gauche = feSpace1DP2(mesh, "U", "Gauche", funSol);
    feSpaceTriP2 U_surface = feSpaceTriP2(mesh, "U", "Surface", funSol);

    std::vector<feSpace *> fespace = {&U_angle, &U_haut, &U_gauche, &U_surface};
    std::vector<feSpace *> feEssBC = {&U_angle, &U_haut, &U_gauche};

    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);

    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);

    // Formes (bi)lineaires
    int nQuad = 5; // TODO : change to deg
    std::vector<feSpace *> spaceDiffusion2D_U = {&U_surface};
    feBilinearForm *diffU =
      new feBilinearForm(spaceDiffusion2D_U, mesh, nQuad, new feSysElm_2D_Diffusion(kd, nullptr));
    std::vector<feSpace *> spaceSource2D_U = {&U_surface};
    feBilinearForm *sourceU =
      new feBilinearForm(spaceSource2D_U, mesh, nQuad, new feSysElm_2D_Source(1.0, funSource));

    std::vector<feBilinearForm *> formMatrices = {diffU};
    std::vector<feBilinearForm *> formResiduals = {diffU, sourceU};

    feNorm *norm = new feNorm(&U_surface, mesh, nQuad, funSol);
    std::vector<feNorm *> norms = {norm};
    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
    linearSystem->initialize();
    feTolerances tol{1e-9, 1e-8, 3};
    solveStationary(&normL2_exacte[2 * iter], tol, metaNumber, linearSystem, formMatrices,
                    formResiduals, sol, norms, mesh);
    linearSystem->finalize();

    // std::string vtkFile = "../../data/visu.vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);

    std::vector<double> estErreur(2, 0.);
    feRecovery *rec = new feRecovery(metaNumber, &U_surface, mesh, sol, estErreur, funSol);

    std::cout << estErreur[0] << " -- " << estErreur[1] << std::endl;
    normL2_recons[2 * iter] = estErreur[0];
    normL2_recons_exacte[2 * iter] = estErreur[1];

    feMetricOptions metricOptions;
    metricOptions.nTargetVertices = 1000 * pow(2, iter);
    metricOptions.hMin = 0.000001;
    metricOptions.hMax = 100;
    metricOptions.computationMethod = 1;
    metricOptions.LpNorm = 2.0;
    metricOptions.nPhi = 501;

    feMetric *metric = new feMetric(rec, metricOptions);
    metric->writeSizeFieldGmsh(meshName, metricMeshName);

    std::string cmdMMGS = "mmgs " + metricMeshName + " -o " + nextMeshName + " -hgrad 100";
    std::string cmdGMSH = "gmsh " + nextMeshName;

    system(cmdMMGS.c_str());
    system(cmdGMSH.c_str());

    meshName = nextMeshName;

    delete rec;
    delete norm;
    delete linearSystem;
    delete sourceU;
    delete diffU;
    delete sol;
    delete metaNumber;
    delete mesh;
  }
  delete funSource;
  delete funSol;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_exacte[2 * i + 1] = -log(normL2_exacte[2 * i] / normL2_exacte[2 * (i - 1)]) /
                               log(sqrt(nNodes[i]) / sqrt(nNodes[i - 1]));
    normL2_recons[2 * i + 1] = -log(normL2_recons[2 * i] / normL2_recons[2 * (i - 1)]) /
                               log(sqrt(nNodes[i]) / sqrt(nNodes[i - 1]));
    normL2_recons_exacte[2 * i + 1] =
      -log(normL2_recons_exacte[2 * i] / normL2_recons_exacte[2 * (i - 1)]) /
      log(sqrt(nNodes[i]) / sqrt(nNodes[i - 1]));
  }
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s\n", "nNodes", "||u-uh||", "taux",
         "||p-uh||", "taux", "||u-p||", "taux");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n", nNodes[i],
           normL2_exacte[2 * i], normL2_exacte[2 * i + 1], normL2_recons[2 * i],
           normL2_recons[2 * i + 1], normL2_recons_exacte[2 * i], normL2_recons_exacte[2 * i + 1]);

  return 0;
}
