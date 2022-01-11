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

#include "feCurvedAdapt.h"

double fSol(const double t, const std::vector<double> x, const std::vector<double> par) {
  return pow(x[0], 4) * pow(x[1], 4);
}

double fSource(const double t, const std::vector<double> x, const std::vector<double> par) {
  double kd = par[0];
  return kd * 12.0 * x[0] * x[0] * x[1] * x[1] * (x[0] * x[0] + x[1] * x[1]);
}

int main(int argc, char **argv) {
#ifdef HAVE_PETSC
  petscInitialize(argc, argv);
#endif

  double kd = 1.0;
  feFunction funSol(fSol, {kd});
  feFunction funSource(fSource, {kd});

  int nIter = 1;
  std::vector<double> normL2(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);

  for(int iter = 0; iter < nIter; ++iter) {
    std::string meshName = "../../data/square" + std::to_string(iter + 2) + "Msh2.msh";

    feMesh2DP1 mesh(meshName, false);
    nElm[iter] = mesh.getNbInteriorElems();

    feSpace1DP2 U_angle(&mesh, "U", "Angle", &funSol);
    feSpace1DP2 U_haut(&mesh, "U", "Haut", &funSol);
    feSpace1DP2 U_gauche(&mesh, "U", "Gauche", &funSol);
    feSpaceTriP2 U_surface(&mesh, "U", "Surface", &funSol);

    std::vector<feSpace *> fespace = {&U_angle, &U_surface, &U_haut, &U_gauche};
    std::vector<feSpace *> feEssBC = {&U_angle, &U_haut, &U_gauche};

    feMetaNumber metaNumber(&mesh, fespace, feEssBC);
    feSolution sol(&mesh, fespace, feEssBC, &metaNumber);

    int dQuad = 5;
    std::vector<feSpace *> spaceDiffusion2D_U = {&U_surface};
    std::vector<feSpace *> spaceSource2D_U = {&U_surface};

    feBilinearForm diffU(spaceDiffusion2D_U, &mesh, dQuad, new feSysElm_2D_Diffusion(kd, nullptr));
    feBilinearForm sourceU(spaceSource2D_U, &mesh, dQuad, new feSysElm_2D_Source(1.0, &funSource));

    std::vector<feBilinearForm *> formMatrices = {&diffU};
    std::vector<feBilinearForm *> formResiduals = {&diffU, &sourceU};

    feNorm norm(&U_surface, &mesh, dQuad, &funSol);
    std::vector<feNorm *> norms = {&norm};

    feLinearSystemPETSc linearSystem(argc, argv, formMatrices, formResiduals, &metaNumber, &mesh);
    linearSystem.initialize();
    feTolerances tol{1e-9, 1e-8, 20};
    solveStationary(&normL2[2 * iter], tol, &metaNumber, &linearSystem, formMatrices, formResiduals,
                    &sol, norms, &mesh);
    linearSystem.finalize();

    // Compute metrics and adapt mesh
    std::vector<double> estErreur(2, 0.);
    feRecovery recovery(&metaNumber, &U_surface, &mesh, &sol, estErreur, &funSol);
    std::vector<feRecovery *> rec = {&recovery};

    // feMetricOptions metricOptions;
    // metricOptions.computationMethod = 1;
    // metricOptions.nTargetVertices = 1000 * pow(2,iter);
    // metricOptions.hMin = 0.000001;
    // metricOptions.hMax = 100;
    // metricOptions.LpNorm = 2.0;
    // metricOptions.nPhi = 501;

    // feCurvedAdapt foo(&mesh, rec);

    // std::cout<<estErreur[0]<<" -- "<<estErreur[1]<<std::endl;
    // normL2_recons[2*iter] = estErreur[0];
    // normL2_recons_exacte[2*iter] = estErreur[1];

    // feMetric *metric = new feMetric(rec, metricOptions);
    // metric->writeSizeFieldGmsh(meshName, metricMeshName);

    // std::string cmdMMGS = "mmgs " + metricMeshName + " -o " + nextMeshName + " -hgrad 100";
    // std::string cmdGMSH = "gmsh " + nextMeshName;

    // system(cmdMMGS.c_str());
    // system(cmdGMSH.c_str());

    // meshName = nextMeshName;
  }

#ifdef HAVE_PETSC
  petscFinalize();
#endif

  return 0;
}
