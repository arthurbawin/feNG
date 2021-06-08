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

double fInlet(const double t, const std::vector<double> &pos, const std::vector<double> &par){
  double y = pos[1];
  double h = 2.0;
  if(y <= h/2.0)
    return 1.0;
  else
    return 2.0;
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par){
  return 0.0;
}


int main(int argc, char** argv){

  // Diamètre du cylindre et nombre de Reynolds
  double r = 0.25;
  double d = 2*r;
  double U = 2.0;
  double rho = 1.0;
  double Re = 100.0;
  double mu = rho*U*d/Re;

  std::vector<double> stokesParam = {rho, mu};
  feFunction *funZero  = new feFunction(fZero, {});
  feFunction *funInlet = new feFunction(fInlet,{});

  int nIter = 5;

  std::string root, meshName, metricMeshName, nextMeshName;
  feMesh2DP1::mapType physicalEntities;

  for(int iter = 0; iter < nIter; ++iter){
    if(iter == 0){
      meshName = "../../data/cyl" + std::to_string(iter+1) + "Msh2.msh";
      size_t lastindex = meshName.find_last_of(".");
      root = meshName.substr(0, lastindex);
    }
    metricMeshName = root + "SizeField" + std::to_string(iter+1) + ".msh";
    nextMeshName = root + "Adapted" + std::to_string(iter+1) + ".msh";

    feMesh2DP1 *mesh;
    if(iter == 0){
      mesh = new feMesh2DP1(meshName, false);
      physicalEntities = mesh->getPhysicalEntitiesDescription();
    } else{
      mesh = new feMesh2DP1(meshName, false, physicalEntities);
    }

    feSpace1DP2  U_bas     (mesh, "U", "Bas",      funZero);
    feSpace1DP2  U_haut    (mesh, "U", "Haut",     funZero);
    feSpace1DP2  U_gauche  (mesh, "U", "Gauche",   funInlet); // Inlet
    feSpace1DP2  U_cylindre(mesh, "U", "Cylindre", funZero);
    feSpaceTriP2 U_surface (mesh, "U", "Surface",  funZero);

    feSpace1DP2  V_bas      = feSpace1DP2(mesh, "V", "Bas",      funZero);
    feSpace1DP2  V_haut     = feSpace1DP2(mesh, "V", "Haut",     funZero);
    feSpace1DP2  V_gauche   = feSpace1DP2(mesh, "V", "Gauche",   funZero);
    feSpace1DP2  V_cylindre = feSpace1DP2(mesh, "V", "Cylindre", funZero);
    feSpaceTriP2 V_surface  = feSpaceTriP2(mesh,"V", "Surface",  funZero);

    // feSpace1DP0  P_point   = feSpace1DP0(mesh, "P", "PointPression", funSolP);
    feSpace1DP1  P_droite   = feSpace1DP1(mesh, "P", "Droite",   funZero);
    feSpaceTriP1 P_surface  = feSpaceTriP1(mesh,"P", "Surface",  funZero);

    std::vector<feSpace*> fespace = {&U_bas, &U_haut, &U_gauche, &U_cylindre, &U_surface,
                                     &V_bas, &V_haut, &V_gauche, &V_cylindre, &V_surface, 
                                     &P_droite, &P_surface};

    std::vector<feSpace*> feEssBC = {&U_gauche, &U_cylindre,
                                     &V_bas, &V_haut, &V_gauche, &V_cylindre, 
                                     &P_droite};

    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // metaNumber->printNumberings();

    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);

    // Formes (bi)lineaires
    int nQuad = 16;
    std::vector<feSpace*> spacesNS2D = {&U_surface, &V_surface, &P_surface};
    feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, nQuad, new feSysElm_2D_NavierStokes(stokesParam, nullptr));

    std::vector<feBilinearForm*> formMatrices  = {NS2D};
    std::vector<feBilinearForm*> formResiduals = {NS2D};

    std::vector<feNorm*> norms = {};
    feLinearSystemPETSc *linearSystem = new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);    
    linearSystem->initialize();
    feTolerances tol{1e-9, 1e-8, 20};
    double fooNorm;
    solveStationary(&fooNorm, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, norms, mesh);
    
    // double t0 = 0.;
    // double t1 = 10.0;
    // int nTimeSteps = 100 * pow(2, iter);
    // sol->initializeTemporalSolution(t0,t1,nTimeSteps);
    // std::vector<double> normL2BDF(3*nTimeSteps,0.0);
    // std::vector<double> normL2DC3(3*nTimeSteps,0.0);
    // // solveBDF2(normL2BDF, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, norms, mesh, fespace);
    // solveDC3(normL2BDF, normL2DC3, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, norms, mesh, fespace);
    // // normL2_U[2*iter] = *std::max_element(normL2BDF.begin(), normL2BDF.end());
    // // maxNormL2BDF[2*iter] = *std::max_element(normL2BDF.begin(), normL2BDF.end());
    // // for (int i = 0; i < nTimeSteps; ++i){
    // //   maxNormL2DC3_U[2*iter] = fmax(maxNormL2DC3_U[2*iter], normL2DC3[3*i]  );
    // //   maxNormL2DC3_V[2*iter] = fmax(maxNormL2DC3_V[2*iter], normL2DC3[3*i+1]);
    // //   maxNormL2DC3_P[2*iter] = fmax(maxNormL2DC3_P[2*iter], normL2DC3[3*i+2]);
    // // }
    
    linearSystem->finalize();

    std::vector<double> estErreur(2, 0.);
    feRecovery *rec = new feRecovery(metaNumber, &U_surface, mesh, sol, nQuad, estErreur, funZero);

    std::cout<<estErreur[0]<<" -- "<<estErreur[1]<<std::endl;
    // normL2_recons[2*iter] = estErreur[0];
    // normL2_recons_exacte[2*iter] = estErreur[1];

    feMetricOptions metricOptions;
    metricOptions.nTargetVertices = 2000;
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

    // std::string vtkFile = "../../data/solutionManufacturee" + std::to_string(iter+1) + ".vtk";
    // std::string vtkFile = "../../data/cylindreNonConfine.vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);

    // normU->computeL2Norm(metaNumber, sol, mesh);
    // normV->computeL2Norm(metaNumber, sol, mesh);
    // normP->computeL2Norm(metaNumber, sol, mesh);
    // normL2_U[2*iter] = normU->getNorm();
    // normL2_V[2*iter] = normV->getNorm();
    // normL2_P[2*iter] = normP->getNorm();

    // delete normU;
    // delete normV;
    // delete normP;
    delete linearSystem;
    delete NS2D;
    delete sol;
    delete metaNumber;
    delete mesh;
  }
  delete funZero;
  delete funInlet;

  // // Calcul du taux de convergence
  // for(int i = 1; i < nIter; ++i){
  //   normL2_U[2*i+1] = - log(normL2_U[2*i]/normL2_U[2*(i-1)])/log(sqrt(nElm[i])/sqrt(nElm[i-1]));
  //   normL2_V[2*i+1] = - log(normL2_V[2*i]/normL2_V[2*(i-1)])/log(sqrt(nElm[i])/sqrt(nElm[i-1]));
  //   normL2_P[2*i+1] = - log(normL2_P[2*i]/normL2_P[2*(i-1)])/log(sqrt(nElm[i])/sqrt(nElm[i-1]));
  // }
  // for(int i = 1; i < nIter; ++i){
  //   maxNormL2DC3_U[2*i+1] = log(maxNormL2DC3_U[2*(i-1)]/maxNormL2DC3_U[2*i])/log(2.);
  //   maxNormL2DC3_V[2*i+1] = log(maxNormL2DC3_V[2*(i-1)]/maxNormL2DC3_V[2*i])/log(2.);
  //   maxNormL2DC3_P[2*i+1] = log(maxNormL2DC3_P[2*(i-1)]/maxNormL2DC3_P[2*i])/log(2.);
  // }
  // printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s\n", "nElm", "||u-uh||", "tauxU", "||v-vh||", "tauxV", "||p-ph||", "tauxP");
  // for(int i = 0; i < nIter; ++i)
  // printf("%12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n", nElm[i], 
  //   normL2_U[2*i], normL2_U[2*i+1],
  //   normL2_V[2*i], normL2_V[2*i+1],
  //   normL2_P[2*i], normL2_P[2*i+1]);

  //  printf("%12s \t %12s \t %12s\t %12s \t %12s\t %12s \t %12s\n", "nElm", "||E-DC3_U||", "Taux DC3_U", "||E-DC3_V||", "Taux DC3_V", "||E-DC3_P||", "Taux DC3_P");
  // for(int i = 0; i < nIter; ++i)
  //   printf("%12d \t %12.6e \t %12.6e\t %12.6e \t %12.6e\t %12.6e \t %12.6e\n", 
  //     nElm[i], maxNormL2DC3_U[2*i], maxNormL2DC3_U[2*i+1], maxNormL2DC3_V[2*i], maxNormL2DC3_V[2*i+1], maxNormL2DC3_P[2*i], maxNormL2DC3_P[2*i+1]);


  return 0;
}
