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

static double h = 16.0;
static double cx = h ;
static double cy = h / 2.;
static double r = 0.5;
static double d = 2 * r;
static double U = 1;
static double rho = 1.;
static double Re = 200;

double fSolU_inlet(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double y = pos[1];
  double coeff = 1.;
  // if (t<2.) coeff = pow(t,2)/4.;
  // else coeff=1.;
  return ((0.125*y - 0.5) + std::min(t/10,1.)* (-0.125*y + 0.5) +1) * coeff ;
}

// double fSolU_inlet(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
//   double y = pos[1];
//   double tEnd = 60.0;
//   double tRamp = fmin(1.0, t / tEnd);

//   return 1.0 + 2.0 * y / h * (1.0 - tRamp);
// }

// double fSolU_inlet_dot(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
//   double y = pos[1];
//   double tEnd = 10.0;
//   double tRamp = fmin(1.0, t / tEnd);

//   return  2.0 / h * (1.0 - tRamp);
// }

double fSolU_inlet_dot(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double y = pos[1];
  if(t< 10.) (-0.125*y + 0.5) / 10. ;
  else 0. ; 
}


double fSolP(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  return 0;
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 0.0;
}

void fgrad(const double t, const std::vector<double> &pos, const std::vector<double> &par,
             std::vector<double> &res) {

  res[0] = 0;
  res[1] = 0;

}

int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double mu = rho * U * d / Re;
  std::vector<double> stokesParam = {rho, mu};

  feFunction *funSolU_inlet = new feFunction(fSolU_inlet, stokesParam);
  feFunction *funSolU_inlet_dot = new feFunction(fSolU_inlet_dot, stokesParam);
  feFunction *funSolP = new feFunction(fSolP, stokesParam);
  feFunction *funZero = new feFunction(fZero, {});
  feVectorFunction *funGrad = new feVectorFunction(fgrad, stokesParam);


  int nIter = 1;
  std::vector<double> intL_u(2 * nIter, 0.0);
  std::vector<double> intL_v(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);

  for(int iter = 0; iter < nIter; ++iter) {
    std::string meshName = "../../data/VonKarman/vonKarman" + std::to_string(iter + 1) + "_V3.msh";
    // std::string meshName = "../../data/VonKarman/TestvonKarman" + std::to_string(iter + 1) + ".msh";

    feMesh2DP1 *mesh = new feMesh2DP1(meshName, false);
    nElm[iter] = mesh->getNbInteriorElems();

    feSpace1DP1 U_walls =  feSpace1DP1(mesh, "U", "walls", funZero);
    feSpace1DP1 U_inlet =  feSpace1DP1(mesh, "U", "inlet", funSolU_inlet); // Inlet
    feSpace1DP1 U_cylindre =  feSpace1DP1(mesh, "U", "cylinder", funZero);
    feSpaceTriP1 U_surface = feSpaceTriP1(mesh, "U", "surface", funZero);

    feSpace1DP1 V_walls =  feSpace1DP1(mesh, "V", "walls", funZero);
    feSpace1DP1 V_inlet =  feSpace1DP1(mesh, "V", "inlet", funZero); // Inlet
    feSpace1DP1 V_cylindre =  feSpace1DP1(mesh, "V", "cylinder", funZero);
    feSpaceTriP1 V_surface =  feSpaceTriP1(mesh, "V", "surface", funZero);

    feSpaceTriP1 P_surface =  feSpaceTriP1(mesh, "P", "surface", funSolP);

    feSpace1DP1 Lu_cylinder = feSpace1DP1(mesh, "Lu", "cylinder", funZero);
    feSpace1DP1 Lv_cylinder = feSpace1DP1(mesh, "Lv", "cylinder", funZero);
  
    feSpace1DP1 Bu_cylinder = feSpace1DP1(mesh, "Bu", "cylinder", funZero);
    feSpace1DP1 Bv_cylinder = feSpace1DP1(mesh, "Bv", "cylinder", funZero);


    std::vector<feSpace *> fespace = { &U_walls, &U_inlet, &U_cylindre, &U_surface, &V_walls, &V_inlet,  &V_cylindre, &V_surface, &P_surface, 
                                       &Lu_cylinder ,&Lv_cylinder , &Bu_cylinder, &Bv_cylinder};
 
    std::vector<feSpace *> feEssBC = { &U_inlet, &V_walls, &V_inlet};

    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
   
    double t0 = 0.;
    double t1 = 20.;
    int nTimeSteps =  200* pow(2, iter);
    TT[iter] = nTimeSteps;
    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);

    // Formes (bi)lineaires
    int nQuad = 10;
    std::vector<feSpace *> spacesNS2D = {&U_surface, &V_surface, &P_surface};
    std::vector<feSpace *> spaceWeak_U = {&U_cylindre ,&Lu_cylinder , &Bu_cylinder};
    std::vector<feSpace *> spaceWeak_V = {&V_cylindre ,&Lv_cylinder , &Bv_cylinder};

    // feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, nQuad, new
    // feSysElm_2D_NavierStokes(stokesParam, funSource));
    feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, nQuad,
                                              new feSysElm_2D_NavierStokes(stokesParam, nullptr));
    feBilinearForm *weakBC_UCylinder = new feBilinearForm(spaceWeak_U, mesh, nQuad,
                                                     new feSysElm_1D_weakBC_edo1(1.0, funZero));
    feBilinearForm *weakBC_VCylinder = new feBilinearForm(spaceWeak_V, mesh, nQuad,
                                                     new feSysElm_1D_weakBC_edo1(1.0, funZero));

    std::vector<feBilinearForm *> formMatrices = {NS2D , weakBC_UCylinder, weakBC_VCylinder};
    std::vector<feBilinearForm *> formResiduals ={NS2D , weakBC_UCylinder, weakBC_VCylinder};

    feNorm *IntLambda_U = new feNorm(&Lu_cylinder, mesh, nQuad, funSolP, funGrad ,"IntLambda");
    feNorm *IntLambda_V = new feNorm(&Lv_cylinder, mesh, nQuad, funSolP, funGrad , "IntLambda");
    //Group by type of fields, not type of solver 
    std::vector<feNorm *> norms = {IntLambda_U, IntLambda_V};
    // Systeme lineaire

    // feLinearSystemPETSc *linearSystem =
    // new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_MKL
    feLinearSystemMklPardiso *linearSystem =
      new feLinearSystemMklPardiso(formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
    // Resolution
    feTolerances tol{1e-5, 1e-5, 20};
    std::string CodeIni = "BDF1/DC"; // Define the way of initialization |"SolEx"->for exact
                                      // solution|  |"BDF1/DCF"->using only initial conditions|
    // BDF1Solver solver(tol, metaNumber, linearSystem, sol, norms, mesh, t0, t1, nTimeSteps);
    BDF2Solver solver(tol, metaNumber, linearSystem, sol, norms, mesh, t0, t1, nTimeSteps, CodeIni);
    solver.makeSteps(nTimeSteps, fespace);
    


    std::vector<double> &normLu = solver.getNorm(0);
    std::vector<double> &normLv = solver.getNorm(1);

    std::cout<<"Int LU"<<std::endl;
    for(int i = 0; i < normLu.size(); i++){
      std::cout<<normLu[i]<<std::endl;
    }

    std::cout<<"Int LV"<<std::endl;
    for(int i = 0; i < normLu.size(); i++){
      std::cout<<normLv[i]<<std::endl;
    }

    intL_u[2 * iter] =
      *std::max_element(normLu.begin(), normLu.end());    
    intL_v[2 * iter] =
      *std::max_element(normLv.begin(), normLv.end()); 


#endif
    std::string vtkFile = "../../data/TestNS_DC3.vtk";
    feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);
    delete linearSystem;
    delete NS2D;
    delete sol;
    delete metaNumber;
    delete mesh;
#endif
  }

  delete funSolP;
  delete funZero;

  printf("\n");
  printf("\n");
  printf("%12s \t %12s \t %12s \t %12s \n", "nSteps", "nElm",
         "||L||_BDF2",  "||L||_DC3");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e \n", TT[i], nElm[i],
            intL_u[2 * i], intL_v[2 * i]);
  printf("\n");
  printf("\n");

#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}
