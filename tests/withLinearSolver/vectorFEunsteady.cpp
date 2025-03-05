
#include "feAPI.h"
#include "../tests.h"

#include <gtest/gtest.h>

/*
  Tests the implementation of scalar FE spaces and weak forms
*/

static int my_argc;
static char** my_argv;

namespace diffusion {

  void uSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    const double a = par[0];
    const double b = par[1];
    res[0] = exp(-t) * (a * pow(x, 6) + pow(x,5) * y*y);
    res[1] = exp(-t) * b * sin(x*y);
  }

  void uSrc_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    const double a = par[0];
    const double b = par[1];
    const double expmt = exp(-t);
    const double dudt[2] = {- expmt * (a * pow(x, 6) + pow(x,5) * y*y),
                            - expmt * b * sin(x*y)};
    const double lapu[2] = {expmt * (30. * a * pow(x, 4) + 20. * x*x*x*y*y + 2. * pow(x, 5)),
                            expmt * (- b * sin(x*y) * (x*x + y*y))};                        
    res[0] = - (dudt[0] - lapu[0]);
    res[1] = - (dudt[1] - lapu[1]);
  }

  feStatus solve(const std::string &meshFile,
                 const timeIntegratorScheme scheme,
                 const double t0,
                 const double t1,
                 const int nTimeSteps,
                 const int order,
                 const int degreeQuadrature,
                 const feNLSolverOptions &NLoptions,
                 int &numInteriorElements,
                 double &L2Error)
  {
    feStatus s;
    double a = 0.1234;
    double b = 1.2345;
    feVectorFunction uSol(uSol_f, {a, b});
    feVectorFunction uSrc(uSrc_f, {a, b});

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();

    feSpace *u = nullptr, *uBord = nullptr;
    s = createFiniteElementSpace(uBord, &mesh, elementType::VECTOR_LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol);
    if(s != FE_STATUS_OK) return s;
    s = createFiniteElementSpace(    u, &mesh, elementType::VECTOR_LAGRANGE, order, "U", "Domaine", degreeQuadrature, &vectorConstant::zero);
    if(s != FE_STATUS_OK) return s;

    std::vector<feSpace*> spaces = {uBord, u};
    std::vector<feSpace*> essentialSpaces = {uBord};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    feBilinearForm *dudt = nullptr, *diff = nullptr, *src = nullptr;
    s = createBilinearForm(dudt, {u}, new feSysElm_TransientVectorMass(&scalarConstant::one));
    if(s != FE_STATUS_OK) return s;
    s = createBilinearForm(diff, {u}, new feSysElm_VectorDiffusion(&scalarConstant::one, &scalarConstant::one));
    if(s != FE_STATUS_OK) return s;
    s = createBilinearForm( src, {u}, new feSysElm_VectorSource(&uSrc));
    if(s != FE_STATUS_OK) return s;
    std::vector<feBilinearForm*> forms = {dudt, diff, src};

    feLinearSystem *system;
    s = createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns());
    if(s != FE_STATUS_OK) return s;

    feNorm *errorU_L2 = nullptr;
    s = createNorm(errorU_L2, VECTOR_L2_ERROR, {u}, &sol, nullptr, &uSol);
    if(s != FE_STATUS_OK) return s;
    std::vector<feNorm *> norms = {};

    TimeIntegrator *solver;
    s = createTimeIntegrator(solver, scheme, NLoptions,
      system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps);
    if(s != FE_STATUS_OK) return s;
    s = solver->makeSteps(nTimeSteps);
    if(s != FE_STATUS_OK) return s;

    L2Error = errorU_L2->compute();

    delete errorU_L2;
    delete solver;
    delete system;
    for(feBilinearForm* f : forms)
      delete f;
    for(feSpace *space : spaces)
      delete space;

    return FE_STATUS_OK;
  }

  feStatus meshConvergence(std::stringstream &resultBuffer,
    timeIntegratorScheme scheme, double t0, double t1, 
    int initial_nTimeSteps, int order, int numMeshes, int degreeQuadrature)
  {
    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};
    {
      std::vector<int> nElm(numMeshes);
      std::vector<double> err(numMeshes, 0.);
      int nTimeSteps = initial_nTimeSteps;

      for(int i = 0; i < numMeshes; ++i, nTimeSteps *= 2)
      {
        std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
        feStatus s = solve(meshFile, scheme, t0, t1, nTimeSteps, order, degreeQuadrature, NLoptions, nElm[i], err[i]);
        if(s != FE_STATUS_OK) return s;
      }
      resultBuffer << "Unsteady vector diffusion - Error on u - Lagrange elements P" << order << std::endl;
      computeAndPrintConvergence(2, numMeshes, err, nElm, resultBuffer);
    }
    return FE_STATUS_OK;
  }
}

TEST(VectorFEUnsteady, DiffusionBDF1)
{
  initialize(my_argc, my_argv);
  setVerbose(0);
  std::string testRoot = "../../../tests/withLinearSolver/vectorFEunsteady_bdf1";
  std::stringstream resultBuffer;
  int degreeQuadrature = 8;
  double t0 = 0.;
  double t1 = 1.;
  int initial_nTimeSteps = 20;
  timeIntegratorScheme scheme = timeIntegratorScheme::BDF1;
  ASSERT_TRUE(diffusion::meshConvergence(resultBuffer, scheme, t0, t1, initial_nTimeSteps, 1, 4, degreeQuadrature) == FE_STATUS_OK);

  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

TEST(VectorFEUnsteady, DiffusionBDF2)
{
  initialize(my_argc, my_argv);
  setVerbose(0);
  std::string testRoot = "../../../tests/withLinearSolver/vectorFEunsteady_bdf2";
  std::stringstream resultBuffer;
  int degreeQuadrature = 8;
  double t0 = 0.;
  double t1 = 1.;
  int initial_nTimeSteps = 10;
  timeIntegratorScheme scheme = timeIntegratorScheme::BDF2;
  ASSERT_TRUE(diffusion::meshConvergence(resultBuffer, scheme, t0, t1, initial_nTimeSteps, 1, 4, degreeQuadrature) == FE_STATUS_OK);

  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

//
// Custom main to capture argc and argv
//
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  my_argc = argc;
  my_argv = argv;
  return RUN_ALL_TESTS();
}