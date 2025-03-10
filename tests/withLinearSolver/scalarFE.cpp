
#include "feAPI.h"
#include "../tests.h"

#include <gtest/gtest.h>

/*
  Tests the implementation of scalar FE spaces and weak forms
*/

static int my_argc;
static char** my_argv;

namespace diffusion {

  double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    return sin(M_PI*x)*sin(M_PI*y);
  }

  double uSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    return -2. * M_PI*M_PI * sin(M_PI*x)*sin(M_PI*y);
  }

  feStatus solve(const std::string &meshFile,
                 const int order,
                 const int degreeQuadrature,
                 const feNLSolverOptions &NLoptions,
                 int &numInteriorElements,
                 double &L2Error)
  {
    feStatus s;
    feFunction uSol(uSol_f);
    feFunction uSrc(uSrc_f);

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();

    feSpace *u = nullptr, *uBord = nullptr;
    s = createFiniteElementSpace(uBord, &mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol);
    if(s != FE_STATUS_OK) return s;
    s = createFiniteElementSpace(    u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &scalarConstant::zero);
    if(s != FE_STATUS_OK) return s;

    std::vector<feSpace*> spaces = {uBord, u};
    std::vector<feSpace*> essentialSpaces = {uBord};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    feBilinearForm *diff = nullptr, *src = nullptr;
    s = createBilinearForm(diff, {u}, new feSysElm_Diffusion<2>(&scalarConstant::one));
    if(s != FE_STATUS_OK) return s;
    s = createBilinearForm(src, {u}, new feSysElm_Source(&uSrc));
    if(s != FE_STATUS_OK) return s;
    std::vector<feBilinearForm*> forms = {diff, src};

    feLinearSystem *system;
    s = createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns());
    if(s != FE_STATUS_OK) return s;

    feNorm *errorU_L2 = nullptr;
    s = createNorm(errorU_L2, L2_ERROR, {u}, &sol, &uSol);
    if(s != FE_STATUS_OK) return s;
    std::vector<feNorm *> norms = {};

    TimeIntegrator *solver;
    s = createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, {nullptr, 1, ""});
    if(s != FE_STATUS_OK) return s;
    s = solver->makeSteps(1);
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

  feStatus meshConvergence(std::stringstream &resultBuffer, int order, int numMeshes, int degreeQuadrature)
  {
    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};
    {
      std::vector<int> nElm(numMeshes);
      std::vector<double> err(numMeshes, 0.);
      for(int i = 0; i < numMeshes; ++i)
      {
        std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
        feStatus s = solve(meshFile, order, degreeQuadrature, NLoptions, nElm[i], err[i]);
        if(s != FE_STATUS_OK) return s;
      }
      resultBuffer << "Diffusion - Error on u - Lagrange elements P" << order << std::endl;
      computeAndPrintConvergence(2, numMeshes, err, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
    }
    return FE_STATUS_OK;
  }
}

namespace advectionDiffusion {

  double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    return sin(x*y);
  }

  void velocity_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    UNUSED(x,y);
    res[0] = par[0];
    res[1] = par[1];
  }

  thread_local std::vector<double> velocityVec(2, 0.);

  double uSrc_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double k = par[2];
    velocity_f(args, par, velocityVec);
    const double cx = velocityVec[0];
    const double cy = velocityVec[1];
    return - ((cx*y+cy*x) * cos(x*y) + k * (x*x+y*y)*sin(x*y));
  }

  feStatus solve(const std::string &meshFile,
                 const int order,
                 const int degreeQuadrature,
                 const feNLSolverOptions &NLoptions,
                 int &numInteriorElements,
                 double &L2Error)
  {
    feStatus s;

    double cx = 1., cy = 2., k = 1.;
    std::vector<double> param = {cx, cy, k};
    feFunction uSol(uSol_f);
    feFunction uSrc(uSrc_f, param);
    feVectorFunction velocity(velocity_f, param);
    feConstantFunction conductivity(k);

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();

    feSpace *u = nullptr, *uBord = nullptr;
    s = createFiniteElementSpace(uBord, &mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol);
    if(s != FE_STATUS_OK) return s;
    s = createFiniteElementSpace(    u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &scalarConstant::zero);
    if(s != FE_STATUS_OK) return s;

    std::vector<feSpace*> spaces = {uBord, u};
    std::vector<feSpace*> essentialSpaces = {uBord};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    feBilinearForm *adv = nullptr, *diff = nullptr, *src = nullptr;
    s = createBilinearForm(adv, {u}, new feSysElm_Advection<2>(&velocity));
    if(s != FE_STATUS_OK) return s;
    s = createBilinearForm(diff, {u}, new feSysElm_Diffusion<2>(&conductivity));
    if(s != FE_STATUS_OK) return s;
    s = createBilinearForm(src, {u}, new feSysElm_Source(&uSrc));
    if(s != FE_STATUS_OK) return s;
    std::vector<feBilinearForm*> forms = {adv, diff, src};
    // adv->setComputeMatrixWithFD(true);

    feLinearSystem *system;
    s = createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns());
    if(s != FE_STATUS_OK) return s;

    feNorm *errorU_L2 = nullptr;
    s = createNorm(errorU_L2, L2_ERROR, {u}, &sol, &uSol);
    if(s != FE_STATUS_OK) return s;
    std::vector<feNorm *> norms = {};

    TimeIntegrator *solver;
    s = createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, {nullptr, 1, ""});
    if(s != FE_STATUS_OK) return s;
    s = solver->makeSteps(1);
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

  feStatus meshConvergence(std::stringstream &resultBuffer, int order, int numMeshes, int degreeQuadrature)
  {
    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};
    {
      std::vector<int> nElm(numMeshes);
      std::vector<double> err(numMeshes, 0.);
      for(int i = 0; i < numMeshes; ++i)
      {
        std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
        feStatus s = solve(meshFile, order, degreeQuadrature, NLoptions, nElm[i], err[i]);
        if(s != FE_STATUS_OK) return s;
      }
      resultBuffer << "Advection-diffusion - Error on u - Lagrange elements P" << order << std::endl;
      computeAndPrintConvergence(2, numMeshes, err, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
    }
    return FE_STATUS_OK;
  }
}

namespace advectionDiffusionReaction {

  double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    return sin(x*y);
  }

  void velocity_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    UNUSED(x,y);
    res[0] = par[0];
    res[1] = par[1];
  }

  thread_local std::vector<double> velocityVec(2, 0.);

  double uSrc_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double k = par[2];
    const double alpha = par[3];
    velocity_f(args, par, velocityVec);
    const double cx = velocityVec[0];
    const double cy = velocityVec[1];
    return - ((cx*y+cy*x) * cos(x*y) + k * (x*x+y*y)*sin(x*y) + alpha*sin(x*y));
  }

  feStatus solve(const std::string &meshFile,
                 const int order,
                 const int degreeQuadrature,
                 const feNLSolverOptions &NLoptions,
                 int &numInteriorElements,
                 double &L2Error)
  {
    feStatus s;

    double cx = 1., cy = 2., k = 1., alpha = 7.;
    std::vector<double> param = {cx, cy, k, alpha};
    feFunction uSol(uSol_f);
    feFunction uSrc(uSrc_f, param);
    feVectorFunction velocity(velocity_f, param);
    feConstantFunction conductivity(k);
    feConstantFunction reactionCoeff(alpha);

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();

    feSpace *u = nullptr, *uBord = nullptr;
    s = createFiniteElementSpace(uBord, &mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol);
    if(s != FE_STATUS_OK) return s;
    s = createFiniteElementSpace(    u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &scalarConstant::zero);
    if(s != FE_STATUS_OK) return s;

    std::vector<feSpace*> spaces = {uBord, u};
    std::vector<feSpace*> essentialSpaces = {uBord};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    feBilinearForm *adv = nullptr, *diff = nullptr, *mass = nullptr, *src = nullptr;
    s = createBilinearForm(adv, {u}, new feSysElm_Advection<2>(&velocity));
    if(s != FE_STATUS_OK) return s;
    s = createBilinearForm(diff, {u}, new feSysElm_Diffusion<2>(&conductivity));
    if(s != FE_STATUS_OK) return s;
    s = createBilinearForm(mass, {u}, new feSysElm_Mass(&reactionCoeff));
    if(s != FE_STATUS_OK) return s;
    s = createBilinearForm(src, {u}, new feSysElm_Source(&uSrc));
    if(s != FE_STATUS_OK) return s;
    std::vector<feBilinearForm*> forms = {adv, diff, mass, src};

    feLinearSystem *system;
    s = createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns());
    if(s != FE_STATUS_OK) return s;

    feNorm *errorU_L2 = nullptr;
    s = createNorm(errorU_L2, L2_ERROR, {u}, &sol, &uSol);
    if(s != FE_STATUS_OK) return s;
    std::vector<feNorm *> norms = {};

    TimeIntegrator *solver;
    s = createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, {nullptr, 1, ""});
    if(s != FE_STATUS_OK) return s;
    s = solver->makeSteps(1);
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

  feStatus meshConvergence(std::stringstream &resultBuffer, int order, int numMeshes, int degreeQuadrature)
  {
    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};
    {
      std::vector<int> nElm(numMeshes);
      std::vector<double> err(numMeshes, 0.);
      for(int i = 0; i < numMeshes; ++i)
      {
        std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
        feStatus s = solve(meshFile, order, degreeQuadrature, NLoptions, nElm[i], err[i]);
        if(s != FE_STATUS_OK) return s;
      }
      resultBuffer << "Advection-diffusion-reaction - Error on u - Lagrange elements P" << order << std::endl;
      computeAndPrintConvergence(2, numMeshes, err, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
    }
    return FE_STATUS_OK;
  }
}

TEST(ScalarFE, Diffusion)
{
  initialize(my_argc, my_argv);
  setVerbose(0);
  std::string testRoot = "../../../tests/withLinearSolver/scalarFE_diffusion";
  std::stringstream resultBuffer;
  int degreeQuadrature = 8;
  ASSERT_TRUE(diffusion::meshConvergence(resultBuffer, 1, 4, degreeQuadrature) == FE_STATUS_OK);
  ASSERT_TRUE(diffusion::meshConvergence(resultBuffer, 2, 4, degreeQuadrature) == FE_STATUS_OK);
  ASSERT_TRUE(diffusion::meshConvergence(resultBuffer, 3, 4, degreeQuadrature) == FE_STATUS_OK);
  ASSERT_TRUE(diffusion::meshConvergence(resultBuffer, 4, 4, degreeQuadrature) == FE_STATUS_OK);
  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

TEST(ScalarFE, AdvectionDiffusion)
{
  initialize(my_argc, my_argv);
  setVerbose(0);
  std::string testRoot = "../../../tests/withLinearSolver/scalarFE_advdiff";
  std::stringstream resultBuffer;
  int degreeQuadrature = 8;
  ASSERT_TRUE(advectionDiffusion::meshConvergence(resultBuffer, 1, 4, degreeQuadrature) == FE_STATUS_OK);
  ASSERT_TRUE(advectionDiffusion::meshConvergence(resultBuffer, 2, 4, degreeQuadrature) == FE_STATUS_OK);
  ASSERT_TRUE(advectionDiffusion::meshConvergence(resultBuffer, 3, 4, degreeQuadrature) == FE_STATUS_OK);
  ASSERT_TRUE(advectionDiffusion::meshConvergence(resultBuffer, 4, 4, degreeQuadrature) == FE_STATUS_OK);
  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

TEST(ScalarFE, AdvectionDiffusionReaction)
{
  initialize(my_argc, my_argv);
  setVerbose(0);
  std::string testRoot = "../../../tests/withLinearSolver/scalarFE_advdiffreact";
  std::stringstream resultBuffer;
  int degreeQuadrature = 8;
  ASSERT_TRUE(advectionDiffusionReaction::meshConvergence(resultBuffer, 1, 4, degreeQuadrature) == FE_STATUS_OK);
  ASSERT_TRUE(advectionDiffusionReaction::meshConvergence(resultBuffer, 2, 4, degreeQuadrature) == FE_STATUS_OK);
  ASSERT_TRUE(advectionDiffusionReaction::meshConvergence(resultBuffer, 3, 4, degreeQuadrature) == FE_STATUS_OK);
  ASSERT_TRUE(advectionDiffusionReaction::meshConvergence(resultBuffer, 4, 4, degreeQuadrature) == FE_STATUS_OK);
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