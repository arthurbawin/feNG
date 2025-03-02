
#include "feAPI.h"
#include "../tests.h"

#include <gtest/gtest.h>

/*
  Tests the implementation of vector FE spaces and weak forms :

  - weakGradient : Solve u = grad(f)
  - vectorLaplacian : Solve lap(u) = f with strong or weak Dirichlet conditions

*/

static int my_argc;
static char** my_argv;

// Weak gradient test case
//
// Solve u = grad(f) on \Omega, where f is a scalar source term
//       
//       u = uExact  on \partial \Omega
//
namespace weakGradient {

  void uSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double a = par[0];
    const double b = par[1];
    const double c = par[2];
    const double d = par[3];
    res[0] = a*(x*y + cos(x)) + b*y + 3.*c*x*x;
    res[1] = a/2.*x*x + b*x + 3.*d*y*y + 2.*d*y*cos(y*y);
  }

  double uScalarSrc_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double a = par[0];
    const double b = par[1];
    const double c = par[2];
    const double d = par[3];
    return - (a/2.*x*x*y + a*sin(x) + b*x*y + c*x*x*x + d*y*y*y + d*sin(y*y));
  }

  int solve(const std::string &meshFile,
            const int order,
            const int degreeQuadrature,
            const feNLSolverOptions &NLoptions,
            int &numInteriorElements,
            double &L2Error)
  {
    double a = 0.1234;
    double b = 1.2345;
    double c = 2.3456;
    double d = 3.4567;
    feVectorFunction uSol(uSol_f, {a, b, c, d});
    feFunction scalarSource(uScalarSrc_f, {a, b, c, d});

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();

    feSpace *u = nullptr, *uBord = nullptr;
    feCheck(createFiniteElementSpace(uBord, &mesh, elementType::VECTOR_LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol));
    feCheck(createFiniteElementSpace(    u, &mesh, elementType::VECTOR_LAGRANGE, order, "U", "Domaine", degreeQuadrature, &vectorConstant::zero));
    
    std::vector<feSpace*> spaces = {uBord, u};
    std::vector<feSpace*> essentialSpaces = {uBord};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    feBilinearForm *mass = nullptr, *gradsource = nullptr;
    feCheck(createBilinearForm(      mass, {u}, new feSysElm_VectorMass(&scalarConstant::one)));
    feCheck(createBilinearForm(gradsource, {u}, new feSysElm_GradSource(&scalarSource)));
    std::vector<feBilinearForm*> forms = {mass, gradsource};

    feLinearSystem *system;
    feCheck(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));

    feNorm *errorU_L2 = nullptr;
    feCheck(createNorm(errorU_L2, VECTOR_L2_ERROR, {u}, &sol, nullptr, &uSol));
    std::vector<feNorm *> norms = {};

    // feExporter *exporter;
    // feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
    // feExportData exportData = {exporter, 1, "vec"};

    TimeIntegrator *solver;
    feCheck(createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, {nullptr, 1, ""}));
    feCheck(solver->makeSteps(1));

    L2Error = errorU_L2->compute();

    delete errorU_L2;
    delete solver;
    delete system;
    for(feBilinearForm* f : forms)
      delete f;
    for(feSpace *s : spaces)
      delete s;

    return 0;
  }

  int meshConvergence(std::stringstream &resultBuffer, int order, int numMeshes, int degreeQuadrature)
  {
    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};
    {
      std::vector<int> nElm(numMeshes);
      std::vector<double> err(numMeshes, 0.);
      for(int i = 0; i < numMeshes; ++i)
      {
        std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
        solve(meshFile, order, degreeQuadrature, NLoptions, nElm[i], err[i]);
      }
      resultBuffer << "Weak gradient - Error on u - Lagrange elements P" << order << std::endl;
      computeAndPrintConvergence(2, numMeshes, err, nElm, resultBuffer);
    }
    return 0;
  }
}

// Vector Laplacian test case
//
// Solve u = grad(f) on \Omega, where f is a scalar source term
//       
//       u = uD      on \partial \Omega
//
// The Dirichlet BC can be enforced strongly, or weakly,
// in which case we solve instead the mixed problem for (u, lambda) := (u, -dudn):
//
// int_Omega grad(u):grad(v) + int_GammaD lambda*v + ... = 0,
//
//                       int_GammaD phiLambda * (uD - u) = 0
namespace vectorLaplacian {

  void uSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double a = par[0];
    const double b = par[1];
    res[0] = a * pow(x, 6) + pow(x,5) * y*y;
    res[1] = b * sin(x*y);
  }

  // Normal derivative on the boundary
  // Using tests on the position (assuming unit square domain) to check on which side we are
  void minus_dudn_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double a = par[0];
    const double b = par[1];
    const double grad_u[2][2] = {{6. * a * pow(x,5) + 5. * pow(x,4) * y*y, 2. * pow(x,5)*y},
                                 {                           b*y*cos(x*y),    b*x*cos(x*y)}};
    double n[2];
    if(fabs(x) < 1e-10) {
      n[0] = -1.;
      n[1] =  0.;
    } else if(fabs(x - 1.) < 1e-10) {
      n[0] = +1.;
      n[1] =  0.;
    } else if(fabs(y) < 1e-10) {
      n[0] =  0.;
      n[1] = -1.;
    } else {
      n[0] =  0.;
      n[1] = +1.;
    }
    res[0] = - (grad_u[0][0] * n[0] + grad_u[0][1] * n[1]);
    res[1] = - (grad_u[1][0] * n[0] + grad_u[1][1] * n[1]);
  }

  void uSrc_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double a = par[0];
    const double b = par[1];
    res[0] = + (30. * a * pow(x, 4) + 20. * x*x*x*y*y + 2. * pow(x, 5));
    res[1] = - b * sin(x*y) * (x*x + y*y);
  }

  int solve(const std::string &meshFile,
            const int order,
            const int degreeQuadrature,
            const feNLSolverOptions &NLoptions,
            int &numInteriorElements,
            double &L2Error,
            double &L2ErrorSecondaryVariable,
            const bool imposeDirichletBCWeakly)
  {
    double a = 0.1234;
    double b = 1.2345;
    feVectorFunction uSol(uSol_f, {a, b});
    feVectorFunction uSrc(uSrc_f, {a, b});
    feVectorFunction minus_dudn(minus_dudn_f, {a, b});

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();

    feSpace *u = nullptr, *uB = nullptr;
    feSpace *lB = nullptr;
    feCheck(createFiniteElementSpace( u, &mesh, elementType::VECTOR_LAGRANGE, order, "U","Domaine", degreeQuadrature, &vectorConstant::zero));
    feCheck(createFiniteElementSpace(uB, &mesh, elementType::VECTOR_LAGRANGE, order, "U",   "Bord", degreeQuadrature, &uSol));

    std::vector<feSpace*> spaces = {u, uB};
    std::vector<feSpace*> essentialSpaces = {};

    if(imposeDirichletBCWeakly)
    {
      // Weak imposition
      feCheck(createFiniteElementSpace(lB, &mesh, elementType::VECTOR_LAGRANGE, order-1, "L",    "Bord", degreeQuadrature, &vectorConstant::zero));
      spaces.push_back(lB);
    } else {
      // Strong imposition of Dirichlet BC
      UNUSED(lB);
      essentialSpaces.push_back(uB);
    }

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

    feBilinearForm *diff = nullptr, *source = nullptr;
    feCheck(createBilinearForm(  diff, {u}, new feSysElm_VectorDiffusion(&scalarConstant::one, &scalarConstant::one) ));
    feCheck(createBilinearForm(source, {u}, new feSysElm_VectorSource(&uSrc) ));
    std::vector<feBilinearForm*> forms = {diff, source};

    feBilinearForm *massUL = nullptr, *massLU = nullptr, *srcLambda = nullptr;
    if(imposeDirichletBCWeakly)
    {
      feCheck(createBilinearForm(   massUL, {uB, lB}, new feSysElm_MixedVectorMass(&scalarConstant::one)));
      feCheck(createBilinearForm(   massLU, {lB, uB}, new feSysElm_MixedVectorMass(&scalarConstant::minusOne)));
      feCheck(createBilinearForm(srcLambda, {lB}, new feSysElm_VectorSource(&uSol)));
      forms.push_back(massUL);
      forms.push_back(massLU);
      forms.push_back(srcLambda);
    } else {
      UNUSED(massUL, massLU, srcLambda);
    }

    feLinearSystem *system;
    feCheck(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));

    feNorm *errorU_L2 = nullptr;
    feCheck(createNorm(errorU_L2, VECTOR_L2_ERROR, {u}, &sol, nullptr, &uSol));
    std::vector<feNorm *> norms = {};

    // feExporter *exporter;
    // feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
    // feExportData exportData = {exporter, 1, "vec"};

    TimeIntegrator *solver;
    feCheck(createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, {nullptr, 1, ""}));
    feCheck(solver->makeSteps(1));

    L2Error = errorU_L2->compute();

    if(imposeDirichletBCWeakly)
    {
      // Evaluate error on the secondary variable lambda = -du/dn
      feNorm *errorL_L2 = nullptr;
      feCheck(createNorm(errorL_L2, VECTOR_L2_ERROR, {lB}, &sol, nullptr, &minus_dudn));
      L2ErrorSecondaryVariable = errorL_L2->compute();
    }

    delete errorU_L2;
    delete solver;
    delete system;
    for(feBilinearForm* f : forms)
      delete f;
    for(feSpace *s : spaces)
      delete s;

    return 0;
  }

  int meshConvergence(std::stringstream &resultBuffer, int order, int numMeshes, int degreeQuadrature)
  {
    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};

    std::vector<int> nElm(numMeshes);
    std::vector<double> err(numMeshes, 0.);
    std::vector<double> errSecondary(numMeshes, 0.);

    // Strongly enforced Dirichlet BC
    bool imposeDirichletBCWeakly = false;
    for(int i = 0; i < numMeshes; ++i)
    {
      std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
      solve(meshFile, order, degreeQuadrature, NLoptions, nElm[i], err[i], errSecondary[i], imposeDirichletBCWeakly);
    }
    resultBuffer << "Vector Laplacian - Strong Dirichlet BC - Error on primary variable - Lagrange elements P" << order << std::endl;
    computeAndPrintConvergence(2, numMeshes, err, nElm, resultBuffer);

    // Weakly enforced Dirichlet BC
    if(order > 1) {
      imposeDirichletBCWeakly = true;
      for(int i = 0; i < numMeshes; ++i)
      {
        std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
        solve(meshFile, order, degreeQuadrature, NLoptions, nElm[i], err[i], errSecondary[i], imposeDirichletBCWeakly);
      }
      resultBuffer << "Vector Laplacian - Weak Dirichlet BC - Error on primary variable - Lagrange elements P" << order << std::endl;
      computeAndPrintConvergence(2, numMeshes, err, nElm, resultBuffer);
      resultBuffer << "Vector Laplacian - Weak Dirichlet BC - Error on secondary variable - Lagrange elements P" << order << std::endl;
      computeAndPrintConvergence(2, numMeshes, errSecondary, nElm, resultBuffer);
    }

    return 0;
  }
}

//
// Weak gradient test case
// Tests the implementation of :
//    - vector-valued mass matrix (dot product between vector-valued test functions)
//    - divergence of vector-valued test functions (GradSource)
//
TEST(VectorFE, WeakGradient)
{
  initialize(my_argc, my_argv);
  setVerbose(0);
  
  std::string testRoot = "../../../tests/withLinearSolver/vectorFE_weakGradient";

  std::stringstream resultBuffer;

  int degreeQuadrature = 8;
  weakGradient::meshConvergence(resultBuffer, 1, 4, degreeQuadrature);
  weakGradient::meshConvergence(resultBuffer, 2, 4, degreeQuadrature);
  weakGradient::meshConvergence(resultBuffer, 3, 4, degreeQuadrature);
  weakGradient::meshConvergence(resultBuffer, 4, 4, degreeQuadrature);

  // Compare the resultBuffer with the file testRoot.output
  // If output file does not exist, a temporary file .tmp_output
  // is created, which may be used as the new .output file for this test.
  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);

  finalize();
}

TEST(VectorFE, VectorLaplacian)
{
  initialize(my_argc, my_argv);
  setVerbose(0);
  
  std::string testRoot = "../../../tests/withLinearSolver/vectorFE_vectorLaplacian";

  std::stringstream resultBuffer;

  int degreeQuadrature = 8;
  vectorLaplacian::meshConvergence(resultBuffer, 1, 4, degreeQuadrature);
  vectorLaplacian::meshConvergence(resultBuffer, 2, 4, degreeQuadrature);
  vectorLaplacian::meshConvergence(resultBuffer, 3, 4, degreeQuadrature);
  vectorLaplacian::meshConvergence(resultBuffer, 4, 4, degreeQuadrature);

  // Compare the resultBuffer with the file testRoot.output
  // If output file does not exist, a temporary file .tmp_output
  // is created, which may be used as the new .output file for this test.
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