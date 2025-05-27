
#include "feAPI.h"
#include "../tests.h"

#include <gtest/gtest.h>

/*
  Tests the implementation of scalar FE spaces and weak forms
*/

static int my_argc;
static char** my_argv;

// struct diffusion {

//   static double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
//   {
//     const double x = args.pos[0];
//     const double y = args.pos[1];
//     return sin(M_PI*x)*sin(M_PI*y);
//   }

//   static double uSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
//   {
//     const double x = args.pos[0];
//     const double y = args.pos[1];
//     return -2. * M_PI*M_PI * sin(M_PI*x)*sin(M_PI*y);
//   }

//   static feStatus solve(const std::string &meshFile,
//                         const int order,
//                         const int degreeQuadrature,
//                         const feNLSolverOptions &NLoptions,
//                         int &numInteriorElements,
//                         double &L2Error)
//   {
//     feFunction uSol(uSol_f);
//     feFunction uSrc(uSrc_f);

//     feMesh2DP1 mesh(meshFile);
//     numInteriorElements = mesh.getNumInteriorElements();

//     feSpace *u = nullptr, *uBord = nullptr;
//     feCheckReturn(createFiniteElementSpace(uBord, &mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol));
//     feCheckReturn(createFiniteElementSpace(    u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &scalarConstant::zero));

//     std::vector<feSpace*> spaces = {uBord, u};
//     std::vector<feSpace*> essentialSpaces = {uBord};

//     feMetaNumber numbering(&mesh, spaces, essentialSpaces);
//     feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
//     feBilinearForm *diff = nullptr, *src = nullptr;
//     feCheckReturn(createBilinearForm(diff, {u}, new feSysElm_Diffusion<2>(&scalarConstant::one)));
//     feCheckReturn(createBilinearForm(src, {u}, new feSysElm_Source(&uSrc)));
//     std::vector<feBilinearForm*> forms = {diff, src};

//     feLinearSystem *system;
//     #if defined(HAVE_MKL)
//       feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
//     #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
//       feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
//     #else
//       feCheckReturn(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
//     #endif

//     feNorm *errorU_L2 = nullptr;
//     feCheckReturn(createNorm(errorU_L2, L2_ERROR, {u}, &sol, &uSol));
//     std::vector<feNorm *> norms = {};

//     TimeIntegrator *solver;
//     feCheckReturn(createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, {nullptr, 1, ""}));
//     feCheckReturn(solver->makeSteps(1));

//     L2Error = errorU_L2->compute();

//     delete errorU_L2;
//     delete solver;
//     delete system;
//     for(feBilinearForm* f : forms)
//       delete f;
//     for(feSpace *space : spaces)
//       delete space;

//     return FE_STATUS_OK;
//   }
// };

// struct nonLinearDiffusion {

//   static double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
//   {
//     const double x = args.pos[0];
//     const double y = args.pos[1];
//     return sin(M_PI*x)*sin(M_PI*y);
//   }

//   static double uSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
//   {
//     const double x = args.pos[0];
//     const double y = args.pos[1];
//     const double cpx = cos(M_PI*x);
//     const double cpy = cos(M_PI*y);
//     const double spx = sin(M_PI*x);
//     const double spy = sin(M_PI*y);
//     return 2.*M_PI*M_PI* (cpx*cpx*spx*spy*spy*spy + cpy*cpy*spx*spx*spx*spy - spx*spy*(spx*spx*spy*spy + 1.));
//   }

//   static double conductivity_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
//   {
//     const double u = args.u;
//     return 1. + u*u;
//   }

//   static double ddu_conductivity_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
//   {
//     const double u = args.u;
//     return 2.*u;
//   }

//   static feStatus solve(const std::string &meshFile,
//                  const int order,
//                  const int degreeQuadrature,
//                  const feNLSolverOptions &NLoptions,
//                  int &numInteriorElements,
//                  double &L2Error)
//   {
//     feFunction uSol(uSol_f);
//     feFunction uSrc(uSrc_f);
//     feFunction conductivity(conductivity_f);
//     feFunction ddu_conductivity(ddu_conductivity_f);

//     feMesh2DP1 mesh(meshFile);
//     numInteriorElements = mesh.getNumInteriorElements();

//     feSpace *u = nullptr, *uBord = nullptr;
//     feCheckReturn(createFiniteElementSpace(uBord, &mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol));
//     feCheckReturn(createFiniteElementSpace(    u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &scalarConstant::zero));

//     std::vector<feSpace*> spaces = {uBord, u};
//     std::vector<feSpace*> essentialSpaces = {uBord};

//     feMetaNumber numbering(&mesh, spaces, essentialSpaces);
//     feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
//     feBilinearForm *diff = nullptr, *src = nullptr;
//     feCheckReturn(createBilinearForm(diff, {u}, new feSysElm_NonlinearDiffusion<2>(&conductivity, &ddu_conductivity)));
//     feCheckReturn(createBilinearForm( src, {u}, new feSysElm_Source(&uSrc)));
//     std::vector<feBilinearForm*> forms = {diff, src};

//     feLinearSystem *system;
//     #if defined(HAVE_MKL)
//       feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
//     #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
//       feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
//     #else
//       feCheckReturn(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
//     #endif

//     feNorm *errorU_L2 = nullptr;
//     feCheckReturn(createNorm(errorU_L2, L2_ERROR, {u}, &sol, &uSol));
//     std::vector<feNorm *> norms = {};

//     TimeIntegrator *solver;
//     feCheckReturn(createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, {nullptr, 1, ""}));
//     feCheckReturn(solver->makeSteps(1));

//     L2Error = errorU_L2->compute();

//     delete errorU_L2;
//     delete solver;
//     delete system;
//     for(feBilinearForm* f : forms)
//       delete f;
//     for(feSpace *space : spaces)
//       delete space;

//     return FE_STATUS_OK;
//   }
// };

// thread_local std::vector<double> velocityVec(2, 0.);

// struct advectionDiffusion {

//   static double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
//   {
//     const double x = args.pos[0];
//     const double y = args.pos[1];
//     return sin(x*y);
//   }

//   static void velocity_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
//   {
//     const double x = args.pos[0];
//     const double y = args.pos[1];
//     UNUSED(x,y);
//     res[0] = par[0];
//     res[1] = par[1];
//   }

//   static double uSrc_f(const feFunctionArguments &args, const std::vector<double> &par)
//   {
//     const double x = args.pos[0];
//     const double y = args.pos[1];
//     const double k = par[2];
//     velocity_f(args, par, velocityVec);
//     const double cx = velocityVec[0];
//     const double cy = velocityVec[1];
//     return - ((cx*y+cy*x) * cos(x*y) + k * (x*x+y*y)*sin(x*y));
//   }

//   static feStatus solve(const std::string &meshFile,
//                  const int order,
//                  const int degreeQuadrature,
//                  const feNLSolverOptions &NLoptions,
//                  int &numInteriorElements,
//                  double &L2Error)
//   {
//     double cx = 1., cy = 2., k = 1.;
//     std::vector<double> param = {cx, cy, k};
//     feFunction uSol(uSol_f);
//     feFunction uSrc(uSrc_f, param);
//     feVectorFunction velocity(velocity_f, param);
//     feConstantFunction conductivity(k);

//     feMesh2DP1 mesh(meshFile);
//     numInteriorElements = mesh.getNumInteriorElements();

//     feSpace *u = nullptr, *uBord = nullptr;
//     feCheckReturn(createFiniteElementSpace(uBord, &mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol));
//     feCheckReturn(createFiniteElementSpace(    u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &scalarConstant::zero));

//     std::vector<feSpace*> spaces = {uBord, u};
//     std::vector<feSpace*> essentialSpaces = {uBord};

//     feMetaNumber numbering(&mesh, spaces, essentialSpaces);
//     feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
//     feBilinearForm *adv = nullptr, *diff = nullptr, *src = nullptr;
//     feCheckReturn(createBilinearForm(adv, {u}, new feSysElm_Advection<2>(&velocity)));
//     feCheckReturn(createBilinearForm(diff, {u}, new feSysElm_Diffusion<2>(&conductivity)));
//     feCheckReturn(createBilinearForm(src, {u}, new feSysElm_Source(&uSrc)));
//     std::vector<feBilinearForm*> forms = {adv, diff, src};
//     // adv->setComputeMatrixWithFD(true);

//     feLinearSystem *system;
//     #if defined(HAVE_MKL)
//       feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
//     #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
//       feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
//     #else
//       feCheckReturn(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
//     #endif

//     feNorm *errorU_L2 = nullptr;
//     feCheckReturn(createNorm(errorU_L2, L2_ERROR, {u}, &sol, &uSol));
//     std::vector<feNorm *> norms = {};

//     TimeIntegrator *solver;
//     feCheckReturn(createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, {nullptr, 1, ""}));
//     feCheckReturn(solver->makeSteps(1));

//     L2Error = errorU_L2->compute();

//     delete errorU_L2;
//     delete solver;
//     delete system;
//     for(feBilinearForm* f : forms)
//       delete f;
//     for(feSpace *space : spaces)
//       delete space;

//     return FE_STATUS_OK;
//   }
// };

struct advectionDiffusionReaction {

  static double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    return sin(x*y);
  }

  static void velocity_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    UNUSED(x,y);
    res[0] = par[0];
    res[1] = par[1];
  }

  // thread_local std::vector<double> velocityVec(2, 0.);

  static double uSrc_f(const feFunctionArguments &args, const std::vector<double> &par)
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

  static feStatus solve(const std::string &meshFile,
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
    #if defined(HAVE_MKL)
      s = createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns());
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      s = createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns());
    #else
      s = createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns());
    #endif
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
};

/*
  Tests the implementation of mixed forms for scalar-valued FE:
    - feSysElm_MixedMass
    - feSysElm_MixedMassPower
    - feSysElm_MixedGradGrad

  Solves a coupled system with unknowns u, v:

    u + v^3 - div(grad(v)) = 0
    u + v   - div(grad(u)) = 0
*/
struct mixedDiffusion {

  static double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    return sin(M_PI*x)*sin(M_PI*y);
  }

  static double vSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    return cos(M_PI*x)*cos(M_PI*y);
  }

  static double uSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x    = args.pos[0];
    const double y    = args.pos[1];
    const double u    = sin(M_PI*x)*sin(M_PI*y);
    const double v    = cos(M_PI*x)*cos(M_PI*y);
    const double lapV = - 2. * M_PI * M_PI * cos(M_PI*x)*cos(M_PI*y);
    return - (u + v*v*v + lapV);
  }

  static double vSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double u    = sin(M_PI*x)*sin(M_PI*y);
    const double v    = cos(M_PI*x)*cos(M_PI*y);
    const double lapU = - 2. * M_PI * M_PI * sin(M_PI*x)*sin(M_PI*y);
    return - (u + v + lapU);
  }

  static feStatus solve(const std::string &meshFile,
                        const int order,
                        const int degreeQuadrature,
                        const feNLSolverOptions &NLoptions,
                        int &numInteriorElements,
                        double &L2Error)
  {
    feFunction uSol(uSol_f);
    feFunction vSol(vSol_f);
    feFunction uSrc(uSrc_f);
    feFunction vSrc(vSrc_f);

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();

    feSpace *u = nullptr, *uB = nullptr;
    feSpace *v = nullptr, *vB = nullptr;
    feCheckReturn(createFiniteElementSpace(uB, &mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol));
    feCheckReturn(createFiniteElementSpace(vB, &mesh, elementType::LAGRANGE, order, "V",    "Bord", degreeQuadrature, &vSol));
    feCheckReturn(createFiniteElementSpace( u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &scalarConstant::zero));
    feCheckReturn(createFiniteElementSpace( v, &mesh, elementType::LAGRANGE, order, "V", "Domaine", degreeQuadrature, &scalarConstant::zero));

    std::vector<feSpace*> spaces = {uB, u, vB, v};
    std::vector<feSpace*> essentialSpaces = {uB, vB};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    feBilinearForm *u_uMass, *vCube, *v_uMass, *v_vMass, *diffU, *diffV, *srcU, *srcV;
    // Equation for u : u + lap(v) + Fu = 0
    feCheckReturn(createBilinearForm(u_uMass,    {u}, new feSysElm_Mass(&scalarConstant::one)));
    feCheckReturn(createBilinearForm(  vCube, {u, v}, new feSysElm_MixedMassPower(&scalarConstant::one, 3.)));
    feCheckReturn(createBilinearForm(  diffV, {u, v}, new feSysElm_MixedGradGrad<2>(&scalarConstant::minusOne)));
    feCheckReturn(createBilinearForm(   srcU,    {u}, new feSysElm_Source(&uSrc)));
    // Equation for v : v + lap(u) + Fv = 0
    feCheckReturn(createBilinearForm(v_vMass,    {v}, new feSysElm_Mass(&scalarConstant::one)));
    feCheckReturn(createBilinearForm(v_uMass, {v, u}, new feSysElm_MixedMass(&scalarConstant::one)));
    feCheckReturn(createBilinearForm(  diffU, {v, u}, new feSysElm_MixedGradGrad<2>(&scalarConstant::minusOne)));
    feCheckReturn(createBilinearForm(   srcV,    {v}, new feSysElm_Source(&vSrc)));
    std::vector<feBilinearForm*> forms = {u_uMass, vCube, v_vMass, v_uMass, diffU, diffV, srcU, srcV};

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
    #else
      feCheckReturn(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
    #endif

    feNorm *errorU_L2 = nullptr;
    feCheckReturn(createNorm(errorU_L2, L2_ERROR, {u}, &sol, &uSol));
    std::vector<feNorm *> norms = {};

    TimeIntegrator *solver;
    feCheckReturn(createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, {nullptr, 1, ""}));
    feCheckReturn(solver->makeSteps(1));

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
};

/*
  Tests the implementation of mixed forms for scalar-valued FE:

    - feSysElm_DiffusionFieldDependentCoeff
    - feSysElm_MixedGradGradFieldDependentCoeff

  Solves a coupled system with unknowns u, v and prescribed w:

    u - div(a(w) * grad(u)) = 0  : tests feSysElm_DiffusionFieldDependentCoeff
    v - div(a(w) * grad(u)) = 0  : tests feSysElm_MixedGradGradFieldDependentCoeff
    w - div(a(w) * grad(u)) = 0  : tests with coeff depending on field
                                   used to define the test functions

    with coefficient a(w) linear (C + w) or nonlinear (C + w^2)

  Be careful that as a diffusion coefficient, a(w) should remain (strictly) positive!
*/


struct mixedDiffusionFieldDependent {

  static double coeff_a(const double w)
  {
    return 2. + w;
    // return 2. + w*w;
    // return 2.; UNUSED(w);
  }

  static double dcoeff_a_dw(const double w)
  {
    return 1.; UNUSED(w);
    // return 2.*w;
    // return 0.; UNUSED(w);
  }

  static double coeff_b(const double val)
  {
    return 10. + val;
  }

  static double dcoeff_b_dval(const double val)
  {
    return 1.; UNUSED(val);
  }

  static double a_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    return coeff_a(args.u);
  }

  static double b_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    return coeff_b(args.u);
  }

  static double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    return sin(M_PI*x)*sin(M_PI*y);
  }

  static double vSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    return cos(M_PI*x)*cos(M_PI*y);
  }

  static double wSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    return sin(M_PI*x)*cos(M_PI*y);
  }

  static double uSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x     = args.pos[0];
    const double y     = args.pos[1];
    const double u     = sin(M_PI*x)*sin(M_PI*y);
    const double w     = sin(M_PI*x)*cos(M_PI*y);
    const double a     = coeff_a(w);
    const double dadw  = dcoeff_a_dw(w);
    const double lapU  = - 2. * M_PI * M_PI * sin(M_PI*x)*sin(M_PI*y);

    const double gradU[2] = {M_PI * cos(M_PI*x)*sin(M_PI*y),
                             M_PI * sin(M_PI*x)*cos(M_PI*y)};
    const double gradW[2] = {M_PI * cos(M_PI*x)*cos(M_PI*y),
                           - M_PI * sin(M_PI*x)*sin(M_PI*y)};
    const double gradW_dot_gradU = gradU[0] * gradW[0] + gradU[1] * gradW[1];

    return - (u - (a * lapU + dadw * gradW_dot_gradU));
  }

  static double vSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x     = args.pos[0];
    const double y     = args.pos[1];
    const double v     = cos(M_PI*x)*cos(M_PI*y);
    const double w     = sin(M_PI*x)*cos(M_PI*y);
    const double a     = coeff_a(w);
    const double dadw  = dcoeff_a_dw(w);
    const double lapU  = - 2. * M_PI * M_PI * sin(M_PI*x)*sin(M_PI*y);

    const double gradU[2] = {M_PI * cos(M_PI*x)*sin(M_PI*y),
                             M_PI * sin(M_PI*x)*cos(M_PI*y)};
    const double gradW[2] = {M_PI * cos(M_PI*x)*cos(M_PI*y),
                           - M_PI * sin(M_PI*x)*sin(M_PI*y)};
    const double gradW_dot_gradU = gradU[0] * gradW[0] + gradU[1] * gradW[1];

    return - (v - (a * lapU + dadw * gradW_dot_gradU));
  }

  static double wSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x     = args.pos[0];
    const double y     = args.pos[1];
    const double w     = sin(M_PI*x)*cos(M_PI*y);
    const double lapU  = - 2. * M_PI * M_PI * sin(M_PI*x)*sin(M_PI*y);

    const double gradU[2] = {M_PI * cos(M_PI*x)*sin(M_PI*y),
                             M_PI * sin(M_PI*x)*cos(M_PI*y)};
    const double gradV[2] = {- M_PI * sin(M_PI*x)*cos(M_PI*y),
                             - M_PI * cos(M_PI*x)*sin(M_PI*y)};
    const double gradW[2] = {M_PI * cos(M_PI*x)*cos(M_PI*y),
                           - M_PI * sin(M_PI*x)*sin(M_PI*y)};

    UNUSED(gradV, gradW);

    const double gradW_dot_gradU = gradU[0] * gradW[0] + gradU[1] * gradW[1];
    const double a     = coeff_a(w);
    const double dadw  = dcoeff_a_dw(w);
    return - (w - (a * lapU + dadw * gradW_dot_gradU));

    // const double gradV_dot_gradU = gradU[0] * gradV[0] + gradU[1] * gradV[1];
    // const double v     = cos(M_PI*x)*cos(M_PI*y);
    // const double a     = coeff_b(v);
    // const double dadv  = dcoeff_b_dval(v);
    // return - (w - (a * lapU + dadv * gradV_dot_gradU));
  }

  static feStatus solve(const std::string &meshFile,
                        const int order,
                        const int degreeQuadrature,
                        const feNLSolverOptions &NLoptions,
                        int &numInteriorElements,
                        std::vector<double> &errors)
  {
    feFunction uSol(uSol_f);
    feFunction vSol(vSol_f);
    feFunction wSol(wSol_f);
    feFunction uSrc(uSrc_f);
    feFunction vSrc(vSrc_f);
    feFunction wSrc(wSrc_f);

    feFunction aCoeff(a_f);
    feFunction bCoeff(b_f);

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();

    feSpace *u = nullptr, *uB = nullptr;
    feSpace *v = nullptr, *vB = nullptr;
    feSpace *w = nullptr, *wB = nullptr;
    feCheckReturn(createFiniteElementSpace(uB, &mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol));
    feCheckReturn(createFiniteElementSpace(vB, &mesh, elementType::LAGRANGE, order, "V",    "Bord", degreeQuadrature, &vSol));
    feCheckReturn(createFiniteElementSpace(wB, &mesh, elementType::LAGRANGE, order, "W",    "Bord", degreeQuadrature, &wSol));

    feCheckReturn(createFiniteElementSpace( u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &scalarConstant::zero));
    feCheckReturn(createFiniteElementSpace( v, &mesh, elementType::LAGRANGE, order, "V", "Domaine", degreeQuadrature, &scalarConstant::zero));
    feCheckReturn(createFiniteElementSpace( w, &mesh, elementType::LAGRANGE, order, "W", "Domaine", degreeQuadrature, &scalarConstant::zero));

    // feCheckReturn(createFiniteElementSpace( u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &uSol));
    // feCheckReturn(createFiniteElementSpace( v, &mesh, elementType::LAGRANGE, order, "V", "Domaine", degreeQuadrature, &vSol));
    // feCheckReturn(createFiniteElementSpace( w, &mesh, elementType::LAGRANGE, order, "W", "Domaine", degreeQuadrature, &wSol));

    std::vector<feSpace*> spaces = {uB, u, vB, v, wB, w};
    std::vector<feSpace*> essentialSpaces = {uB, vB, wB};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

    // u - div(a(w) * grad(u)) = 0  : tests feSysElm_DiffusionFieldDependentCoeff
    // v - div(a(w) * grad(u)) = 0  : tests feSysElm_MixedGradGradFieldDependentCoeff
    // w - div(a(w) * grad(u)) = 0  : tests with coeff depending on field
    
    feBilinearForm *uMass, *vMass, *wMass, *diffU_U, *diffU_V, *diffU_W, *srcU, *srcV, *srcW;
    // Equation for u
    feCheckReturn(createBilinearForm(  uMass,       {u}, new feSysElm_Mass(&scalarConstant::one)));
    feCheckReturn(createBilinearForm(diffU_U,    {u, w}, new feSysElm_DiffusionFieldDependentCoeff<2>(&aCoeff)));
    feCheckReturn(createBilinearForm(   srcU,       {u}, new feSysElm_Source(&uSrc)));
    // Equation for v
    feCheckReturn(createBilinearForm(  vMass,       {v}, new feSysElm_Mass(&scalarConstant::one)));
    feCheckReturn(createBilinearForm(diffU_V, {v, u, w}, new feSysElm_MixedGradGradFieldDependentCoeff<2>(&aCoeff)));
    feCheckReturn(createBilinearForm(   srcV,       {v}, new feSysElm_Source(&vSrc)));
    // Equation for w
    feCheckReturn(createBilinearForm(  wMass,       {w}, new feSysElm_Mass(&scalarConstant::one)));
    feCheckReturn(createBilinearForm(diffU_W, {w, u, w}, new feSysElm_MixedGradGradFieldDependentCoeff<2>(&aCoeff)));
    feCheckReturn(createBilinearForm(   srcW,       {w}, new feSysElm_Source(&wSrc)));

    // std::vector<feBilinearForm*> forms = {uMass, vMass, wMass, diffU_U, diffU_V, diffU_W, srcU, srcV, srcW};
    std::vector<feBilinearForm*> forms = {uMass, diffU_U, srcU, vMass, diffU_V, srcV, wMass, diffU_W, srcW};

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
    #else
      feCheckReturn(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
    #endif

    std::vector<feNorm *> norms = {};

    TimeIntegrator *solver;
    feCheckReturn(createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, {nullptr, 1, ""}));
    feCheckReturn(solver->makeSteps(1));

    feNorm *errorU_L2, *errorV_L2, *errorW_L2;
    feCheckReturn(createNorm(errorU_L2, L2_ERROR, {u}, &sol, &uSol));
    feCheckReturn(createNorm(errorV_L2, L2_ERROR, {v}, &sol, &vSol));
    feCheckReturn(createNorm(errorW_L2, L2_ERROR, {w}, &sol, &wSol));

    errors.resize(3);
    errors[0] = errorU_L2->compute();
    errors[1] = errorV_L2->compute();
    errors[2] = errorW_L2->compute();

    delete errorU_L2;
    delete errorV_L2;
    delete errorW_L2;
    delete solver;
    delete system;
    for(feBilinearForm* f : forms)
      delete f;
    for(feSpace *space : spaces)
      delete space;

    return FE_STATUS_OK;
  }
};

template<typename problemType>
feStatus meshConvergence(std::stringstream &resultBuffer,
                         const int order,
                         const int numMeshes,
                         const int degreeQuadrature,
                         const std::string &message)
{
  feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 20, 4, 1e-1};
  {
    std::vector<int> nElm(numMeshes);
    std::vector<std::vector<double>> err(numMeshes);

    for(int i = 0; i < numMeshes; ++i)
    {
      // std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
      // std::string meshFile = "../../../data/transfiniteSquare" + std::to_string(i+1) + ".msh";
      std::string meshFile = "../data/transfiniteSquare" + std::to_string(i+1) + ".msh";
      feCheckReturn(problemType::solve(meshFile, order, degreeQuadrature, NLoptions, nElm[i], err[i]));
    }

    for(auto vec : err)
      for (auto val : vec)
        feInfo("err = %f", val);

    size_t nFields = err[0].size();
    for(size_t i = 0; i < nFields; ++i)
    {
      std::vector<double> fieldError(numMeshes, 0.);
      for(int j = 0; j < numMeshes; ++j)
      {
        fieldError[j] = err[j][i];
      }
      resultBuffer << message << " field " << i << std::endl;
      // computeAndPrintConvergence(2, numMeshes, err, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      computeAndPrintConvergence(2, numMeshes, fieldError, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
    }
  }
  return FE_STATUS_OK;
}

// TEST(ScalarFE, Diffusion)
// {
//   initialize(my_argc, my_argv);
//   // setVerbose(0);
//   std::string testRoot = "../../../tests/withLinearSolver/scalarFE_diffusion";
//   std::stringstream resultBuffer;
//   int degreeQuadrature = 8;
//   for(int order = 1; order <= 4; ++order)
//   {
//     std::string message = "Diffusion - Error on u - Lagrange elements P" + std::to_string(order);
//     ASSERT_TRUE(meshConvergence<diffusion>(resultBuffer, order, 4, degreeQuadrature, message) == FE_STATUS_OK);
//   }
//   EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
//   finalize();
// }

// TEST(ScalarFE, NonLinearDiffusion)
// {
//   initialize(my_argc, my_argv);
//   // setVerbose(0);
//   std::string testRoot = "../../../tests/withLinearSolver/scalarFE_nonlineardiffusion";
//   std::stringstream resultBuffer;
//   int degreeQuadrature = 8;
//   for(int order = 1; order <= 4; ++order)
//   {
//     std::string message = "Diffusion - Error on u - Lagrange elements P" + std::to_string(order);
//     ASSERT_TRUE(meshConvergence<nonLinearDiffusion>(resultBuffer, order, 4, degreeQuadrature, message) == FE_STATUS_OK);
//   }
//   EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
//   finalize();
// }

// TEST(ScalarFE, AdvectionDiffusion)
// {
//   initialize(my_argc, my_argv);
//   setVerbose(0);
//   std::string testRoot = "../../../tests/withLinearSolver/scalarFE_advdiff";
//   std::stringstream resultBuffer;
//   int degreeQuadrature = 8;
//   for(int order = 1; order <= 4; ++order)
//   {
//     std::string message = "Advection-diffusion - Error on u - Lagrange elements P" + std::to_string(order);
//     ASSERT_TRUE(meshConvergence<advectionDiffusion>(resultBuffer, order, 4, degreeQuadrature, message) == FE_STATUS_OK);
//   }
//   EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
//   finalize();
// }

// TEST(ScalarFE, AdvectionDiffusionReaction)
// {
//   initialize(my_argc, my_argv);
//   setVerbose(0);
//   std::string testRoot = "../../../tests/withLinearSolver/scalarFE_advdiffreact";
//   std::stringstream resultBuffer;
//   int degreeQuadrature = 8;
//   for(int order = 1; order <= 4; ++order)
//   {
//     std::string message = "Advection-diffusion-reaction - Error on u - Lagrange elements P" + std::to_string(order);
//     ASSERT_TRUE(meshConvergence<advectionDiffusionReaction>(resultBuffer, order, 4, degreeQuadrature, message) == FE_STATUS_OK);
//   }
//   EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
//   finalize();
// }

// TEST(ScalarFE, MixedDiffusion)
// {
//   initialize(my_argc, my_argv);
//   setVerbose(1);
//   std::string testRoot = "../../../tests/withLinearSolver/scalarFE_mixedDiffusion";
//   std::stringstream resultBuffer;
//   int degreeQuadrature = 8;
//   for(int order = 1; order <= 4; ++order)
//   {
//     std::string message = "Mixed diffusion - Error on u - Lagrange elements P" + std::to_string(order);
//     ASSERT_TRUE(meshConvergence<mixedDiffusion>(resultBuffer, order, 4, degreeQuadrature, message) == FE_STATUS_OK);
//   }
//   EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
//   finalize();
// }

TEST(ScalarFE, MixedDiffusionFieldDependent)
{
  initialize(my_argc, my_argv);
  setVerbose(1);
  std::string testRoot = "../../../tests/withLinearSolver/scalarFE_mixedDiffusion";
  std::stringstream resultBuffer;
  int degreeQuadrature = 16;
  for(int order = 1; order <= 4; ++order)
  {
    std::string message = "Mixed diffusion - Error on u - Lagrange elements P" + std::to_string(order);
    ASSERT_TRUE(meshConvergence<mixedDiffusionFieldDependent>(resultBuffer, order, 4, degreeQuadrature, message) == FE_STATUS_OK);
  }
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