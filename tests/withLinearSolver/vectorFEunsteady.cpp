
#include "feAPI.h"
#include "../tests.h"

#include <gtest/gtest.h>

/*
  Tests the implementation of unsteady weak forms for vector-valued spaces
*/

static int my_argc;
static char** my_argv;

static std::vector<timeIntegratorScheme> schemes = {timeIntegratorScheme::BDF1,
                                                    timeIntegratorScheme::BDF2,
                                                    timeIntegratorScheme::DC2BDF1,
                                                    timeIntegratorScheme::DC3BDF1};

static std::vector<std::string> schemesStr = {"BDF1", "BDF2", "DC2BDF1", "DC3BDF1"};

//
// Solve the linear vector-valued ODE : u'(t) = (-sin t, cos t)
//                                       u(0) = (1, 0)
// with solution u(t) = (cos t, sin t).
struct linearODE {

  static double FxT(const double t) { return cos(t); }
  static double FyT(const double t) { return sin(t); }
  static double DFxDT(const double t) { return -sin(t); }
  static double DFyDT(const double t) { return  cos(t); }

  // static double FxT(const double t) { return t*t*t*t; }
  // static double FyT(const double t) { return t*t; }
  // static double DFxDT(const double t) { return 4.*t*t*t; }
  // static double DFyDT(const double t) { return 2.*t; }

  // static double FxT(const double t) { return exp(-t); }
  // static double FyT(const double t) { return exp(-t); }
  // static double DFxDT(const double t) { return -exp(-t); }
  // static double DFyDT(const double t) { return -exp(-t); }

  static void uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/, std::vector<double> &res)
  {
    const double t = args.t;
    res[0] = FxT(t);
    res[1] = FyT(t);
  }

  static void uSolDot_f(const feFunctionArguments &args, const std::vector<double> &/*par*/, std::vector<double> &res)
  {
    const double t = args.t;
    res[0] = DFxDT(t);
    res[1] = DFyDT(t);
  }

  static void uSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/, std::vector<double> &res)
  {
    const double t = args.t;                      
    res[0] = - (DFxDT(t));
    res[1] = - (DFyDT(t));
  }

  static feStatus solve(const std::string &meshFile,
                 const timeIntegratorScheme scheme,
                 const double t0,
                 const double t1,
                 const int nTimeSteps,
                 const int order,
                 const int degreeQuadrature,
                 const feNLSolverOptions &NLoptions,
                 double &error)
  {
    feVectorFunction uSol(uSol_f), uSolDot(uSolDot_f), uSrc(uSrc_f);

    feMesh2DP1 mesh(meshFile);

    feSpace *u = nullptr;
    feCheckReturn(createFiniteElementSpace(u, &mesh, elementType::VECTOR_LAGRANGE, order, "U", "Domaine", degreeQuadrature, &uSol));
    std::vector<feSpace*> spaces = {u}, essentialSpaces = {};

    // feSpace *u = nullptr, *uB = nullptr;
    // feCheckReturn(createFiniteElementSpace(u, &mesh, elementType::VECTOR_LAGRANGE, order, "U", "Domaine", degreeQuadrature, &uSol));
    // feCheckReturn(createFiniteElementSpace(uB, &mesh, elementType::VECTOR_LAGRANGE, order, "U", "Bord", degreeQuadrature, &uSol));
    // std::vector<feSpace*> spaces = {u, uB}, essentialSpaces = {uB};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    feBilinearForm *dudt = nullptr, *src = nullptr;
    feCheckReturn(createBilinearForm(dudt, {u}, new feSysElm_TransientVectorMass<2>(&scalarConstant::one)));
    feCheckReturn(createBilinearForm( src, {u}, new feSysElm_VectorSource<2>(&uSrc)));
    std::vector<feBilinearForm*> forms = {dudt, src};

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
    #else
      feCheckReturn(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
    #endif

    feNorm *errorU = nullptr, *errurUdot = nullptr;
    feCheckReturn(createNorm(errorU, VECTOR_LINF_ERROR, {u}, &sol, nullptr, &uSol));
    feCheckReturn(createNorm(errurUdot, VECTOR_LINF_UDOT_ERROR, {u}, &sol, nullptr, &uSolDot));
    std::vector<feNorm *> norms = {errorU, errurUdot};

    TimeIntegrator *solver;
    feCheckReturn(createTimeIntegrator(solver, scheme, NLoptions,
      system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));

    // Start BDF2 with manufactured solution
    if(scheme == timeIntegratorScheme::BDF2)
      static_cast<BDF2Integrator*>(solver)->setStartingMethod(BDF2Starter::InitialCondition);

    feCheckReturn(solver->makeSteps(nTimeSteps));

    std::vector<std::vector<double>> postProc = solver->getPostProcessingData();
    error = 0.;
    for(int j = 0; j < nTimeSteps+1; ++j) {
      feInfo("Error at step %2d (t = %1.3e) = %1.6e (U) = %1.6e (Udot)", j, postProc[0][j], postProc[1][j], postProc[2][j]);
      error = fmax(error, postProc[1][j]);
    }

    delete errorU;
    delete solver;
    delete system;
    for(feBilinearForm* f : forms)
      delete f;
    for(feSpace *space : spaces)
      delete space;

    return FE_STATUS_OK;
  }
};

//
// Solve dudt - lap(u) + f = 0 for manufactured solution:
//
// u = exp(-t) * (ax^3 + xy^2, bxy).
//
// The "spatial" part of u is captured by cubic polynomials.
//
struct diffusion {

  // static double FT(const double t) { return t*t*t*t*t; }
  // static double DFDT(const double t) { return 5.*t*t*t*t; }

  // static double FT(const double t) { return exp(-t); }
  // static double DFDT(const double t) { return -exp(-t); }

  // static double FT(const double t) { return cos(t); }
  // static double DFDT(const double t) { return -sin(t); }

  static double FT(const double t) { return 2./(1.+t); }
  static double DFDT(const double t) { return -2./((1.+t)*(1.+t)); }

  static void uSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    const double a = par[0];
    const double b = par[1];
    // res[0] = exp(-t) * (a * pow(x, 6) + pow(x,5) * y*y);
    // res[1] = exp(-t) * b * sin(x*y);
    // res[0] = FT(t) * (a * pow(x, 3) + x * y*y);
    // res[1] = FT(t) * b * x*y;
    res[0] = FT(t) * (a*y + x);
    res[1] = FT(t) * b * x*y;
  }

  static void duSoldt_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    const double a = par[0];
    const double b = par[1];
    // res[0] = exp(-t) * (a * pow(x, 6) + pow(x,5) * y*y);
    // res[1] = exp(-t) * b * sin(x*y);
    // res[0] = DFDT(t) * (a * pow(x, 3) + x * y*y);
    // res[1] = DFDT(t) * b * x*y;
    res[0] = DFDT(t) * (a*y + x);
    res[1] = DFDT(t) * b * x*y;
  }

  static void uSrc_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    const double a = par[0];
    const double b = par[1];
    // const double expmt = exp(-t);
    // const double dudt[2] = {- expmt * (a * pow(x, 6) + pow(x,5) * y*y),
    //                         - expmt * b * sin(x*y)};
    // const double lapu[2] = {expmt * (30. * a * pow(x, 4) + 20. * x*x*x*y*y + 2. * pow(x, 5)),
    //                         expmt * (- b * sin(x*y) * (x*x + y*y))};                        
    // const double dudt[2] = {DFDT(t) * (a * pow(x, 3) + x * y*y),
    //                         DFDT(t) * b * x*y};
    // const double lapu[2] = {FT(t) * (6.*a*x + 2.*x), 0.};     
    const double dudt[2] = {DFDT(t) * (a * y + x),
                            DFDT(t) * b * x*y};
    const double lapu[2] = {0., 0.};                        
    res[0] = - (dudt[0] - lapu[0]);
    res[1] = - (dudt[1] - lapu[1]);
  }

  static feStatus solve(const std::string &meshFile,
                 const timeIntegratorScheme scheme,
                 const double t0,
                 const double t1,
                 const int nTimeSteps,
                 const int order,
                 const int degreeQuadrature,
                 const feNLSolverOptions &NLoptions,
                 double &error)
  {
    double a = 0.;
    double b = 0.;
    feVectorFunction uSol(uSol_f, {a, b});
    feVectorFunction duSoldt(duSoldt_f, {a, b});
    feVectorFunction uSrc(uSrc_f, {a, b});

    feMesh2DP1 mesh(meshFile);

    feSpace *u = nullptr, *uBord = nullptr;
    feCheckReturn(createFiniteElementSpace(uBord, &mesh, elementType::VECTOR_LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol));
    feCheckReturn(createFiniteElementSpace(    u, &mesh, elementType::VECTOR_LAGRANGE, order, "U", "Domaine", degreeQuadrature, &uSol));

    std::vector<feSpace*> spaces = {uBord, u};
    std::vector<feSpace*> essentialSpaces = {uBord};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    feBilinearForm *dudt = nullptr, *diff = nullptr, *src = nullptr;
    feCheckReturn(createBilinearForm(dudt, {u}, new feSysElm_TransientVectorMass<2>(&scalarConstant::one)));
    feCheckReturn(createBilinearForm(diff, {u}, new feSysElm_VectorDiffusion<2>(&scalarConstant::one, &scalarConstant::one)));
    feCheckReturn(createBilinearForm( src, {u}, new feSysElm_VectorSource<2>(&uSrc)));
    std::vector<feBilinearForm*> forms = {dudt, diff, src};

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
    #else
      feCheckReturn(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
    #endif

    // system->setDisplayRHSInConsole(true);

    feNorm *errorU = nullptr, *errurUdot;
    feCheckReturn(createNorm(errorU, VECTOR_LINF_ERROR, {u}, &sol, nullptr, &uSol));
    feCheckReturn(createNorm(errurUdot, VECTOR_LINF_UDOT_ERROR, {u}, &sol, nullptr, &duSoldt));
    std::vector<feNorm *> norms = {errorU, errurUdot};

    TimeIntegrator *solver;
    feCheckReturn(createTimeIntegrator(solver, scheme, NLoptions,
      system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));

    // Start BDF2 with manufactured solution
    if(scheme == timeIntegratorScheme::BDF2)
      static_cast<BDF2Integrator*>(solver)->setStartingMethod(BDF2Starter::InitialCondition);

    feCheckReturn(solver->makeSteps(nTimeSteps));

    int cnt = 0;
    for(auto &val : sol.getSolution()) {
      (cnt++ >= numbering.getNbUnknowns()) ? feInfo("Sol = %+-1.16e (ess)", val) : feInfo("Sol = %+-1.16e", val);
    }

    std::vector<std::vector<double>> postProc = solver->getPostProcessingData();
    error = 0.;
    for(int j = 0; j < nTimeSteps+1; ++j) {
      feInfo("Error at step %2d (t = %1.3e) = %1.16e (U) = %1.16e (Udot)", j, postProc[0][j], postProc[1][j], postProc[2][j]);
      error = fmax(error, postProc[1][j]);
    }

    delete errorU;
    delete solver;
    delete system;
    for(feBilinearForm* f : forms)
      delete f;
    for(feSpace *space : spaces)
      delete space;

    return FE_STATUS_OK;
  }
};

//
// Solve dudt - lap(u) + f = 0 for manufactured solution:
//
// u = exp(-t) * (ax^3 + xy^2, bxy)
//
// with weakly enforced Dirichlet boundary conditions.
//
struct diffusionWeakBC {

  static void uSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    const double a = par[0];
    const double b = par[1];
    res[0] = exp(-t) * (a * pow(x, 3) + x * y*y);
    res[1] = exp(-t) * b * x*y;
  }

  static void duSoldt_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    const double a = par[0];
    const double b = par[1];
    res[0] = - (- exp(-t) * (a * pow(x, 3) + x * y*y));
    res[1] = - (- exp(-t) * b * x*y);
  }

  static void uSrc_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    const double a = par[0];
    const double b = par[1];
    const double expmt = exp(-t);
    const double dudt[2] = {- expmt * (a * pow(x, 3) + x * y*y),
                            - expmt * b * x*y};
    const double lapu[2] = {expmt * (6.*a*x + 2.*x), 0.};                        
    res[0] = - (dudt[0] - lapu[0]);
    res[1] = - (dudt[1] - lapu[1]);
  }

  static feStatus solve(const std::string &meshFile,
                 const timeIntegratorScheme scheme,
                 const double t0,
                 const double t1,
                 const int nTimeSteps,
                 const int order,
                 const int degreeQuadrature,
                 const feNLSolverOptions &NLoptions,
                 double &error)
  {
    double a = 0.1234;
    double b = 1.2345;
    feVectorFunction uSol(uSol_f, {a, b});
    feVectorFunction duSoldt(duSoldt_f, {a, b});
    feVectorFunction uSrc(uSrc_f, {a, b});

    feMesh2DP1 mesh(meshFile);

    // l : Lagrange multiplier to couple de PDE with the weak BC
    // vBord : auxiliary field associated to u on the boundary
    feSpace *u, *uBord, *l, *vBord;
    feCheckReturn(createFiniteElementSpace(    u, &mesh, elementType::VECTOR_LAGRANGE, order, "U", "Domaine", degreeQuadrature, &uSol));    
    feCheckReturn(createFiniteElementSpace(uBord, &mesh, elementType::VECTOR_LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol));    
    feCheckReturn(createFiniteElementSpace(    l, &mesh, elementType::VECTOR_LAGRANGE, order, "L",    "Bord", degreeQuadrature, &vectorConstant::zero));
    feCheckReturn(createFiniteElementSpace(vBord, &mesh, elementType::VECTOR_LAGRANGE, order, "U_aux","Bord", degreeQuadrature, &uSol));

    // There is no essential space here
    std::vector<feSpace*> spaces = {u, uBord, l, vBord}, essentialSpaces = {};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    std::vector<feBilinearForm*> forms = {};

    // Equation 1: PDE for u with Lagrange multiplier l
    feBilinearForm *dudt = nullptr, *diff = nullptr, *src = nullptr, *massLU = nullptr;
    feCheckReturn(createBilinearForm(  dudt, {u}, new feSysElm_TransientVectorMass<2>(&scalarConstant::one)));
    feCheckReturn(createBilinearForm(  diff, {u}, new feSysElm_VectorDiffusion<2>(&scalarConstant::one, &scalarConstant::one)));
    feCheckReturn(createBilinearForm(   src, {u}, new feSysElm_VectorSource<2>(&uSrc)));
    feCheckReturn(createBilinearForm(massLU, {uBord, l}, new feSysElm_MixedVectorMass<2>(&scalarConstant::one)));
    forms.push_back(dudt);
    forms.push_back(diff);
    forms.push_back(src);
    forms.push_back(massLU);

    // Equation 2: Coupling uBord - vBord (equation for l)
    feBilinearForm *massU = nullptr, *massV = nullptr;
    feCheckReturn(createBilinearForm(massU, {l, uBord}, new feSysElm_MixedVectorMass<2>(&scalarConstant::one)));
    feCheckReturn(createBilinearForm(massV, {l, vBord}, new feSysElm_MixedVectorMass<2>(&scalarConstant::minusOne)));
    forms.push_back(massU);
    forms.push_back(massV);

    // Equation 3: Coupling vBord - d/dt of BC (equation for vBord)
    feBilinearForm *dvdt = nullptr, *srcv = nullptr;
    feCheckReturn(createBilinearForm(dvdt, {vBord}, new feSysElm_TransientVectorMass<2>(&scalarConstant::one)));
    feCheckReturn(createBilinearForm(srcv, {vBord}, new feSysElm_VectorSource<2>(&duSoldt)));
    forms.push_back(dvdt);
    forms.push_back(srcv);

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
    #else
      feCheckReturn(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
    #endif

    feNorm *errorU = nullptr;
    feCheckReturn(createNorm(errorU, VECTOR_LINF_ERROR, {u}, &sol, nullptr, &uSol));
    std::vector<feNorm *> norms = {errorU};

    TimeIntegrator *solver;
    feCheckReturn(createTimeIntegrator(solver, scheme, NLoptions,
      system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));
    feCheckReturn(solver->makeSteps(nTimeSteps));

    std::vector<std::vector<double>> postProc = solver->getPostProcessingData();
    error = 0.;
    for(int j = 0; j < nTimeSteps+1; ++j) {
      feInfo("Error at step %2d (t = %1.3e) = %1.6e", j, postProc[0][j], postProc[1][j]);
      error = fmax(error, postProc[1][j]);
    }

    delete errorU;
    delete solver;
    delete system;
    for(feBilinearForm* f : forms)
      delete f;
    for(feSpace *space : spaces)
      delete space;

    return FE_STATUS_OK;
  }
};

struct navierStokesMMS {

  // A "separated" function f(t) such that u(x,t) = f(t) * w(x)
  static double FT(const double t) { return t*t*t*t*t; }
  static double DFDT(const double t) { return 5.*t*t*t*t; }

  // static double FT(const double t) { return exp(-t); }
  // static double DFDT(const double t) { return -exp(-t); }

  // static double FT(const double t) { return cos(t); }
  // static double DFDT(const double t) { return -sin(t); }

  // static double FT(const double t) { return 1./(1.+t); }
  // static double DFDT(const double t) { return -1./((1.+t)*(1.+t)); }

  static void uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    res[0] = FT(t) * x*y;
    res[1] = FT(t) * -y*y/2.;
  }

  static double pSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    return x*y;
    // const double t = args.t;
    // return sin(t) * x*y;
  }

  static void uSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/, std::vector<double> &res)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;

    const double gradp[2] = {y, x};
    // const double gradp[2] = {sin(t) * y, sin(t) * x};

    const double lap_u[2] = {FT(t) * 0.,
                             FT(t) *-1.};
    const double u[2] = {FT(t)*x*y,
                         FT(t)*(-y*y/2.)};

    const double dudt[2] = {DFDT(t) * x*y,
                            DFDT(t) * (-y*y/2.)};

    // Convention: defined as gradu_ij = du_j/dx_i to agree with directional derivatives
    const double gradu[2][2] = {{FT(t)*y, FT(t)*0.},
                                {FT(t)*x, FT(t)*(-y)}};

    const double uDotGradu[2] = {(u[0]*gradu[0][0] + u[1]*gradu[1][0]),
                                 (u[0]*gradu[0][1] + u[1]*gradu[1][1])};
                                 
    res[0] = - (dudt[0] + uDotGradu[0] + gradp[0] - lap_u[0]);
    res[1] = - (dudt[1] + uDotGradu[1] + gradp[1] - lap_u[1]);
  }

  static feStatus solve(const std::string &meshFile,
                 const timeIntegratorScheme scheme,
                 const double t0,
                 const double t1,
                 const int nTimeSteps,
                 const int order,
                 const int degreeQuadrature,
                 const feNLSolverOptions &NLoptions,
                 const bool divergenceFormulation,
                 int &numInteriorElements,
                 double &errorU,
                 double &errorP)
  {
    feVectorFunction uSol(uSol_f);
    feVectorFunction uSrc(uSrc_f);
    feFunction pSol(pSol_f);

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();

    feSpace *u = nullptr, *uB = nullptr, *p = nullptr, *pB = nullptr;
    feCheck(createFiniteElementSpace( u, &mesh, elementType::VECTOR_LAGRANGE,   order, "U", "Domaine", degreeQuadrature, &uSol));
    feCheck(createFiniteElementSpace(uB, &mesh, elementType::VECTOR_LAGRANGE,   order, "U", "Bord"   , degreeQuadrature, &uSol));
    feCheck(createFiniteElementSpace( p, &mesh, elementType::LAGRANGE       , order-1, "P", "Domaine", degreeQuadrature, &pSol));
    feCheck(createFiniteElementSpace(pB, &mesh, elementType::LAGRANGE       , order-1, "P", "PointPression", degreeQuadrature, &pSol));
    std::vector<feSpace*> spaces = {u, uB, p, pB};
    std::vector<feSpace*> essentialSpaces = {uB, pB};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

    // Momentum
    feBilinearForm *dudt = nullptr, *convU = nullptr, *source = nullptr;
    feCheck(createBilinearForm( dudt,    {u}, new feSysElm_TransientVectorMass<2>(&scalarConstant::one)));
    feCheck(createBilinearForm(convU,    {u}, new feSysElm_VectorConvectiveAcceleration<2>(&scalarConstant::one)));
    feCheck(createBilinearForm(source,   {u}, new feSysElm_VectorSource<2>(&uSrc)));
    std::vector<feBilinearForm*> forms = {dudt, convU, source};
    // std::vector<feBilinearForm*> forms = {convU, source};

    feBilinearForm *divSigma = nullptr, *diffU = nullptr, *gradP = nullptr;
    if(divergenceFormulation) {
      UNUSED(diffU, gradP);
      feCheck(createBilinearForm(divSigma, {u, p}, new feSysElm_DivergenceNewtonianStress<2>(&scalarConstant::minusOne, &scalarConstant::one)));
      forms.push_back(divSigma);
    } else {
      // Laplacian formulation
      UNUSED(divSigma);
      feCheck(createBilinearForm(gradP, {u, p}, new feSysElm_MixedGradient<2>(&scalarConstant::one)));
      feCheck(createBilinearForm(diffU,    {u}, new feSysElm_VectorDiffusion<2>(&scalarConstant::one, &scalarConstant::one)));
      forms.push_back(diffU);
      forms.push_back(gradP);
    }

    // Continuity
    feBilinearForm *divU = nullptr;
    feCheck(createBilinearForm(divU, {p, u}, new feSysElm_MixedDivergence<2>(&scalarConstant::one)));
    forms.push_back(divU);

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      feCheck(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheck(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
    #else
      feCheck(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
    #endif

    feNorm *errorUinf, *errorPinf;
    feCheckReturn(createNorm(errorUinf, VECTOR_LINF_ERROR, {u}, &sol, nullptr, &uSol));
    feCheckReturn(createNorm(errorPinf,        LINF_ERROR, {p}, &sol, &pSol));
    std::vector<feNorm *> norms = {errorUinf, errorPinf};

    TimeIntegrator *solver;
    feCheck(createTimeIntegrator(solver, scheme, NLoptions, system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));
    feCheck(solver->makeSteps(nTimeSteps));

    std::vector<std::vector<double>> postProc = solver->getPostProcessingData();
    errorU = 0.;
    errorP = 0.;
    for(int j = 0; j < nTimeSteps+1; ++j) {
      feInfo("Error at step %2d (t = %1.3e) = %1.6e (U) = %1.6e (P)", j, postProc[0][j], postProc[1][j], postProc[2][j]);
      errorU = fmax(errorU, postProc[1][j]);
      errorP = fmax(errorP, postProc[2][j]);
    }

    delete errorUinf;
    delete errorPinf;
    delete solver;
    delete system;
    for(feBilinearForm* f : forms)
      delete f;
    for(feSpace *s : spaces)
      delete s;

    return FE_STATUS_OK;
  }

  static int meshConvergence(std::stringstream &resultBuffer,
                             const timeIntegratorScheme scheme,
                             const std::string &schemeStr,
                             const int order, 
                             const int numMeshes,
                             const int degreeQuadrature)
  {
    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};

    std::vector<int> nElm(numMeshes);
    std::vector<int> timeSteps(numMeshes);
    std::vector<double> errU(numMeshes, 0.);
    std::vector<double> errP(numMeshes, 0.);

    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 10;

    bool divergenceFormulation = true;
    for(int i = 0; i < numMeshes; ++i, nTimeSteps *= 2)
    {
      timeSteps[i] = nTimeSteps;
      // std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
      // std::string meshFile = "../../../data/square1.msh";
      std::string meshFile = "../../../data/stokes1.msh";
      solve(meshFile, scheme, t0, t1, nTimeSteps, order, degreeQuadrature, NLoptions, divergenceFormulation, nElm[i], errU[i], errP[i]);
    }
    resultBuffer << "Taylor-Hood elements P" << order << "-P" << order-1 << std::endl; 
    resultBuffer << "Navier-Stokes MMS - " << schemeStr << " Divergence formulation - Error on u" << std::endl;
    // computeAndPrintConvergence(2, numMeshes, errU, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
    computeAndPrintConvergence(1, numMeshes, errU, timeSteps, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
    resultBuffer << "Navier-Stokes MMS - " << schemeStr << " Divergence formulation - Error on p" << std::endl;
    // computeAndPrintConvergence(2, numMeshes, errP, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
    computeAndPrintConvergence(1, numMeshes, errP, timeSteps, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

    // nTimeSteps = 10;
    // divergenceFormulation = false;
    // for(int i = 0; i < numMeshes; ++i, nTimeSteps *= 2)
    // {
    //   std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
    //   solve(meshFile, timeIntegratorScheme::BDF1, t0, t1, nTimeSteps, order, degreeQuadrature, NLoptions, divergenceFormulation, nElm[i], errU[i], errP[i]);
    // }
    // resultBuffer << "Taylor-Hood elements P" << order << "-P" << order-1 << std::endl; 
    // resultBuffer << "Navier-Stokes MMS - Laplacian formulation - Error on u" << std::endl;
    // computeAndPrintConvergence(2, numMeshes, errU, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
    // resultBuffer << "Navier-Stokes MMS - Laplacian formulation - Error on p" << std::endl;
    // computeAndPrintConvergence(2, numMeshes, errP, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

    return 0;
  }
};

//
// Convergence study - Refine only the time step
//
template<typename problemType>
feStatus timeStepConvergence(std::stringstream &resultBuffer,
                             const timeIntegratorScheme scheme,
                             const double t0,
                             const double t1, 
                             const int initial_nTimeSteps,
                             const int order,
                             const int N,
                             const int degreeQuadrature,
                             const std::string &message)
{
  feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};
  {
    std::vector<int> timeSteps(N);
    std::vector<double> err(N, 0.);
    int nTimeSteps = initial_nTimeSteps;

    for(int i = 0; i < N; ++i, nTimeSteps *= 2)
    {
      timeSteps[i] = nTimeSteps;
      // std::string meshFile = "../../../data/square1.msh";
      std::string meshFile = "../../../data/dummy.msh";
      feCheckReturn(problemType::solve(meshFile, scheme, t0, t1, nTimeSteps, order, degreeQuadrature, NLoptions, err[i]));
    }
    // resultBuffer << "Unsteady diffusion - Error on u - Lagrange elements P" << order << std::endl;
    resultBuffer << message << std::endl;
    computeAndPrintConvergence(1, N, err, timeSteps, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  }
  return FE_STATUS_OK;
}

// //
// // Convergence study - Refine mesh and time step
// //
// template<typename problemType>
// feStatus meshAndTimeStepConvergence(std::stringstream &resultBuffer,
//                                     const timeIntegratorScheme scheme,
//                                     const double t0,
//                                     const double t1, 
//                                     const int initial_nTimeSteps,
//                                     const int order,
//                                     const int N,
//                                     const int degreeQuadrature)
// {
//   feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};
//   {
//     std::vector<int> nElm(N);
//     std::vector<int> timeSteps(N);
//     std::vector<double> err(N, 0.);
//     int nTimeSteps = initial_nTimeSteps;

//     for(int i = 0; i < N; ++i, nTimeSteps *= 2)
//     {
//       timeSteps[i] = nTimeSteps;
//       // std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
//       std::string meshFile = "../../../data/square1.msh";
//       feCheckReturn(problemType::solve(meshFile, scheme, t0, t1, nTimeSteps, order, degreeQuadrature, NLoptions, nElm[i], err[i]));
//     }
//     resultBuffer << "Unsteady vector diffusion - Error on u - Lagrange elements P" << order << std::endl;
//     // computeAndPrintConvergence(2, N, err, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
//     computeAndPrintConvergence(1, N, err, timeSteps, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
//   }
//   return FE_STATUS_OK;
// }

TEST(VectorFEUnsteady, LinearODE)
{
  initialize(my_argc, my_argv);
  // setVerbose(0);
  std::string testRoot = "../../../tests/withLinearSolver/vectorFEunsteady_linearODE";
  std::stringstream resultBuffer;
  int degreeQuadrature = 2;
  double t0 = 0.;
  double t1 = 1.;
  int initial_nTimeSteps = 10;

  for(size_t i = 0; i < schemes.size(); ++i)
  {
    std::string message = "Linear ODE - Error on vector-valued solution u - Time integrator: " + schemesStr[i];
    ASSERT_TRUE(timeStepConvergence<linearODE>(resultBuffer, schemes[i], t0, t1, initial_nTimeSteps, 1, 5, degreeQuadrature, message) == FE_STATUS_OK);
  }

  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

// TEST(VectorFEUnsteady, Diffusion)
// {
//   initialize(my_argc, my_argv);
//   // setVerbose(0);
//   std::string testRoot = "../../../tests/withLinearSolver/vectorFEunsteady_diffusion";
//   std::stringstream resultBuffer;
//   int degreeQuadrature = 4;
//   double t0 = 0.;
//   double t1 = 1.;
//   int initial_nTimeSteps = 1;

//   // for(size_t i = 0; i < schemes.size(); ++i)
//   for(size_t i = 2; i < 3; ++i)
//   // for(size_t i = 3; i < 4; ++i)
//   {
//     std::string message = "Diffusion PDE w/ strong BC - Error on vector-valued solution u - Time integrator: " + schemesStr[i];
//     ASSERT_TRUE(timeStepConvergence<diffusion>(resultBuffer, schemes[i], t0, t1, initial_nTimeSteps, 1, 1, degreeQuadrature, message) == FE_STATUS_OK);
//   }

//   EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
//   finalize();
// }

// TEST(VectorFEUnsteady, DiffusionWeakBC)
// {
//   initialize(my_argc, my_argv);
//   // setVerbose(0);
//   std::string testRoot = "../../../tests/withLinearSolver/vectorFEunsteady_diffusionWeakBC";
//   std::stringstream resultBuffer;
//   int degreeQuadrature = 8;
//   double t0 = 0.;
//   double t1 = 1.;
//   int initial_nTimeSteps = 10;

//   for(size_t i = 0; i < schemes.size(); ++i)
//   {
//     std::string message = "Diffusion PDE w/ weak BC - Error on vector-valued solution u - Time integrator: " + schemesStr[i];
//     ASSERT_TRUE(timeStepConvergence<diffusionWeakBC>(resultBuffer, schemes[i], t0, t1, initial_nTimeSteps, 3, 4, degreeQuadrature, message) == FE_STATUS_OK);
//   }

//   EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
//   finalize();
// }

// TEST(VectorFEUnsteady, NSMMS)
// {
//   initialize(my_argc, my_argv);
//   // setVerbose(0);
//   std::string testRoot = "../../../tests/withLinearSolver/vectorFEunsteady_nsmms";
//   std::stringstream resultBuffer;
//   int degreeQuadrature = 8;
//   // double t0 = 0.;
//   // double t1 = 1.;
//   // int initial_nTimeSteps = 10;

//   for(size_t i = 0; i < schemes.size(); ++i)
//   {
//     // std::string message = " - Error on vector-valued solution u - Time integrator: " + schemesStr[i];
//     // ASSERT_TRUE(timeStepConvergence<diffusion>(resultBuffer, schemes[i], t0, t1, initial_nTimeSteps, 3, 4, degreeQuadrature, message) == FE_STATUS_OK);
//     navierStokesMMS::meshConvergence(resultBuffer, schemes[i], schemesStr[i], 3, 6, degreeQuadrature);
//   }

//   EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
//   finalize();
// }

//
// Custom main to capture argc and argv
//
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  my_argc = argc;
  my_argv = argv;
  return RUN_ALL_TESTS();
}