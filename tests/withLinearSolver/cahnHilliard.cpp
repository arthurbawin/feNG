
#include "feAPI.h"
#include "feParameters.h"
#include "CHNS_Solver.h"
#include "../tests.h"

#include <gtest/gtest.h>

/*
  Tests the implementation of the Cahn-Hilliard and CH-Navier-Stokes weak forms.
*/

static int my_argc;
static char** my_argv;

///////////////////////////////////////////////////////////////////////////////
// Modify here for another MMS

// #define TEST_TIME_CONVERGENCE

#define CONSTANT_MOBILITY(M, phi) M
#define DCONSTANT_MOBILITY_DPHI(M, phi) 0. 

// Time dependency.
// Note: the phi^3 term in the potential equation is of order FT*FT*FT,
// so a BDF2 cannot capture this exactly even for linear FT(t).
//
// #define FT(t) 1.
// #define DFDT(t) 0.

// #define FT(t) t
// #define DFDT(t) 1.

// #define FT(t) t*t
// #define DFDT(t) 2.*t

#define FT(t) t*t*t*t
#define DFDT(t) 4.*t*t*t

// #define FT(t) exp(-t)
// #define DFDT(t) -exp(-t)

// #define FT(t) sin(t)
// #define DFDT(t) cos(t)

#if !defined(TEST_TIME_CONVERGENCE)
  // MMS 1 for spatial convergence
  //
  // #define CU 1.
  // #define CV 1.
  // #define U(x,y,t) CU * FT(t) * x*x*x*x * y*y*y*y
  // #define V(x,y,t) CV * FT(t) * -4./5. * x*x*x * y*y*y*y*y

  // #define DUDT(x,y,t) CU * DFDT(t) * x*x*x*x * y*y*y*y
  // #define DVDT(x,y,t) CV * DFDT(t) * -4./5. * x*x*x * y*y*y*y*y

  // #define DUDX(x,y,t) CU * FT(t) *      4. *   x*x*x *   y*y*y*y
  // #define DUDY(x,y,t) CU * FT(t) *      4. * x*x*x*x *     y*y*y
  // #define DVDX(x,y,t) CV * FT(t) * -12./5. *     x*x * y*y*y*y*y
  // #define DVDY(x,y,t) CV * FT(t) *     -4. *   x*x*x *   y*y*y*y

  // #define UXX(x,y,t) CU * FT(t) *  12.    *          x*x * y*y*y*y
  // #define UXY(x,y,t) CU * FT(t) *  16.    *        x*x*x * y*y*y
  // #define UYY(x,y,t) CU * FT(t) *  12.    *      x*x*x*x * y*y
  // #define VXX(x,y,t) CV * FT(t) * -24./5. *            x * y*y*y*y*y
  // #define VXY(x,y,t) CV * FT(t) * -12.    *          x*x * y*y*y*y
  // #define VYY(x,y,t) CV * FT(t) * -16.    *        x*x*x * y*y*y

  // // #define P(x,y,t) x*x*y*y - 1./9. // Zero-average on [0,1]^2
  // #define P(x,y,t) x*x*y*y // Zero-average on [0,1]^2
  // #define DPDX(x,y,t) 2.*x*y*y
  // #define DPDY(x,y,t) 2.*x*x*y
  // #define DPDXX(x,y,t) 2.*y*y
  // #define DPDYY(x,y,t) 2.*x*x

  #define U(x,y,t) FT(t) * sin(M_PI*x) * sin(M_PI*y)
  #define V(x,y,t) FT(t) * cos(M_PI*x) * cos(M_PI*y)

  #define DUDT(x,y,t) DFDT(t) * sin(M_PI*x) * sin(M_PI*y)
  #define DVDT(x,y,t) DFDT(t) * cos(M_PI*x) * cos(M_PI*y)

  #define DUDX(x,y,t)   FT(t) * M_PI * cos(M_PI*x) * sin(M_PI*y)
  #define DUDY(x,y,t)   FT(t) * M_PI * sin(M_PI*x) * cos(M_PI*y)
  #define DVDX(x,y,t) - FT(t) * M_PI * sin(M_PI*x) * cos(M_PI*y)
  #define DVDY(x,y,t) - FT(t) * M_PI * cos(M_PI*x) * sin(M_PI*y)

  #define UXX(x,y,t) - FT(t) * M_PI * M_PI * sin(M_PI*x) * sin(M_PI*y)
  #define UXY(x,y,t)   FT(t) * M_PI * M_PI * cos(M_PI*x) * cos(M_PI*y)
  #define UYY(x,y,t) - FT(t) * M_PI * M_PI * sin(M_PI*x) * sin(M_PI*y)
  #define VXX(x,y,t) - FT(t) * M_PI * M_PI * cos(M_PI*x) * cos(M_PI*y)
  #define VXY(x,y,t)   FT(t) * M_PI * M_PI * sin(M_PI*x) * sin(M_PI*y)
  #define VYY(x,y,t) - FT(t) * M_PI * M_PI * cos(M_PI*x) * cos(M_PI*y)

  #define P(x,y,t)                     sin(M_PI*x) * cos(M_PI*y)
  #define DPDX(x,y,t)           M_PI * cos(M_PI*x) * cos(M_PI*y)
  #define DPDY(x,y,t)         - M_PI * sin(M_PI*x) * sin(M_PI*y)
  #define DPDXX(x,y,t) - M_PI * M_PI * sin(M_PI*x) * cos(M_PI*y)
  #define DPDYY(x,y,t) - M_PI * M_PI * sin(M_PI*x) * cos(M_PI*y)

  #define PHI(x,y,t)                     FT(t) * sin(M_PI*x) * sin(M_PI*y)
  #define DPHIDT(x,y,t)                DFDT(t) * sin(M_PI*x) * sin(M_PI*y)
  #define DPHIDX(x,y,t)           M_PI * FT(t) * cos(M_PI*x) * sin(M_PI*y)
  #define DPHIDY(x,y,t)           M_PI * FT(t) * sin(M_PI*x) * cos(M_PI*y)
  #define DPHIDXX(x,y,t) - M_PI * M_PI * FT(t) * sin(M_PI*x) * sin(M_PI*y)
  #define DPHIDYY(x,y,t) - M_PI * M_PI * FT(t) * sin(M_PI*x) * sin(M_PI*y)

  #define MU(x,y,t)                     sin(M_PI*x) * sin(M_PI*y)
  #define DMUDX(x,y,t)           M_PI * cos(M_PI*x) * sin(M_PI*y)
  #define DMUDY(x,y,t)           M_PI * sin(M_PI*x) * cos(M_PI*y)
  #define DMUDXX(x,y,t) - M_PI * M_PI * sin(M_PI*x) * sin(M_PI*y)
  #define DMUDYY(x,y,t) - M_PI * M_PI * sin(M_PI*x) * sin(M_PI*y)

#else 
  // MMS 2 for time convergence
  // Spatial component captured by P2-P1.
  #define U(x,y,t) FT(t) * (x*x - 2.*x + 3.*y*y + y - x*y + 1)
  #define V(x,y,t) FT(t) * (3.*x*x + y*y/2. + 2.*y - 1. - 2.*x*y)

  #define DUDT(x,y,t) DFDT(t) * (x*x - 2.*x + 3.*y*y + y - x*y + 1)
  #define DVDT(x,y,t) DFDT(t) * (3.*x*x + y*y/2. + 2.*y - 1. - 2.*x*y)

  #define DUDX(x,y,t) FT(t) * (2.*x - 2. - y)
  #define DUDY(x,y,t) FT(t) * (6.*y + 1. - x)
  #define DVDX(x,y,t) FT(t) * (6.*x - 2.*y)
  #define DVDY(x,y,t) FT(t) * (y + 2. - 2.*x)

  #define UXX(x,y,t) FT(t) * 2.
  #define UXY(x,y,t) FT(t) * (-1.)
  #define UYY(x,y,t) FT(t) * 6.
  #define VXX(x,y,t) FT(t) * 6.
  #define VXY(x,y,t) FT(t) * (-2.)
  #define VYY(x,y,t) FT(t) * 1.

  #define P(x,y,t) 10. * (x + y - 1.)
  #define DPDX(x,y,t) 10.
  #define DPDY(x,y,t) 10.
  #define DPDXX(x,y,t) 0.
  #define DPDYY(x,y,t) 0.

  #define PHI(x,y,t)      FT(t) * (x + y - 1.)
  #define DPHIDT(x,y,t) DFDT(t) * (x + y - 1.)
  #define DPHIDX(x,y,t)   FT(t) * (1.)
  #define DPHIDY(x,y,t)   FT(t) * (1.)
  #define DPHIDXX(x,y,t)  0.
  #define DPHIDYY(x,y,t)  0.

  #define MU(x,y,t)      (x + y - 1.)
  #define DMUDX(x,y,t)   (1.)
  #define DMUDY(x,y,t)   (1.)
  #define DMUDXX(x,y,t)  0.
  #define DMUDYY(x,y,t)  0.
#endif
///////////////////////////////////////////////////////////////////////////////
// These should not be changed (-:
#define UVEC(x,y,t) {U(x,y,t), V(x,y,t)}
#define DUDTVEC(x,y,t) {DUDT(x,y,t), DVDT(x,y,t)}

#define DIVU(x,y,t) DUDX(x,y,t) + DVDY(x,y,t)

// Convention: defined as gradu_ij = du_j/dx_i to agree with directional derivatives
#define GRADU(x,y,t) {{DUDX(x,y,t), DVDX(x,y,t)},    \
                      {DUDY(x,y,t), DVDY(x,y,t)}}    \

#define SYM_GRADU(x,y,t) {{2. * DUDX(x,y,t), DUDY(x,y,t) + DVDX(x,y,t)},    \
                          {DUDY(x,y,t) + DVDX(x,y,t), 2. * DVDY(x,y,t)}}    \

// Divergence of (grad(u) + grad(u)^T) / 2
#define DIV_RATE_STRAIN(x,y,t) {UXX(x,y,t) + 0.5 * (UYY(x,y,t) + VXY(x,y,t)), \
                                VYY(x,y,t) + 0.5 * (UXY(x,y,t) + VXX(x,y,t))}

#define GRAD_DIV_U(x,y,t) {UXX(x,y,t) + VXY(x,y,t), UXY(x,y,t) + VYY(x,y,t)}

#define GRADP(x,y,t) {DPDX(x,y,t), DPDY(x,y,t)}
#define LAP_P(x,y,t) DPDXX(x,y,t) + DPDYY(x,y,t)

#define GRADPHI(x,y,t) {DPHIDX(x,y,t), DPHIDY(x,y,t)}
#define LAP_PHI(x,y,t) DPHIDXX(x,y,t) + DPHIDYY(x,y,t)

#define GRADMU(x,y,t) {DMUDX(x,y,t), DMUDY(x,y,t)}
#define LAP_MU(x,y,t) DMUDXX(x,y,t) + DMUDYY(x,y,t)
///////////////////////////////////////////////////////////////////////////////

#define DECLARE_VARIABLES(args, par)   \
  const double x   = args.pos[0];      \
  const double y   = args.pos[1];      \
  const double t   = args.t;           \
  const double phi = args.u;           \
                                       \
  UNUSED(x, y, t, phi, par);           \

void uSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
{
  DECLARE_VARIABLES(args, par);
  res[0] = U(x,y,t);
  res[1] = V(x,y,t);
}

void grad_uSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
{
  DECLARE_VARIABLES(args, par);
  const double gradu[2][2] = GRADU(x,y,t);
  res[0] = gradu[0][0];
  res[1] = gradu[0][1];
  res[2] = gradu[1][0];
  res[3] = gradu[1][1];
}

double pSol_f(const feFunctionArguments &args, const std::vector<double> &par)
{
  DECLARE_VARIABLES(args, par);
  return P(x,y,t);
}

void grad_pSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
{
  DECLARE_VARIABLES(args, par);
  const double gradp[2] = GRADP(x,y,t);
  res[0] = gradp[0];
  res[1] = gradp[1];
}

double phiSol_f(const feFunctionArguments &args, const std::vector<double> &par)
{
  DECLARE_VARIABLES(args, par);
  return PHI(x,y,t);
}

void grad_phiSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
{
  DECLARE_VARIABLES(args, par);
  const double gradphi[2] = GRADPHI(x,y,t);
  res[0] = gradphi[0];
  res[1] = gradphi[1];
}

double muSol_f(const feFunctionArguments &args, const std::vector<double> &par)
{
  DECLARE_VARIABLES(args, par);
  return MU(x,y,t);
}

void grad_muSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
{
  DECLARE_VARIABLES(args, par);
  const double gradmu[2] = GRADMU(x,y,t);
  res[0] = gradmu[0];
  res[1] = gradmu[1];
}

//
// Unsteady Cahn-Hilliard system without Navier-Stokes
// Lagrange elements of equal order for (phi,mu) (phase marker and chemical potential)
//
namespace cahnHilliard {

  double phiSol_f(const feFunctionArguments &args, const std::vector<double> & /* par */)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    return exp(-t) * sin(M_PI*x) * sin(M_PI*y);
  }

  double muSol_f(const feFunctionArguments &args, const std::vector<double> & /* par */)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    return sin(M_PI*x) * sin(M_PI*y);
  }

  double phiSource_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    const double sx = sin(M_PI*x);
    const double cx = cos(M_PI*x);
    const double sy = sin(M_PI*y);
    const double cy = cos(M_PI*y);
    const double dphidt = - exp(-t) * sx*sy;
    const double lapMu  =              -2.*M_PI*M_PI*sx*sy;
    const double uDotGradPhi = exp(-t) * M_PI*(cx*sy + sx*cy);
    return - dphidt - uDotGradPhi + lapMu;
  }

  double muSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    const double epsilon = par[0];
    const double sx = sin(M_PI*x);
    const double sy = sin(M_PI*y);
    const double emt = exp(-t);
    const double mu = sx*sy;
    const double phi     = emt * sx*sy;
    const double phiCube = emt*emt*emt*sx*sx*sx*sy*sy*sy;
    const double lapPhi  = emt * (-2.*M_PI*M_PI*sx*sy);
    return - mu + phiCube - phi - epsilon*epsilon*lapPhi;
  }

  int solve(const std::string &meshFile,
            const int order,
            const int degreeQuadrature,
            const feNLSolverOptions &NLoptions,
            int &numInteriorElements,
            double &L2ErrorPhi,
            double &L2ErrorMu)
  {
    const double epsilon = 1.;
    feFunction phiSol(phiSol_f);
    feFunction  muSol( muSol_f);
    feFunction phiSource(phiSource_f);
    feFunction  muSource( muSource_f, {epsilon});

    feConstantVectorFunction cst_velocity({1., 1.});
    feConstantFunction minusEpsInterfaceSquared(-epsilon*epsilon);
    feConstantFunction mobility(1.);

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();

    // Phase marker
    feSpace *phi = nullptr, *phiBord = nullptr, *u = nullptr;
    // Chemical potential
    feSpace *mu = nullptr,  *muBord = nullptr;

    feCheck(createFiniteElementSpace(    phi, &mesh, elementType::LAGRANGE, order, "Phi", "Domaine", degreeQuadrature, &phiSol));
    feCheck(createFiniteElementSpace(phiBord, &mesh, elementType::LAGRANGE, order, "Phi",    "Bord", degreeQuadrature, &phiSol));
    feCheck(createFiniteElementSpace(      u, &mesh, elementType::VECTOR_LAGRANGE, 1, "U",   "Domaine", degreeQuadrature, &cst_velocity));
    feCheck(createFiniteElementSpace(     mu, &mesh, elementType::LAGRANGE, order,  "Mu",  "Domaine", degreeQuadrature, &muSol));
    feCheck(createFiniteElementSpace( muBord, &mesh, elementType::LAGRANGE, order,  "Mu",     "Bord", degreeQuadrature, &muSol));

    std::vector<feSpace*> spaces = {phi, mu, phiBord, muBord, u};
    std::vector<feSpace*> essentialSpaces = {phiBord, muBord, u};
    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

    // Tracer convection with imposed velocity field u = (1,1)
    feBilinearForm *dphidt = nullptr, *convPhi = nullptr, *diffMu = nullptr, *phiS = nullptr;
    feCheck(createBilinearForm( dphidt,     {phi}, new feSysElm_TransientMass(&scalarConstant::one)));
    feCheck(createBilinearForm(convPhi,  {phi, u}, new feSysElm_TracerConvection<2>(&scalarConstant::one)));
    feCheck(createBilinearForm( diffMu, {phi, mu}, new feSysElm_MixedGradGrad<2>(&mobility)));
    feCheck(createBilinearForm(   phiS,     {phi}, new feSysElm_Source(&phiSource)));

    // Chemical potential
    feBilinearForm *muMass = nullptr, *phiCube = nullptr, *phiMass = nullptr, *diffPhi = nullptr, *muS = nullptr;
    feCheck(createBilinearForm( muMass,      {mu}, new feSysElm_Mass(&scalarConstant::one)));
    feCheck(createBilinearForm(phiCube, {mu, phi}, new feSysElm_MixedMassPower(&scalarConstant::minusOne, 3)));
    feCheck(createBilinearForm(phiMass, {mu, phi}, new feSysElm_MixedMass(&scalarConstant::one)));
    feCheck(createBilinearForm(diffPhi, {mu, phi}, new feSysElm_MixedGradGrad<2>(&minusEpsInterfaceSquared)));
    feCheck(createBilinearForm(    muS,      {mu}, new feSysElm_Source(&muSource)));

    std::vector<feBilinearForm*> forms = {dphidt, convPhi, diffMu, phiS, muMass, phiCube, phiMass, diffPhi, muS};

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      feCheck(createLinearSystem(system, MKLPARDISO, forms, &numbering));
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheck(createLinearSystem(system, PETSC_MUMPS, forms, &numbering));
    #else
      feCheck(createLinearSystem(system, PETSC, forms, &numbering));
    #endif

    std::vector<feNorm *> norms = {};

    double t0 = 0.;
    double t1 = 0.1;
    int nTimeSteps = 20;
    TimeIntegrator *solver;
    feCheck(createTimeIntegrator(solver, timeIntegratorScheme::BDF2, NLoptions,
      system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));

    // Solve
    feCheck(solver->makeSteps(nTimeSteps));

    feNorm *errorPhi_L2 = nullptr, *errorMu_L2 = nullptr;
    feCheck(createNorm(errorPhi_L2, L2_ERROR, {phi}, &sol, &phiSol));
    feCheck(createNorm(errorMu_L2 , L2_ERROR, {mu}, &sol, &muSol));
    L2ErrorPhi = errorPhi_L2->compute();
    L2ErrorMu  = errorMu_L2->compute();

    delete errorPhi_L2;
    delete errorMu_L2;
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
    std::vector<double> errPhi(numMeshes, 0.);
    std::vector<double> errMu(numMeshes, 0.);

    for(int i = 0; i < numMeshes; ++i)
    {
      std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
      cahnHilliard::solve(meshFile, order, degreeQuadrature, NLoptions, nElm[i], errPhi[i], errMu[i]);
    }
    resultBuffer << "(phi,mu) : Equal order elements P" << order << "-P" << order << std::endl; 
    resultBuffer << "Cahn-Hilliard MMS - Error on phase marker phi" << std::endl;
    computeAndPrintConvergence(2, numMeshes, errPhi, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
    resultBuffer << "Cahn-Hilliard MMS - Error on chemical potential mu" << std::endl;
    computeAndPrintConvergence(2, numMeshes, errMu, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

    return 0;
  }
}

//
// Parameters for Abels formulation of CHNS
//
#define DECLARE_PARAMETERS(args, par)   \
                                        \
  const double rhoA       = par[0];     \
  const double rhoB       = par[1];     \
  const double viscA      = par[2];     \
  const double viscB      = par[3];     \
  const double M          = par[4];     \
  const double gamma      = par[5];     \
  const double epsilon    = par[6];     \
  const double lambda     = par[7];     \
                                        \
  UNUSED(args, rhoA, rhoB,              \
    viscA, viscB, M, gamma,             \
    epsilon, lambda);                   \


struct CHNS_VolumeAveraged_Abels {

  //
  // Source terms
  //
  static void uSource_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);
    DECLARE_PARAMETERS(args, par);

    const double u[2]       = UVEC(x,y,t);
    const double dudt[2]    = DUDTVEC(x,y,t);
    const double phiRef     = PHI(x,y,t);
    const double mu         = MU(x,y,t);

    const double Mphi = CONSTANT_MOBILITY(M, phiRef);
    const double rho  = ( rhoA -  rhoB)/2. * phiRef + ( rhoA +  rhoB)/2.;
    const double visc = (viscA - viscB)/2. * phiRef + (viscA + viscB)/2.;

    const double dviscdphi = (viscA - viscB)/2.;

    const double gradu[2][2]    = GRADU(x,y,t);
    const double symgradu[2][2] = SYM_GRADU(x,y,t);
    const double gradp[2]       = GRADP(x,y,t);
    const double gradphi[2]     = GRADPHI(x,y,t);
    const double gradmu[2]      = GRADMU(x,y,t);

    const double uDotGradu[2]      = {     u[0]*gradu[0][0] +      u[1]*gradu[1][0],
                                           u[0]*gradu[0][1] +      u[1]*gradu[1][1]};
    const double gradMudotGradU[2] = {gradmu[0]*gradu[0][0] + gradmu[1]*gradu[1][0],
                                      gradmu[0]*gradu[0][1] + gradmu[1]*gradu[1][1]};
    const double f[2] = {0., -9.81};

    const double gradPhidotSymGradu[2] = {gradphi[0]*symgradu[0][0] + gradphi[1]*symgradu[1][0],
                                          gradphi[0]*symgradu[0][1] + gradphi[1]*symgradu[1][1]};

    // This is actually div(symGradU/2), should fix it to be consistent
    const double divd_incompressible[2] = DIV_RATE_STRAIN(x,y,t);

    //
    // Main terms
    //
    const double acceleration[2]  = {rho * (dudt[0] + uDotGradu[0] - f[0]),
                                     rho * (dudt[1] + uDotGradu[1] - f[1])};

    const double diffusiveFlux[2] = {(rhoA - rhoB)/2. * Mphi * gradMudotGradU[0],
                                     (rhoA - rhoB)/2. * Mphi * gradMudotGradU[1]};

    const double div_d[2] = {dviscdphi * gradPhidotSymGradu[0] + visc * (2. * divd_incompressible[0]),
                             dviscdphi * gradPhidotSymGradu[1] + visc * (2. * divd_incompressible[1])};

    // res[0] = - (acceleration[0] + diffusiveFlux[0] - div_d[0] + gradp[0] - mu * gradphi[0]);
    // res[1] = - (acceleration[1] + diffusiveFlux[1] - div_d[1] + gradp[1] - mu * gradphi[1]);

    UNUSED(mu);
    res[0] = - (acceleration[0] - diffusiveFlux[0] - div_d[0] + gradp[0] + phiRef * gradmu[0]);
    res[1] = - (acceleration[1] - diffusiveFlux[1] - div_d[1] + gradp[1] + phiRef * gradmu[1]);
  }

  static double pSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    const double div_u = DIVU(x,y,t);
    return - (div_u);
  }

  static double phiSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    DECLARE_PARAMETERS(args, par);

    const double phiRef      = PHI(x,y,t);
    const double dphidt      = DPHIDT(x,y,t);
    const double u[2]        = UVEC(x,y,t);
    const double gradphi[2]  = GRADPHI(x,y,t);
    const double gradMu[2]   = GRADMU(x,y,t);
    const double lapMu       = LAP_MU(x,y,t);

    UNUSED(phiRef);
    const double Mphi   = CONSTANT_MOBILITY(M, phiRef);
    const double dMdphi = DCONSTANT_MOBILITY_DPHI(M,phiRef);

    const double uDotGradPhi = u[0]*gradphi[0] + u[1]*gradphi[1];
    const double gradPhi_dot_gradMu = gradphi[0] * gradMu[0] + gradphi[1] * gradMu[1];

    const double div_mobility_gradMu = Mphi * lapMu + dMdphi * gradPhi_dot_gradMu;

    return - (dphidt + uDotGradPhi - div_mobility_gradMu);
  }

  static double muSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    DECLARE_PARAMETERS(args, par);

    const double mu      = MU(x,y,t);
    const double phiRef  = PHI(x,y,t);
    const double phiCube = PHI(x,y,t) * PHI(x,y,t) * PHI(x,y,t);
    const double lapPhi  = LAP_PHI(x,y,t);

    return - (mu - lambda / (epsilon*epsilon) * (phiCube - phiRef) + lambda * lapPhi);
  }

  // CHNS_VolumeAveraged_Abels
  static feStatus solve(const std::string &meshFile,
                        const double t0,
                        const double t1,
                        const int timeStepFactor,
                        const std::vector<int> &orders,
                        const int degreeQuadrature,
                        int &numInteriorElements,
                        int &numVertices,
                        std::vector<double> &L2_errors,
                        std::vector<double> &H1_errors,
                        int &timesteps)
  {
    const double rhoA    = 1000.;
    const double rhoB    = 1.;
    const double viscA   = 100.;
    const double viscB   = 1.;
    const double M       = 1e-6;
    const double gamma   = 0.01;
    const double epsilon = 1e-3;

    // const double rhoA    = 1.;
    // const double rhoB    = 1.;
    // const double viscA   = 1.;
    // const double viscB   = 1.;
    // const double M       = 1.;
    // const double gamma   = 1.;
    // const double epsilon = 1.;

    const double lambda = 3. / (2.*sqrt(2.)) * epsilon * gamma;

    std::vector<double> params = {rhoA, rhoB, viscA, viscB, M, gamma, epsilon, lambda};

    Parameters::CHNS CHNSparameters;
    CHNSparameters.formulation = Parameters::CHNS::Formulation::abels;
    CHNSparameters.mobilityType = Parameters::CHNS::MobilityType::constant;
    CHNSparameters.fluids[0].density           = rhoA;
    CHNSparameters.fluids[0].dynamic_viscosity = viscA;
    CHNSparameters.fluids[1].density           = rhoB;
    CHNSparameters.fluids[1].dynamic_viscosity = viscB;
    CHNSparameters.gravity        = -9.81;
    CHNSparameters.surfaceTension = gamma;
    CHNSparameters.epsilon        = epsilon;
    CHNSparameters.mobility       = M;
    CHNSparameters.limitPhaseMarker = false;

    feVectorFunction uSol(  uSol_f, params);
    feFunction       pSol(  pSol_f, params);
    feFunction     phiSol(phiSol_f, params);
    feFunction      muSol( muSol_f, params);

    feVectorFunction   grad_uSol(  grad_uSol_f, params);
    feVectorFunction   grad_pSol(  grad_pSol_f, params);
    feVectorFunction grad_phiSol(grad_phiSol_f, params);
    feVectorFunction  grad_muSol( grad_muSol_f, params);

    feVectorFunction   uSource(  uSource_f, params);
    feFunction         pSource(  pSource_f, params);
    feFunction       phiSource(phiSource_f, params);
    feFunction        muSource( muSource_f, params);

    const elementType      lag    = elementType::LAGRANGE;
    const elementType      vlag   = elementType::VECTOR_LAGRANGE;

    VectorFEDescriptor   u("U", "Domaine", vlag, orders[0], degreeQuadrature, &uSol);
    ScalarFEDescriptor   p("P", "Domaine", lag, orders[1], degreeQuadrature, &pSol);
    ScalarFEDescriptor phi("Phi", "Domaine", lag, orders[2], degreeQuadrature, &phiSol);
    ScalarFEDescriptor  mu("Mu", "Domaine", lag, orders[3], degreeQuadrature, &muSol);

    std::vector<FEDescriptor> fieldDescriptors = {u, p, phi, mu};

    VectorFEDescriptor uB("U", "Bord", vlag, orders[0], degreeQuadrature, &uSol);
    ScalarFEDescriptor pP("P", "PointPression", lag, orders[1], degreeQuadrature, &pSol);
    ScalarFEDescriptor phiB("Phi", "Bord", lag, orders[2], degreeQuadrature, &phiSol);
    ScalarFEDescriptor muB("Mu", "Bord", lag, orders[3], degreeQuadrature, &muSol);

    BoundaryConditions::VectorDirichlet uBC(uB);
    BoundaryConditions::ScalarDirichlet pBC(pP), phiBC(phiB), muBC(muB);

    std::vector<const BoundaryConditions::BoundaryCondition*> BC = {&uBC, &pBC, &phiBC, &muBC};

    Parameters::NonLinearSolver newton_param;

    Parameters::TimeIntegration time_param;
    time_param.t_initial  = t0;
    time_param.t_final    = t1;

    #if defined(TEST_TIME_CONVERGENCE)
      time_param.nTimeSteps = 10 * pow(2, timeStepFactor);
    #else 
      time_param.nTimeSteps = 10 * pow(2, timeStepFactor);
    #endif

    timesteps = time_param.nTimeSteps;

    time_param.method =
      Parameters::TimeIntegration::TimeIntegrationMethod::bdf2;
    time_param.bdf2starter =
      // Parameters::TimeIntegration::BDF2Starter::bdf1;
      Parameters::TimeIntegration::BDF2Starter::exactSolution;


    CHNS_Solver *solver;
    feCheck(createCHNS_Solver(solver,
                              CHNSparameters,
                              newton_param,
                              time_param,
                              fieldDescriptors,
                              BC,
                              nullptr,
                              &pSource,
                              &uSource,
                              &phiSource,
                              &muSource));

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();
    numVertices = mesh.getNumVertices();

    std::vector<feSpace *> spaces, essentialSpaces;
    solver->createSpaces(&mesh, spaces, essentialSpaces);

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution   sol(numbering.getNbDOFs(), spaces, essentialSpaces);

    feExporter *exporter = nullptr;
    feExportData exportData = {exporter, 1, "sol"};

    feNorm *uL2, *uH1, *pL2, *pH1, *phiL2, *phiH1, *muL2, *muH1;
    feCheckReturn(createNorm(uL2, VECTOR_L2_ERROR, {spaces[0]}, &sol, nullptr, &uSol));
    feCheckReturn(createNorm(uH1, VECTOR_SEMI_H1_ERROR, {spaces[0]}, &sol, nullptr, &grad_uSol));
    feCheckReturn(createNorm(pL2, L2_ERROR, {spaces[1]}, &sol, &pSol));
    feCheckReturn(createNorm(pH1, SEMI_H1_ERROR, {spaces[1]}, &sol, nullptr, &grad_pSol));
    feCheckReturn(createNorm(phiL2, L2_ERROR, {spaces[2]}, &sol, &phiSol));
    feCheckReturn(createNorm(phiH1, SEMI_H1_ERROR, {spaces[2]}, &sol, nullptr, &grad_phiSol));
    feCheckReturn(createNorm(muL2, L2_ERROR, {spaces[3]}, &sol, &muSol));
    feCheckReturn(createNorm(muH1, SEMI_H1_ERROR, {spaces[3]}, &sol, nullptr, &grad_muSol));

    std::vector<feNorm*> norms = {uL2, pL2, phiL2, muL2, uH1, pH1, phiH1, muH1};
    solver->setNorms(norms);

    feCheck(solver->solve(
      &mesh, &sol, &numbering, spaces, exportData));

    std::vector<std::vector<double>> postProcData = solver->getPostProcessingData();

    for(unsigned int i = 0; i < L2_errors.size(); ++i)
    {
      L2_errors[i] = 0.;
      H1_errors[i] = 0.;
    }

    // First post proc vector is the time
    const double dt = postProcData[0][1] - postProcData[0][0];

    // Compute L1 norm in time of the computed spatial errors
    for(size_t i = 0; i < postProcData[0].size(); ++i)
    {
      L2_errors[0] += dt * postProcData[1][i]; // u 
      L2_errors[1] += dt * postProcData[2][i]; // p 
      L2_errors[2] += dt * postProcData[3][i]; // phi
      L2_errors[3] += dt * postProcData[4][i]; // mu

      H1_errors[0] += dt * postProcData[5][i];
      H1_errors[1] += dt * postProcData[6][i];
      H1_errors[2] += dt * postProcData[7][i];
      H1_errors[3] += dt * postProcData[8][i];
    }

    for(feSpace *s : spaces)
      delete s;
    delete solver;

    return FE_STATUS_OK;
  }
};

//
// Parameters for mass-averaged (ten Eikelder & Brunk) formulation of CHNS
//
#undef DECLARE_PARAMETERS
#define DECLARE_PARAMETERS(args, par)   \
  const double rhoA       = par[0];     \
  const double rhoB       = par[1];     \
  const double viscA      = par[2];     \
  const double viscB      = par[3];     \
  const double M          = par[4];     \
  const double gamma      = par[5];     \
  const double epsilon    = par[6];     \
  const double beta       = par[7];     \
  const double tau        = par[8];     \
  const double alpha      = par[9];     \
                                        \
  UNUSED(rhoA, rhoB,                    \
    viscA, viscB, M, gamma,             \
    epsilon, beta, tau, alpha);         \

//
// This mobility has to match the one in the model (-:
//
#define MOBILITY(M, phi) M * (1. - phi * phi)
#define DMOBILITY_DPHI(M, phi) -2. * M * phi

struct CHNS_Brunk_tenEikelder {

  static void uSource_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);
    DECLARE_PARAMETERS(args, par);

    const double u[2]       = UVEC(x,y,t);
    const double dudt[2]    = DUDTVEC(x,y,t);
    const double phiRef     = PHI(x,y,t);

    const double rho  = ( rhoA -  rhoB)/2. * phiRef + ( rhoA +  rhoB)/2.;
    const double visc = (viscA - viscB)/2. * phiRef + (viscA + viscB)/2.;

    const double drhodphi  = ( rhoA -  rhoB)/2.;
    const double dviscdphi = (viscA - viscB)/2.;

    const double gradu[2][2] = GRADU(x,y,t);
    const double gradp[2]    = GRADP(x,y,t);
    const double gradphi[2]  = GRADPHI(x,y,t);
    const double gradmu[2]   = GRADMU(x,y,t);

    const double uDotGradu[2] = {u[0]*gradu[0][0] + u[1]*gradu[1][0],
                                 u[0]*gradu[0][1] + u[1]*gradu[1][1]};
    const double f[2] = {0., -9.81};

    const double dphidt      = DPHIDT(x,y,t);
    const double div_u       = DIVU(x,y,t);
    const double uDotGradPhi = u[0]*gradphi[0] + u[1]*gradphi[1];
    const double div_rho_u   = rho * div_u + drhodphi * uDotGradPhi;
    const double drhodt      = drhodphi * dphidt;
    const double consMass    = drhodt + div_rho_u;

    const double dim     = 2.;
    const double d[2][2] = {{2. * gradu[0][0] - (2./dim) * div_u, (gradu[0][1] + gradu[1][0])},
                           {(gradu[0][1] + gradu[1][0]) , 2. * gradu[1][1] - (2./dim) * div_u}};

    const double gradPhidotd[2] = {gradphi[0]*d[0][0] + gradphi[1]*d[1][0],
                                   gradphi[0]*d[0][1] + gradphi[1]*d[1][1]};

    const double divd_incompressible[2] = DIV_RATE_STRAIN(x,y,t);
    const double grad_div_u[2]          = GRAD_DIV_U(x,y,t);
    const double div_d[2]               = {visc * (2. * divd_incompressible[0] - (2./dim) * grad_div_u[0]) + dviscdphi*gradPhidotd[0],
                                           visc * (2. * divd_incompressible[1] - (2./dim) * grad_div_u[1]) + dviscdphi*gradPhidotd[1]};

    res[0] = - (rho * (dudt[0] + uDotGradu[0] - f[0]) + u[0]/2. * consMass - div_d[0] + gradp[0] + phiRef*gradmu[0]);
    res[1] = - (rho * (dudt[1] + uDotGradu[1] - f[1]) + u[1]/2. * consMass - div_d[1] + gradp[1] + phiRef*gradmu[1]);
  }

  static double pSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    DECLARE_PARAMETERS(args, par);

    const double div_u = DIVU(x,y,t);
    const double lapP  = LAP_P(x,y,t);
    const double lapMu = LAP_MU(x,y,t);

    const double phiRef = PHI(x,y,t);
    const double Mphi   = MOBILITY(M, phiRef);
    const double dMdphi = DMOBILITY_DPHI(M,phiRef);

    const double gradPhi[2] = GRADPHI(x,y,t);
    const double gradMu[2]  = GRADMU(x,y,t);
    const double gradP[2]   = GRADP(x,y,t);
    const double gradPhi_dot_gradMuPlusAlphaP = gradPhi[0] * (gradMu[0] + alpha * gradP[0])
                                              + gradPhi[1] * (gradMu[1] + alpha * gradP[1]);

    const double div_M_grad_mu_plus_alpha_p = Mphi * (lapMu + alpha * lapP) + dMdphi * gradPhi_dot_gradMuPlusAlphaP;

    return - (div_u - alpha * div_M_grad_mu_plus_alpha_p);
  }

  static double phiSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    DECLARE_PARAMETERS(args, par);

    const double phiRef      = PHI(x,y,t);
    const double u[2]        = UVEC(x,y,t);
    const double gradphi[2]  = GRADPHI(x,y,t);
    const double dphidt      = DPHIDT(x,y,t);
    const double lapP        = LAP_P(x,y,t);
    const double lapMu       = LAP_MU(x,y,t);
    const double uDotGradPhi = u[0]*gradphi[0] + u[1]*gradphi[1];
    const double divu        = DIVU(x,y,t);
    const double div_phiu    = phiRef * divu + uDotGradPhi;

    const double Mphi   = MOBILITY(M, phiRef);
    const double dMdphi = DMOBILITY_DPHI(M,phiRef);

    const double gradPhi[2] = GRADPHI(x,y,t);
    const double gradMu[2]  = GRADMU(x,y,t);
    const double gradP[2]   = GRADP(x,y,t);
    const double gradPhi_dot_gradMuPlusAlphaP = gradPhi[0] * (gradMu[0] + alpha * gradP[0])
                                              + gradPhi[1] * (gradMu[1] + alpha * gradP[1]);

    return - (dphidt + div_phiu - (Mphi * (lapMu + alpha * lapP) + dMdphi * gradPhi_dot_gradMuPlusAlphaP));
  }

  static double muSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    DECLARE_PARAMETERS(args, par);

    const double mu      = MU(x,y,t);
    const double phiRef  = PHI(x,y,t);
    const double phiCube = PHI(x,y,t) * PHI(x,y,t) * PHI(x,y,t);
    const double lapPhi  = LAP_PHI(x,y,t);

    return - (mu + tau * lapPhi - beta * (phiCube - phiRef));
  }

  // CHNS_Brunk_tenEikelder
  static feStatus solve(const std::string &meshFile,
                        const double t0,
                        const double t1,
                        const int timeStepFactor,
                        const std::vector<int> &orders,
                        const int degreeQuadrature,
                        int &numInteriorElements,
                        int &numVertices,
                        std::vector<double> &L2_errors,
                        std::vector<double> &H1_errors,
                        int &timesteps)
  {
    const double rhoA    = 1000.;
    const double rhoB    = 1.;
    const double viscA   = 100.;
    const double viscB   = 1.;
    const double M       = 1e-6;
    const double gamma   = 0.01;
    const double epsilon = 1e-3;

    // const double rhoA    = 1.;
    // const double rhoB    = 1.;
    // const double viscA   = 1.;
    // const double viscB   = 1.;
    // const double M       = 1.;
    // const double gamma   = 1.;
    // const double epsilon = 1.;

    const double alpha = (rhoB - rhoA) / (rhoA + rhoB);
    const double beta   = 3. / (2. * sqrt(2.)) * gamma / epsilon;
    const double tau    = 3. / (2. * sqrt(2.)) * gamma * epsilon;

    std::vector<double> params = {rhoA, rhoB, viscA, viscB, M, gamma, epsilon, beta, tau, alpha};

    Parameters::CHNS CHNSparameters;
    CHNSparameters.formulation = Parameters::CHNS::Formulation::massAveraged;
    CHNSparameters.mobilityType = Parameters::CHNS::MobilityType::degenerate;
    // CHNSparameters.mobilityType = Parameters::CHNS::MobilityType::constant;
    CHNSparameters.fluids[0].density           = rhoA;
    CHNSparameters.fluids[0].dynamic_viscosity = viscA;
    CHNSparameters.fluids[1].density           = rhoB;
    CHNSparameters.fluids[1].dynamic_viscosity = viscB;
    CHNSparameters.gravity        = -9.81;
    CHNSparameters.surfaceTension = gamma;
    CHNSparameters.epsilon        = epsilon;
    CHNSparameters.mobility       = M;
    CHNSparameters.limitPhaseMarker = false;

    feVectorFunction uSol(  uSol_f, params);
    feFunction       pSol(  pSol_f, params);
    feFunction     phiSol(phiSol_f, params);
    feFunction      muSol( muSol_f, params);

    feVectorFunction   grad_uSol(  grad_uSol_f, params);
    feVectorFunction   grad_pSol(  grad_pSol_f, params);
    feVectorFunction grad_phiSol(grad_phiSol_f, params);
    feVectorFunction  grad_muSol( grad_muSol_f, params);

    feVectorFunction   uSource(  uSource_f, params);
    feFunction         pSource(  pSource_f, params);
    feFunction       phiSource(phiSource_f, params);
    feFunction        muSource( muSource_f, params); 

    const elementType      lag    = elementType::LAGRANGE;
    const elementType      vlag   = elementType::VECTOR_LAGRANGE;

    VectorFEDescriptor   u("U", "Domaine", vlag, orders[0], degreeQuadrature, &uSol);
    ScalarFEDescriptor   p("P", "Domaine", lag, orders[1], degreeQuadrature, &pSol);
    ScalarFEDescriptor phi("Phi", "Domaine", lag, orders[2], degreeQuadrature, &phiSol);
    ScalarFEDescriptor  mu("Mu", "Domaine", lag, orders[3], degreeQuadrature, &muSol);

    std::vector<FEDescriptor> fieldDescriptors = {u, p, phi, mu};

    // We have to be careful for the boundary conditions of this formulation:
    // integration by parts of the grad(mu + alpha*p) terms and discarding the boundary
    // term assumes that either 
    //    - grad(mu + alpha*p) cdot n = 0 on the boundaries, or
    //    - the test functions of the equation in which this term appears vanish on the
    //      boundaries, that is, the equation variable (here, p or mu) is set on the boundary.
    // Physical applications (in Brunk and ten Eikelder's preprint)
    // assume grad(mu + alpha*p) cdot n = 0.
    // Here, mu is set on the boundaries, so we also set p on the boundaries.
    // We can't set only a pressure point, because then the pressure equation is wrong.
    VectorFEDescriptor uB("U", "Bord", vlag, orders[0], degreeQuadrature, &uSol);
    // ScalarFEDescriptor pP("P", "PointPression", lag, orders[1], degreeQuadrature, &pSol);
    ScalarFEDescriptor pP("P", "Bord", lag, orders[1], degreeQuadrature, &pSol);
    ScalarFEDescriptor phiB("Phi", "Bord", lag, orders[2], degreeQuadrature, &phiSol);
    ScalarFEDescriptor muB("Mu", "Bord", lag, orders[3], degreeQuadrature, &muSol);

    BoundaryConditions::VectorDirichlet uBC(uB);
    BoundaryConditions::ScalarDirichlet pBC(pP), phiBC(phiB), muBC(muB);

    std::vector<const BoundaryConditions::BoundaryCondition*> BC = {&uBC, &pBC, &phiBC, &muBC};

    Parameters::NonLinearSolver newton_param;

    Parameters::TimeIntegration time_param;
    time_param.t_initial  = t0;
    time_param.t_final    = t1;

    #if defined(TEST_TIME_CONVERGENCE)
      time_param.nTimeSteps = 10 * pow(2, timeStepFactor);
    #else 
      time_param.nTimeSteps = 10 * pow(2, timeStepFactor);
    #endif

    timesteps = time_param.nTimeSteps;

    time_param.method =
      Parameters::TimeIntegration::TimeIntegrationMethod::bdf2;
    time_param.bdf2starter =
      Parameters::TimeIntegration::BDF2Starter::exactSolution;

    CHNS_Solver *solver;
    feCheck(createCHNS_Solver(solver,
                              CHNSparameters,
                              newton_param,
                              time_param,
                              fieldDescriptors,
                              BC,
                              nullptr,
                              &pSource,
                              &uSource,
                              &phiSource,
                              &muSource));

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();
    numVertices = mesh.getNumVertices();

    std::vector<feSpace *> spaces, essentialSpaces;
    solver->createSpaces(&mesh, spaces, essentialSpaces);

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution   sol(numbering.getNbDOFs(), spaces, essentialSpaces);

    feExporter *exporter = nullptr;
    feExportData exportData = {exporter, 1, "sol"};

    feNorm *uL2, *uH1, *pL2, *pH1, *phiL2, *phiH1, *muL2, *muH1;
    feCheckReturn(createNorm(uL2, VECTOR_L2_ERROR, {spaces[0]}, &sol, nullptr, &uSol));
    feCheckReturn(createNorm(uH1, VECTOR_SEMI_H1_ERROR, {spaces[0]}, &sol, nullptr, &grad_uSol));
    feCheckReturn(createNorm(pL2, L2_ERROR, {spaces[1]}, &sol, &pSol));
    feCheckReturn(createNorm(pH1, SEMI_H1_ERROR, {spaces[1]}, &sol, nullptr, &grad_pSol));
    feCheckReturn(createNorm(phiL2, L2_ERROR, {spaces[2]}, &sol, &phiSol));
    feCheckReturn(createNorm(phiH1, SEMI_H1_ERROR, {spaces[2]}, &sol, nullptr, &grad_phiSol));
    feCheckReturn(createNorm(muL2, L2_ERROR, {spaces[3]}, &sol, &muSol));
    feCheckReturn(createNorm(muH1, SEMI_H1_ERROR, {spaces[3]}, &sol, nullptr, &grad_muSol));

    std::vector<feNorm*> norms = {uL2, pL2, phiL2, muL2, uH1, pH1, phiH1, muH1};
    solver->setNorms(norms);

    feCheck(solver->solve(
      &mesh, &sol, &numbering, spaces, exportData));

    std::vector<std::vector<double>> postProcData = solver->getPostProcessingData();

    for(unsigned int i = 0; i < L2_errors.size(); ++i)
    {
      L2_errors[i] = 0.;
      H1_errors[i] = 0.;
    }

    // First post proc vector is the time
    const double dt = postProcData[0][1] - postProcData[0][0];

    // Compute L1 norm in time of the computed spatial errors
    for(size_t i = 0; i < postProcData[0].size(); ++i)
    {
      L2_errors[0] += dt * postProcData[1][i]; // u 
      L2_errors[1] += dt * postProcData[2][i]; // p 
      L2_errors[2] += dt * postProcData[3][i]; // phi
      L2_errors[3] += dt * postProcData[4][i]; // mu

      H1_errors[0] += dt * postProcData[5][i];
      H1_errors[1] += dt * postProcData[6][i];
      H1_errors[2] += dt * postProcData[7][i];
      H1_errors[3] += dt * postProcData[8][i];
    }

    // // Error at last time step only
    // L2_errors[0] = postProcData[1].back(); // u 
    // L2_errors[1] = postProcData[2].back(); // p 
    // L2_errors[2] = postProcData[3].back(); // phi
    // L2_errors[3] = postProcData[4].back(); // mu

    // H1_errors[0] = postProcData[5].back();
    // H1_errors[1] = postProcData[6].back();
    // H1_errors[2] = postProcData[7].back();
    // H1_errors[3] = postProcData[8].back();

    for(feSpace *s : spaces)
      delete s;
    delete solver;

    return FE_STATUS_OK;
  }
}; // CHNS_Brunk_tenEikelder

//
// Parameters for volume-averaged (Khanwale et al.) formulation of CHNS
//
#undef DECLARE_PARAMETERS
#define DECLARE_PARAMETERS(args, par)   \
  const double rhoA   = par[0];     \
  const double rhoB   = par[1];     \
  const double viscA  = par[2];     \
  const double viscB  = par[3];     \
  const double Re     = par[4];     \
  const double Pe     = par[5];     \
  const double Cn     = par[6];     \
  const double We     = par[7];     \
  const double Fr     = par[8];     \
  const double M      = par[9];     \
                                    \
  UNUSED(rhoA, rhoB,  viscA,        \
    viscB, Re, Pe, Cn, We, Fr, M);  \

struct CHNS_KhanwaleEtAl {

  static void uSource_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);
    DECLARE_PARAMETERS(args, par);

    res[0] = 0.;
    res[1] = 0.;
  }

  static double pSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    DECLARE_PARAMETERS(args, par);

    const double div_u       = DIVU(x,y,t);

    const double phiRef      = PHI(x,y,t);
    const double drhodphi    = ( rhoA -  rhoB)/2.;
    const double rho         = ( rhoA -  rhoB)/2. * phiRef + ( rhoA +  rhoB)/2.;

    const double dphidt      = DPHIDT(x,y,t);
    const double drhodt      = drhodphi * dphidt;

    const double u[2]        = UVEC(x,y,t);
    const double gradphi[2]  = GRADPHI(x,y,t);
    const double uDotGradPhi = u[0]*gradphi[0] + u[1]*gradphi[1];
    const double div_rho_u   = rho * div_u + drhodphi * uDotGradPhi;

    const double lapMu       = LAP_MU(x,y,t);
    const double div_j       = (rhoB - rhoA) / (2. * rhoA * Cn) * lapMu;

    return - (drhodt + div_rho_u + 1./Pe * div_j);
    return - (div_u  + drhodt + div_rho_u + 1./Pe * div_j);
  }

  static double phiSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    DECLARE_PARAMETERS(args, par);

    const double phiRef      = PHI(x,y,t);
    const double u[2]        = UVEC(x,y,t);
    const double gradphi[2]  = GRADPHI(x,y,t);
    const double dphidt      = DPHIDT(x,y,t);
    const double lapMu       = LAP_MU(x,y,t);
    const double uDotGradPhi = u[0]*gradphi[0] + u[1]*gradphi[1];
    const double divu        = DIVU(x,y,t);
    const double div_phiu    = phiRef * divu + uDotGradPhi;


    return - (dphidt + div_phiu - M / (Pe * Cn) * lapMu);
  }

  static double muSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    DECLARE_PARAMETERS(args, par);

    const double mu      = MU(x,y,t);
    const double phiRef  = PHI(x,y,t);
    const double phiCube = PHI(x,y,t) * PHI(x,y,t) * PHI(x,y,t);
    const double lapPhi  = LAP_PHI(x,y,t);

    return - (mu - (phiCube - phiRef) + Cn * Cn * lapPhi);
  }

  // CHNS_KhanwaleEtAl
  static feStatus solve(const std::string &meshFile,
                        const double t0,
                        const double t1,
                        const int timeStepFactor,
                        const std::vector<int> &orders,
                        const int degreeQuadrature,
                        int &numInteriorElements,
                        int &numVertices,
                        std::vector<double> &L2_errors,
                        std::vector<double> &H1_errors,
                        int &timesteps)
  {
    const double rhoA    = 1.;
    const double rhoB    = 1.;
    const double viscA   = 1.;
    const double viscB   = 1.;
    const double Re      = 1.;
    const double Pe      = 1.;
    const double Cn      = 10.;
    const double We      = 1.;
    const double Fr      = 1.;
    const double M       = 1.;

    std::vector<double> params = {rhoA, rhoB, viscA, viscB, Re, Pe, Cn, We, Fr, M};

    Parameters::CHNS CHNSparameters;
    CHNSparameters.formulation = Parameters::CHNS::Formulation::khanwale;
    CHNSparameters.mobilityType = Parameters::CHNS::MobilityType::constant;

    CHNSparameters.fluids[0].density           = rhoA;
    CHNSparameters.fluids[0].dynamic_viscosity = viscA;
    CHNSparameters.fluids[1].density           = rhoB;
    CHNSparameters.fluids[1].dynamic_viscosity = viscB;
    
    CHNSparameters.Re = Re;
    CHNSparameters.Pe = Pe;
    CHNSparameters.Cn = Cn;
    CHNSparameters.We = We;
    CHNSparameters.Fr = Fr;
      
    CHNSparameters.mobility = M;

    CHNSparameters.limitPhaseMarker = false;

    feVectorFunction uSol(  uSol_f, params);
    feFunction       pSol(  pSol_f, params);
    feFunction     phiSol(phiSol_f, params);
    feFunction      muSol( muSol_f, params);

    feVectorFunction   grad_uSol(  grad_uSol_f, params);
    feVectorFunction   grad_pSol(  grad_pSol_f, params);
    feVectorFunction grad_phiSol(grad_phiSol_f, params);
    feVectorFunction  grad_muSol( grad_muSol_f, params);

    feVectorFunction   uSource(  uSource_f, params);
    feFunction         pSource(  pSource_f, params);
    feFunction       phiSource(phiSource_f, params);
    feFunction        muSource( muSource_f, params); 

    const elementType      lag    = elementType::LAGRANGE;
    const elementType      vlag   = elementType::VECTOR_LAGRANGE;

    VectorFEDescriptor   u("U", "Domaine", vlag, orders[0], degreeQuadrature, &uSol);
    ScalarFEDescriptor   p("P", "Domaine", lag, orders[1], degreeQuadrature, &pSol);
    ScalarFEDescriptor phi("Phi", "Domaine", lag, orders[2], degreeQuadrature, &phiSol);
    ScalarFEDescriptor  mu("Mu", "Domaine", lag, orders[3], degreeQuadrature, &muSol);

    std::vector<FEDescriptor> fieldDescriptors = {u, p, phi, mu};

    VectorFEDescriptor uB("U", "Domaine", vlag, orders[0], degreeQuadrature, &uSol);
    ScalarFEDescriptor pP("P", "Domaine", lag, orders[1], degreeQuadrature, &pSol);

    // VectorFEDescriptor uB("U", "Bord", vlag, orders[0], degreeQuadrature, &uSol);
    // ScalarFEDescriptor pP("P", "PointPression", lag, orders[1], degreeQuadrature, &pSol);
    // ScalarFEDescriptor pP("P", "Bord", lag, orders[1], degreeQuadrature, &pSol);
    ScalarFEDescriptor phiB("Phi", "Bord", lag, orders[2], degreeQuadrature, &phiSol);
    ScalarFEDescriptor muB("Mu", "Bord", lag, orders[3], degreeQuadrature, &muSol);

    BoundaryConditions::VectorDirichlet uBC(uB);
    BoundaryConditions::ScalarDirichlet pBC(pP), phiBC(phiB), muBC(muB);

    std::vector<const BoundaryConditions::BoundaryCondition*> BC = {&uBC, &pBC, &phiBC, &muBC};

    Parameters::NonLinearSolver newton_param;

    Parameters::TimeIntegration time_param;
    time_param.t_initial  = t0;
    time_param.t_final    = t1;

    #if defined(TEST_TIME_CONVERGENCE)
      time_param.nTimeSteps = 10 * pow(2, timeStepFactor);
    #else 
      time_param.nTimeSteps = 10 * pow(2, timeStepFactor);
    #endif

    timesteps = time_param.nTimeSteps;

    time_param.method =
      Parameters::TimeIntegration::TimeIntegrationMethod::bdf1;
    time_param.bdf2starter =
      Parameters::TimeIntegration::BDF2Starter::exactSolution;

    CHNS_Solver *solver;
    feCheck(createCHNS_Solver(solver,
                              CHNSparameters,
                              newton_param,
                              time_param,
                              fieldDescriptors,
                              BC,
                              nullptr,
                              &pSource,
                              &uSource,
                              &phiSource,
                              &muSource));

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();
    numVertices = mesh.getNumVertices();

    std::vector<feSpace *> spaces, essentialSpaces;
    solver->createSpaces(&mesh, spaces, essentialSpaces);

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution   sol(numbering.getNbDOFs(), spaces, essentialSpaces);

    feExporter *exporter = nullptr;
    feExportData exportData = {exporter, 1, "sol"};

    feNorm *uL2, *uH1, *pL2, *pH1, *phiL2, *phiH1, *muL2, *muH1;
    feCheckReturn(createNorm(uL2, VECTOR_L2_ERROR, {spaces[0]}, &sol, nullptr, &uSol));
    feCheckReturn(createNorm(uH1, VECTOR_SEMI_H1_ERROR, {spaces[0]}, &sol, nullptr, &grad_uSol));
    feCheckReturn(createNorm(pL2, L2_ERROR, {spaces[1]}, &sol, &pSol));
    feCheckReturn(createNorm(pH1, SEMI_H1_ERROR, {spaces[1]}, &sol, nullptr, &grad_pSol));
    feCheckReturn(createNorm(phiL2, L2_ERROR, {spaces[2]}, &sol, &phiSol));
    feCheckReturn(createNorm(phiH1, SEMI_H1_ERROR, {spaces[2]}, &sol, nullptr, &grad_phiSol));
    feCheckReturn(createNorm(muL2, L2_ERROR, {spaces[3]}, &sol, &muSol));
    feCheckReturn(createNorm(muH1, SEMI_H1_ERROR, {spaces[3]}, &sol, nullptr, &grad_muSol));

    std::vector<feNorm*> norms = {uL2, pL2, phiL2, muL2, uH1, pH1, phiH1, muH1};
    solver->setNorms(norms);

    feCheck(solver->solve(
      &mesh, &sol, &numbering, spaces, exportData));

    std::vector<std::vector<double>> postProcData = solver->getPostProcessingData();

    for(unsigned int i = 0; i < L2_errors.size(); ++i)
    {
      L2_errors[i] = 0.;
      H1_errors[i] = 0.;
    }

    // First post proc vector is the time
    const double dt = postProcData[0][1] - postProcData[0][0];

    // Compute L1 norm in time of the computed spatial errors
    for(size_t i = 0; i < postProcData[0].size(); ++i)
    {
      L2_errors[0] += dt * postProcData[1][i]; // u 
      L2_errors[1] += dt * postProcData[2][i]; // p 
      L2_errors[2] += dt * postProcData[3][i]; // phi
      L2_errors[3] += dt * postProcData[4][i]; // mu

      H1_errors[0] += dt * postProcData[5][i];
      H1_errors[1] += dt * postProcData[6][i];
      H1_errors[2] += dt * postProcData[7][i];
      H1_errors[3] += dt * postProcData[8][i];
    }

    // Error at last time step only
    L2_errors[0] = postProcData[1].back(); // u 
    L2_errors[1] = postProcData[2].back(); // p 
    L2_errors[2] = postProcData[3].back(); // phi
    L2_errors[3] = postProcData[4].back(); // mu

    H1_errors[0] = postProcData[5].back();
    H1_errors[1] = postProcData[6].back();
    H1_errors[2] = postProcData[7].back();
    H1_errors[3] = postProcData[8].back();

    for(feSpace *s : spaces)
      delete s;
    delete solver;

    return FE_STATUS_OK;
  }
}; // CHNS_KhanwaleEtAl

template<typename problemType>
feStatus meshConvergence(std::stringstream &resultBuffer,
                         const std::string &meshFileRoot,
                         const std::vector<int> &orders,
                         const int numMeshes,
                         const int degreeQuadrature)
{
  std::vector<int> nElm(numMeshes);
  std::vector<int> timesteps(numMeshes);
  std::vector<int> nVertices(numMeshes);

  const unsigned int nFields = 4;
  std::vector<std::vector<double>> allL2_errors(nFields, std::vector<double>(numMeshes, 0.));
  std::vector<std::vector<double>> allH1_errors(nFields, std::vector<double>(numMeshes, 0.));

  bool writeErrorToFile = false;
  std::string outputDirectory = "./tests/withLinearSolver/";
  std::ofstream errorFileL2, errorFileH1, timeErrorFile;
  if(writeErrorToFile) {
  #if defined(TEST_TIME_CONVERGENCE)
    timeErrorFile.open(outputDirectory + "timeErrorFile.txt");
  #else
    errorFileL2.open(outputDirectory + "errorL2.txt");
    errorFileH1.open(outputDirectory + "errorH1.txt");
  #endif
  }

  const double t0 = 0.;
  const double t1 = 1.;

  for(int i = 0; i < numMeshes; ++i)
  {
    const int timeStepFactor = i;

    #if defined(TEST_TIME_CONVERGENCE)
      std::string meshFile = meshFileRoot + "1.msh";
    #else
      std::string meshFile = meshFileRoot + std::to_string(i+1) + ".msh";
    #endif

    std::vector<double> L2_errors(nFields), H1_errors(nFields);

    feCheckReturn(problemType::solve(meshFile, t0, t1, timeStepFactor, orders,
                       degreeQuadrature, nElm[i], nVertices[i],
                       L2_errors, H1_errors, timesteps[i]));

    for(unsigned int iField = 0; iField < nFields; ++iField)
    {
      allL2_errors[iField][i] = L2_errors[iField]; 
      allH1_errors[iField][i] = H1_errors[iField]; 
    }

    if(writeErrorToFile) {
    #if defined(TEST_TIME_CONVERGENCE)
      //
      // Print semi H1 errors for u, phi, mu and L2 error for p.
      //
      timeErrorFile << timesteps[i] << "\t" << nElm[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << allH1_errors[0][i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << allL2_errors[1][i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << allH1_errors[2][i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << allH1_errors[3][i] << "\n";
      timeErrorFile.flush();
    #else
      errorFileL2 << nVertices[i] << "\t" << nElm[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << allL2_errors[0][i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << allL2_errors[1][i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << allL2_errors[2][i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << allL2_errors[3][i] << "\n";
      errorFileL2.flush();
      errorFileH1 << nVertices[i] << "\t" << nElm[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << allH1_errors[0][i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << allH1_errors[1][i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << allH1_errors[2][i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << allH1_errors[3][i] << "\n";
      errorFileH1.flush();
    #endif
    }
  }

  if(writeErrorToFile) {
  #if defined(TEST_TIME_CONVERGENCE)
    timeErrorFile.close();
  #else
    errorFileL2.close();
    errorFileH1.close();
  #endif
  }

  resultBuffer << "Cahn-Hilliard Navier-Stokes MMS" << std::endl;
  resultBuffer << "   (u,p) : Taylor-Hood elements P" << orders[0] << "-P" << orders[1] << std::endl; 
  resultBuffer << "(phi,mu) : Equal order elements P" << orders[2] << "-P" << orders[3] << std::endl;
  resultBuffer << "Initial time : " << t0 << std::endl;
  resultBuffer << "Final   time : " << t1 << std::endl;
  resultBuffer << "N timesteps for first mesh : " << timesteps[0] << std::endl;
  resultBuffer << "N timesteps for last  mesh : " << timesteps.back() << std::endl;

  resultBuffer << "==================== L2 error norms ==================" << std::endl;
  resultBuffer << "U Velocity" << std::endl;
  computeAndPrintConvergence(2, numMeshes, allL2_errors[0], nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  resultBuffer << "P Pressure" << std::endl;
  computeAndPrintConvergence(2, numMeshes, allL2_errors[1], nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  resultBuffer << "Phi Phase marker" << std::endl;
  computeAndPrintConvergence(2, numMeshes, allL2_errors[2], nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  resultBuffer << "Mu Chemical potential" << std::endl;
  computeAndPrintConvergence(2, numMeshes, allL2_errors[3], nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

  resultBuffer << "==================== H1 error norms ==================" << std::endl;
  resultBuffer << "U Velocity" << std::endl;
  computeAndPrintConvergence(2, numMeshes, allH1_errors[0], nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  resultBuffer << "P Pressure" << std::endl;
  computeAndPrintConvergence(2, numMeshes, allH1_errors[1], nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  resultBuffer << "Phi Phase marker" << std::endl;
  computeAndPrintConvergence(2, numMeshes, allH1_errors[2], nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  resultBuffer << "Mu Chemical potential" << std::endl;
  computeAndPrintConvergence(2, numMeshes, allH1_errors[3], nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

  return FE_STATUS_OK;
}

TEST(CahnHilliard, CH)
{
  initialize(my_argc, my_argv);
  setVerbose(0);
  
  std::string testRoot = "../../../tests/withLinearSolver/CH_MMS";

  std::stringstream resultBuffer;

  int degreeQuadrature = 8;
  int order = 1;
  cahnHilliard::meshConvergence(resultBuffer, order, 4, degreeQuadrature);
  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

TEST(CahnHilliard, CHNS_VolumeAveraged_Abels)
{
  initialize(my_argc, my_argv);
  setVerbose(1);
  
  std::string testRoot = "../../../tests/withLinearSolver/CHNS_Abels";
  std::string meshFileRoot = "../../../data/mmsMeshes/square";

  std::stringstream resultBuffer;

  int degreeQuadrature = 8;
  std::vector<int> orders = {2, 1, 1, 1};
  ASSERT_TRUE(meshConvergence<CHNS_VolumeAveraged_Abels>(resultBuffer, 
    meshFileRoot, orders, 3, degreeQuadrature) == FE_STATUS_OK);
  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

TEST(CahnHilliard, CHNS_Brunk_tenEikelder)
{
  initialize(my_argc, my_argv);
  setVerbose(1);
  
  std::string testRoot = "../../../tests/withLinearSolver/CHNS_BtE";
  std::string meshFileRoot = "../../../data/mmsMeshes/square";

  std::stringstream resultBuffer;

  int degreeQuadrature = 8;
  std::vector<int> orders = {2, 1, 1, 1};
  ASSERT_TRUE(meshConvergence<CHNS_Brunk_tenEikelder>(resultBuffer,
    meshFileRoot, orders, 3, degreeQuadrature) == FE_STATUS_OK);
  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

TEST(CahnHilliard, CHNS_KhanwaleEtAl)
{
  initialize(my_argc, my_argv);
  setVerbose(1);
  
  std::string testRoot = "../../../tests/withLinearSolver/CHNS_Khanwale";
  std::string meshFileRoot = "../../../data/mmsMeshes/square";

  std::stringstream resultBuffer;

  int degreeQuadrature = 8;
  std::vector<int> orders = {2, 1, 1, 1};
  ASSERT_TRUE(meshConvergence<CHNS_KhanwaleEtAl>(resultBuffer,
    meshFileRoot, orders, 3, degreeQuadrature) == FE_STATUS_OK);
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