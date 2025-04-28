
#include "feAPI.h"
#include "../tests.h"

#include <gtest/gtest.h>

/*
  Tests the implementation of the Cahn-Hilliard and CH-Navier-Stokes related weak forms:
    - CHNS_Momentum
    - TracerConvection
    - MixedGradGrad
    - MixedMassPower
    - MixedMass
*/

static int my_argc;
static char** my_argv;

//
// Unsteady Cahn-Hilliard system without Navier-Stokes
// Lagrange elements of equal order for (phi,mu) (phase marker and chemical potential)
//
namespace cahnHilliard {

  namespace MMS {

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

    void uSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/, std::vector<double> &res)
    {
      const double x = args.pos[0];
      const double y = args.pos[1];
      const double minus_dpdx[2] = {-2.*x*y*y, -2.*x*x*y};
      const double lap_u[2] = {12. * (x*x*y*y*y*y + x*x*x*x*y*y),
                              -4./5. * (6. *x*y*y*y*y*y + 20. * x*x*x*y*y*y)};
      res[0] = - (minus_dpdx[0] + lap_u[0]);
      res[1] = - (minus_dpdx[1] + lap_u[1]);
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
        feCheck(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
        feCheck(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
      #else
        feCheck(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
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
        cahnHilliard::MMS::solve(meshFile, order, degreeQuadrature, NLoptions, nElm[i], errPhi[i], errMu[i]);
      }
      resultBuffer << "(phi,mu) : Equal order elements P" << order << "-P" << order << std::endl; 
      resultBuffer << "Cahn-Hilliard MMS - Error on phase marker phi" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errPhi, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      resultBuffer << "Cahn-Hilliard MMS - Error on chemical potential mu" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errMu, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

      return 0;
    }

  } // namespace MMS
}

namespace cahnHilliardNavierStokes {

  namespace MMS {

    void uSol_f(const feFunctionArguments &args, const std::vector<double> & /* par */, std::vector<double> &res)
    {
      const double x = args.pos[0];
      const double y = args.pos[1];
      const double t = args.t;
      const double x3 = x*x*x;
      const double x4 = x*x*x*x;
      const double y4 = y*y*y*y;
      const double y5 = y*y*y*y*y;
      res[0] = exp(-t) * x4*y4;
      res[1] = exp(-t) * (-4./5. * x3*y5);
    }

    double pSol_f(const feFunctionArguments &args, const std::vector<double> & /* par */)
    {
      const double x = args.pos[0];
      const double y = args.pos[1];
      return x*x*y*y;
    }

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

    double rho_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      // Assumes phi is in [-1,1]
      // Property of A when phi = 1, of B when phi = -1
      const double phi = args.u;
      const double rhoA = par[0];
      const double rhoB = par[1];
      return (rhoA - rhoB)/2. * phi + (rhoA + rhoB)/2.;
    }

    double drhodphi_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      UNUSED(args); // rho(phi) is linear here
      const double rhoA = par[0];
      const double rhoB = par[1];
      return (rhoA - rhoB)/2.;
    }

    double viscosity_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      // Assumes phi is in [-1,1]
      // Property of A when phi = 1, of B when phi = -1
      const double phi = args.u;
      const double viscA = par[0];
      const double viscB = par[1];
      return (viscA - viscB)/2. * phi + (viscA + viscB)/2.;
    }

    double dviscdphi_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      UNUSED(args); // visc(phi) is linear here
      const double viscA = par[0];
      const double viscB = par[1];
      return (viscA - viscB)/2.;
    }

    void uSource_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
    {
      const double x = args.pos[0];
      const double y = args.pos[1];
      const double t = args.t;

      const double x2 = x*x;
      const double x3 = x*x*x;
      const double x4 = x*x*x*x;

      const double y2 = y*y;
      const double y3 = y*y*y;
      const double y4 = y*y*y*y;
      const double y5 = y*y*y*y*y;

      const double rhoA  = par[0];
      const double rhoB  = par[1];
      const double viscA = par[2];
      const double viscB = par[3];
      const double M     = par[4];
      const double gamma = par[5];
      const double epsilon = par[6];
      bool isStationary = par[7] > 0;

      const double emt = exp(-t);
      const double sx = sin(M_PI*x);
      const double sy = sin(M_PI*y);

      const double u[2] = {emt * x4*y4, emt * (-4./5. * x3*y5)};
      const double phi = emt * sx * sy;
      const double mu  =       sx * sy;

      const double gradphi[2] = {emt * M_PI * cos(M_PI*x) * sin(M_PI*y),
                           emt * M_PI * sin(M_PI*x) * cos(M_PI*y)};
      const double gradmu[2]  =       {M_PI * cos(M_PI*x) * sin(M_PI*y),
                                 M_PI * sin(M_PI*x) * cos(M_PI*y)};

      const double rho  = ( rhoA -  rhoB)/2. * phi + ( rhoA +  rhoB)/2.;
      const double visc = (viscA - viscB)/2. * phi + (viscA + viscB)/2.;

      const double dudt[2] = {-emt * x4*y4,
                        -emt * (-4./5. * x3*y5)};

      // Convention: defined as gradu_ij = du_j/dx_i to agree with directional derivatives
      const double gradu[2][2] = {{emt *  4.*x3*y4, emt * -12.*x2*y5/5.},
                            {emt *  4.*x4*y3, emt * - 4.*x3*y4}};

      const double gradp[2] = {2.*x*y2, 2*x2*y};
      const double divd[2] = {  6.*emt*x4*y2 +   6.*emt*x2*y4,
                        - 8.*emt*x3*y3 - (12.*emt*x *y5)/5.};
      const double uDotGradu[2] = {u[0]*gradu[0][0] + u[1]*gradu[1][0],
                             u[0]*gradu[0][1] + u[1]*gradu[1][1]};
      const double gMudgu[2] = {gradmu[0]*gradu[0][0] + gradmu[1]*gradu[1][0],
                          gradmu[0]*gradu[0][1] + gradmu[1]*gradu[1][1]};
      const double f[2] = {0., -9.81};

      const double drhodphi = (rhoA-rhoB)/2.;
      const double korteweg = - gamma/epsilon;
      const double dviscdphi = (viscA-viscB)/2.;

      const double d[2][2] = {{gradu[0][0], (gradu[0][1] + gradu[1][0])/2.},
                        {(gradu[0][1] + gradu[1][0])/2., gradu[1][1]}};
      const double gradPhidotd[2] = {gradphi[0]*d[0][0] + gradphi[1]*d[1][0],
                               gradphi[0]*d[0][1] + gradphi[1]*d[1][1]};

      if(isStationary) {
        // For stationary test, do not include dudt because it will not be computed by the solver (BDF coefficient is 0).
        // If dudt is included, it adds du/dt at t = 0 to the source term, which is wrong with a stationary solver.
        UNUSED(dudt);
        res[0] = - (rho * (          uDotGradu[0] - f[0]) + M * drhodphi * gMudgu[0] + gradp[0] - 2.*(dviscdphi*gradPhidotd[0] + visc*divd[0]) + korteweg*mu*gradphi[0]);
        res[1] = - (rho * (          uDotGradu[1] - f[1]) + M * drhodphi * gMudgu[1] + gradp[1] - 2.*(dviscdphi*gradPhidotd[1] + visc*divd[1]) + korteweg*mu*gradphi[1]);
      } else {
        res[0] = - (rho * (dudt[0] + uDotGradu[0] - f[0]) + M * drhodphi * gMudgu[0] + gradp[0] - 2.*(dviscdphi*gradPhidotd[0] + visc*divd[0]) + korteweg*mu*gradphi[0]);
        res[1] = - (rho * (dudt[1] + uDotGradu[1] - f[1]) + M * drhodphi * gMudgu[1] + gradp[1] - 2.*(dviscdphi*gradPhidotd[1] + visc*divd[1]) + korteweg*mu*gradphi[1]);
      }
    }

    double pSource_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      const double x = args.pos[0];
      const double y = args.pos[1];
      const double t = args.t;

      const double x2 = x*x;
      const double x3 = x*x*x;
      const double x4 = x*x*x*x;
      const double y3 = y*y*y;
      const double y4 = y*y*y*y;
      const double y5 = y*y*y*y*y;
      const double sx = sin(M_PI*x);
      const double sy = sin(M_PI*y);
      const double emt = exp(-t);
      const double phi = emt * sx*sy;
      const double gradphi[2] = {emt * M_PI * cos(M_PI*x) * sin(M_PI*y),
                           emt * M_PI * sin(M_PI*x) * cos(M_PI*y)};
      const double u[2] = {emt * x4*y4, emt * (-4./5. * x3*y5)};
      const double gradu[2][2] = {{emt *  4.*x3*y4, emt * -12.*x2*y5/5.},
                            {emt *  4.*x4*y3, emt * - 4.*x3*y4}};

      const double rhoA  = par[0];
      const double rhoB  = par[1];
      const double rho      = (rhoA - rhoB)/2. * phi + (rhoA + rhoB)/2.;
      const double drhodphi = (rhoA - rhoB)/2.;

      const double div_u = gradu[0][0] + gradu[1][1]; // = 0
      const double divrhou = drhodphi * (gradphi[0]*u[0] + gradphi[1]*u[1]) + rho * div_u;

      return -divrhou;
    }

    double phiSource_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      bool isStationary = par[0] > 0;
      const double M = par[1];
      const double x = args.pos[0];
      const double y = args.pos[1];
      const double t = args.t;
      const double x3 = x*x*x;
      const double x4 = x*x*x*x;
      const double y4 = y*y*y*y;
      const double y5 = y*y*y*y*y;
      const double emt = exp(-t);
      const double sx = sin(M_PI*x);
      const double sy = sin(M_PI*y);

      const double u[2] = {emt * x4*y4, emt * (-4./5. * x3*y5)};
      const double gradphi[2] = {emt * M_PI * cos(M_PI*x) * sin(M_PI*y),
                           emt * M_PI * sin(M_PI*x) * cos(M_PI*y)};

      const double dphidt = - emt * sx*sy;
      const double lapMu  = -2.*M_PI*M_PI*sx*sy;
      const double uDotGradPhi = u[0]*gradphi[0] + u[1]*gradphi[1];

      if(isStationary) {
        // For stationary test, do not include dphidt
        UNUSED(dphidt);
        return - (uDotGradPhi - M*lapMu);
      } else {
        return - (dphidt + uDotGradPhi - M*lapMu);
      }
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

      return - (mu - (phiCube - phi) + epsilon*epsilon*lapPhi);
    }

    int solve(const std::string &meshFile,
              const int orderU,
              const int orderPhiMu,
              const int degreeQuadrature,
              const feNLSolverOptions &NLoptions,
              int &numInteriorElements,
              int &numVertices,
              double &L2ErrorU,
              double &L2ErrorP,
              double &L2ErrorPhi,
              double &L2ErrorMu)
              // ,
              // double &H1ErrorU,
              // double &H1ErrorP,
              // double &H1ErrorPhi,
              // double &H1ErrorMu
              // )
    {
      feVectorFunction uSol(uSol_f);
      feFunction   pSol(pSol_f);
      feFunction phiSol(phiSol_f);
      feFunction  muSol( muSol_f);

      const bool isStationary = false;
      const double rhoA = 10.;
      const double rhoB = 1.;
      const double viscA = 10.;
      const double viscB = 1.;
      const double mobility_val = 1e-1;
      const double gamma = 1e-1;
      const double epsilon = 1.;
      const double stationaryVal = isStationary ? 1. : -1.;
      feVectorFunction uSource(uSource_f, {rhoA, rhoB, viscA, viscB, mobility_val, gamma, epsilon, stationaryVal});
      feFunction pSource(pSource_f, {rhoA, rhoB});
      feFunction phiSource(phiSource_f, {stationaryVal, mobility_val});
      feFunction  muSource( muSource_f, {epsilon});

      feConstantFunction minusEpsInterfaceSquared(-epsilon*epsilon);
      feConstantFunction mobility(mobility_val);
      feFunction density(rho_f, {rhoA, rhoB});
      feFunction drhodphi(drhodphi_f, {rhoA, rhoB});
      feConstantFunction coeffKorteweg(-gamma/epsilon);
      feFunction viscosity(viscosity_f, {viscA, viscB});
      feFunction dviscdphi(dviscdphi_f, {viscA, viscB});
      feConstantVectorFunction volumeForce({0, -9.81});

      feMesh2DP1 mesh(meshFile);
      numInteriorElements = mesh.getNumInteriorElements();
      numVertices = mesh.getNumVertices();

      // Velocity and pressure
      feSpace *u = nullptr, *uB = nullptr, *p = nullptr, *pP = nullptr;
      // Phase marker and chemical potential
      feSpace *phi = nullptr, *phiB = nullptr, *mu = nullptr,  *muB = nullptr;

      // Initialize with exact initial condition
      // Use Taylor-Hood elements for (u,p) and equal order for (phi,u)
      int dQ = degreeQuadrature;
      feCheck(createFiniteElementSpace( u, &mesh, elementType::VECTOR_LAGRANGE, orderU, "U",       "Domaine", dQ, &uSol));
      feCheck(createFiniteElementSpace(uB, &mesh, elementType::VECTOR_LAGRANGE, orderU, "U",          "Bord", dQ, &uSol));
      feCheck(createFiniteElementSpace( p, &mesh, elementType::LAGRANGE,      orderU-1, "P",       "Domaine", dQ, &pSol));
      feCheck(createFiniteElementSpace(pP, &mesh, elementType::LAGRANGE,      orderU-1, "P", "PointPression", dQ, &pSol));

      feCheck(createFiniteElementSpace( phi, &mesh, elementType::LAGRANGE, orderPhiMu, "Phi", "Domaine", dQ, &phiSol));
      feCheck(createFiniteElementSpace(phiB, &mesh, elementType::LAGRANGE, orderPhiMu, "Phi",    "Bord", dQ, &phiSol));
      feCheck(createFiniteElementSpace(  mu, &mesh, elementType::LAGRANGE, orderPhiMu,  "Mu", "Domaine", dQ, &muSol));
      feCheck(createFiniteElementSpace( muB, &mesh, elementType::LAGRANGE, orderPhiMu,  "Mu",    "Bord", dQ, &muSol));

      std::vector<feSpace*> spaces = {u, uB, p, pP, phi, phiB, mu, muB};
      std::vector<feSpace*> essentialSpaces = {uB, pP, phiB, muB};
      feMetaNumber numbering(&mesh, spaces, essentialSpaces);
      feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

      // Momentum
      feBilinearForm *CHNS_momentum = nullptr, *uS = nullptr;
      feCheck(createBilinearForm(CHNS_momentum, {u, p, phi, mu}, new feSysElm_CHNS_Momentum<2>(&density, &drhodphi, &viscosity, &dviscdphi, &mobility, &coeffKorteweg, &volumeForce)));
      feCheck(createBilinearForm(           uS,             {u}, new feSysElm_VectorSource<2>(&uSource)));

      // Continuity
      feBilinearForm *divU = nullptr;
      feCheck(createBilinearForm(divU, {p, u}, new feSysElm_MixedDivergence<2>(&scalarConstant::one)));

      // Tracer convection
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

      std::vector<feBilinearForm*> forms = {CHNS_momentum, uS, divU, dphidt,
        convPhi, diffMu, phiS, muMass, phiCube, phiMass, diffPhi, muS};

      feLinearSystem *system;
      #if defined(HAVE_MKL)
        feCheck(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
        feCheck(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
      #else
        feCheck(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
      #endif

      std::vector<feNorm *> norms = {};

      double t0 = 0.;
      // double t1 = 0.1;
      double t1 = 0.0001;
      int nTimeSteps = 10;
      TimeIntegrator *solver;
      feCheck(createTimeIntegrator(solver, timeIntegratorScheme::BDF2, NLoptions,
        system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));

      // Start BDF2 with manufactured solution
      static_cast<BDF2Integrator*>(solver)->setStartingMethod(BDF2Starter::InitialCondition);

      // Solve
      feCheck(solver->makeSteps(nTimeSteps));

      feNorm *errorU_L2 = nullptr, *errorP_L2 = nullptr, *errorPhi_L2 = nullptr, *errorMu_L2 = nullptr;
      feCheck(createNorm(errorU_L2, VECTOR_L2_ERROR, {u}, &sol, nullptr, &uSol));
      feCheck(createNorm(errorP_L2, L2_ERROR, {p}, &sol, &pSol));
      feCheck(createNorm(errorPhi_L2, L2_ERROR, {phi}, &sol, &phiSol));
      feCheck(createNorm(errorMu_L2 , L2_ERROR, {mu}, &sol, &muSol));
      L2ErrorU   = errorU_L2->compute();
      L2ErrorP   = errorP_L2->compute();
      L2ErrorPhi = errorPhi_L2->compute();
      L2ErrorMu  = errorMu_L2->compute();

      // feNorm *errorU_H1 = nullptr, *errorP_H1 = nullptr, *errorPhi_H1 = nullptr, *errorMu_H1 = nullptr;
      // feCheck(createNorm(errorU_H1, VECTOR_H1_ERROR, {u}, &sol, nullptr, &uSol));
      // feCheck(createNorm(errorP_H1, H1_ERROR, {p}, &sol, &pSol));
      // feCheck(createNorm(errorPhi_H1, H1_ERROR, {phi}, &sol, &phiSol));
      // feCheck(createNorm(errorMu_H1 , H1_ERROR, {mu}, &sol, &muSol));
      // H1ErrorU   = errorU_H1->compute();
      // H1ErrorP   = errorP_H1->compute();
      // H1ErrorPhi = errorPhi_H1->compute();
      // H1ErrorMu  = errorMu_H1->compute();

      delete errorU_L2;
      delete errorP_L2;
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

    int meshConvergence(std::stringstream &resultBuffer, int orderU, int orderPhiMu, int numMeshes, int degreeQuadrature)
    {
      feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};

      std::vector<int> nElm(numMeshes);
      std::vector<int> nVertices(numMeshes);
      std::vector<double> errU(numMeshes, 0.);
      std::vector<double> errP(numMeshes, 0.);
      std::vector<double> errPhi(numMeshes, 0.);
      std::vector<double> errMu(numMeshes, 0.);
      // std::vector<double> errU_H1(numMeshes, 0.);
      // std::vector<double> errP_H1(numMeshes, 0.);
      // std::vector<double> errPhi_H1(numMeshes, 0.);
      // std::vector<double> errMu_H1(numMeshes, 0.);

      bool writeErrorToFile = false;
      std::string outputDirectory = "./";
      std::ofstream errorFileL2, errorFileH1;
      if(writeErrorToFile) {
        errorFileL2.open(outputDirectory + "errorL2.txt");
        errorFileH1.open(outputDirectory + "errorH1.txt");
      }

      for(int i = 0; i < numMeshes; ++i)
      {
        std::string meshFile = "../../../data/stokes" + std::to_string(i+1) + ".msh";
        cahnHilliardNavierStokes::MMS::solve(meshFile, orderU, orderPhiMu,
          degreeQuadrature, NLoptions, nElm[i], nVertices[i],
          errU[i], errP[i], errPhi[i], errMu[i]);
          // errU_H1[i], errP_H1[i], errPhi_H1[i], errMu_H1[i]);

        if(writeErrorToFile) {
          errorFileL2 << nVertices[i] << "\t" << nElm[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errU[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errP[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errPhi[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errMu[i] << "\n";
          // errorFileL2.flush();
        }
      }

      if(writeErrorToFile) {
        errorFileL2.close();
        errorFileH1.close();
      }

      resultBuffer << "   (u,p) : Taylor-Hood elements P" << orderU << "-P" << orderU-1 << std::endl; 
      resultBuffer << "(phi,mu) : Equal order elements P" << orderPhiMu << "-P" << orderPhiMu << std::endl; 
      resultBuffer << "Cahn-Hilliard Navier-Stokes MMS - Error on velocity" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errU, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      resultBuffer << "Cahn-Hilliard Navier-Stokes MMS - Error on pressure" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errP, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      resultBuffer << "Cahn-Hilliard Navier-Stokes MMS - Error on phase marker phi" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errPhi, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      resultBuffer << "Cahn-Hilliard Navier-Stokes MMS - Error on chemical potential mu" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errMu, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

      return 0;
    }

  } // namespace MMS
}


TEST(CahnHilliard, CH)
{
  initialize(my_argc, my_argv);
  setVerbose(0);
  
  std::string testRoot = "../../../tests/withLinearSolver/CH_MMS";

  std::stringstream resultBuffer;

  int degreeQuadrature = 8;
  int order = 1;
  cahnHilliard::MMS::meshConvergence(resultBuffer, order, 4, degreeQuadrature);
  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

TEST(CahnHilliard, CHNS)
{
  initialize(my_argc, my_argv);
  setVerbose(1);
  
  std::string testRoot = "../../../tests/withLinearSolver/CHNS_MMS";

  std::stringstream resultBuffer;

  int degreeQuadrature = 8;
  int orderU = 2;
  int orderPhiMu = 1;
  cahnHilliardNavierStokes::MMS::meshConvergence(resultBuffer, orderU, orderPhiMu, 4, degreeQuadrature);
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