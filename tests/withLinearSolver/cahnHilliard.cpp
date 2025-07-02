
#include "feAPI.h"
#include "feParameters.h"
#include "CHNS_Solver.h"
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

      UNUSED(gamma, epsilon);
      const double coeffKorteweg = par[8];

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
      // const double korteweg = - gamma/epsilon;
      const double dviscdphi = (viscA-viscB)/2.;

      const double d[2][2] = {{gradu[0][0], (gradu[0][1] + gradu[1][0])/2.},
                        {(gradu[0][1] + gradu[1][0])/2., gradu[1][1]}};
      const double gradPhidotd[2] = {gradphi[0]*d[0][0] + gradphi[1]*d[1][0],
                               gradphi[0]*d[0][1] + gradphi[1]*d[1][1]};

      if(isStationary) {
        // For stationary test, do not include dudt because it will not be computed by the solver (BDF coefficient is 0).
        // If dudt is included, it adds du/dt at t = 0 to the source term, which is wrong with a stationary solver.
        UNUSED(dudt);
        res[0] = - (rho * (          uDotGradu[0] - f[0]) + M * drhodphi * gMudgu[0] + gradp[0] - 2.*(dviscdphi*gradPhidotd[0] + visc*divd[0]) + coeffKorteweg*mu*gradphi[0]);
        res[1] = - (rho * (          uDotGradu[1] - f[1]) + M * drhodphi * gMudgu[1] + gradp[1] - 2.*(dviscdphi*gradPhidotd[1] + visc*divd[1]) + coeffKorteweg*mu*gradphi[1]);
      } else {
        res[0] = - (rho * (dudt[0] + uDotGradu[0] - f[0]) + M * drhodphi * gMudgu[0] + gradp[0] - 2.*(dviscdphi*gradPhidotd[0] + visc*divd[0]) + coeffKorteweg*mu*gradphi[0]);
        res[1] = - (rho * (dudt[1] + uDotGradu[1] - f[1]) + M * drhodphi * gMudgu[1] + gradp[1] - 2.*(dviscdphi*gradPhidotd[1] + visc*divd[1]) + coeffKorteweg*mu*gradphi[1]);
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
      // const double epsilon = par[0];
      const double coeffPhiCube = par[0];
      const double coeffPhi     = par[1];
      const double coeffDiffPhi = par[2];
      const double sx = sin(M_PI*x);
      const double sy = sin(M_PI*y);
      const double emt = exp(-t);
      const double mu = sx*sy;
      const double phi     = emt * sx*sy;
      const double phiCube = emt*emt*emt*sx*sx*sx*sy*sy*sy;
      const double lapPhi  = emt * (-2.*M_PI*M_PI*sx*sy);

      return - (mu + coeffPhiCube * phiCube + coeffPhi * phi + coeffDiffPhi * lapPhi);
    }

    int solve(const std::string &meshFile,
              const int timeStepFactor,
              const int orderU,
              const int orderPhiMu,
              const int degreeQuadrature,
              const feNLSolverOptions &NLoptions,
              int &numInteriorElements,
              int &numVertices,
              double &L2ErrorU,
              double &L2ErrorP,
              double &L2ErrorPhi,
              double &L2ErrorMu,
              int &timesteps)
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

      // // coeffKorteweg = new feConstantFunction(-_gamma/_epsilon);
      // // coeffPhiCube  = new feConstantFunction(-1.);
      // // coeffPhi      = new feConstantFunction(+1.);
      // // coeffDiffPhi  = new feConstantFunction(-_epsilon*_epsilon);

      // // const double lambda = 3. * _epsilon * _gamma / (2.*sqrt(2.));
      // const double lambda = _epsilon * _gamma;
      // coeffKorteweg = new feConstantFunction(-1.);
      // coeffPhiCube  = new feConstantFunction(- lambda / (_epsilon * _epsilon));
      // coeffPhi      = new feConstantFunction(+ lambda / (_epsilon * _epsilon));
      // coeffDiffPhi  = new feConstantFunction(- lambda);

      // // Mémoire de Zeineb:
      // const double kortewegCoeff = -gamma/epsilon;
      // const double phiCubeCoeff  = -1.;
      // const double phiCoeff      = +1.;
      // const double diffPhiCoeff  = -epsilon*epsilon;

      // Mémoire de Pierre:
      const double lambda = 3. * epsilon * gamma / (2. * sqrt(2.));
      const double kortewegCoeff = -1.;
      // const double phiCubeCoeff  = - (epsilon * epsilon) / lambda;
      // const double phiCoeff      = + (epsilon * epsilon) / lambda;
      const double phiCubeCoeff  = - lambda / (epsilon * epsilon);
      const double phiCoeff      = + lambda / (epsilon * epsilon);
      const double diffPhiCoeff  = - lambda;

      feConstantFunction coeffKorteweg(kortewegCoeff);
      feConstantFunction coeffPhiCube(phiCubeCoeff);
      feConstantFunction coeffPhi(phiCoeff);
      feConstantFunction coeffDiffPhi(diffPhiCoeff);

      feVectorFunction uSource(uSource_f, {rhoA, rhoB, viscA, viscB, mobility_val, gamma, epsilon, stationaryVal, kortewegCoeff});
      feFunction pSource(pSource_f, {rhoA, rhoB});
      feFunction phiSource(phiSource_f, {stationaryVal, mobility_val});
      feFunction  muSource( muSource_f, {phiCubeCoeff, phiCoeff, -diffPhiCoeff}); // Minus sign because the coeff is set after integration by parts

      feConstantFunction mobility(mobility_val);
      feFunction density(rho_f, {rhoA, rhoB});
      feFunction drhodphi(drhodphi_f, {rhoA, rhoB});
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

      // SUPG for tracer
      // feBilinearForm *supgTracer = nullptr;
      // feCheckReturn(createBilinearForm(supgTracer, {phi, u, mu}, new feSysElm_CHNS_Tracer_SUPG<2>(&scalarConstant::one, &mobility)));

      // Chemical potential
      feBilinearForm *muMass = nullptr, *phiCube = nullptr, *phiMass = nullptr, *diffPhi = nullptr, *muS = nullptr;
      feCheck(createBilinearForm( muMass,      {mu}, new feSysElm_Mass(&scalarConstant::one)));
      feCheck(createBilinearForm(phiCube, {mu, phi}, new feSysElm_MixedMassPower(&coeffPhiCube, 3)));
      feCheck(createBilinearForm(phiMass, {mu, phi}, new feSysElm_MixedMass(&coeffPhi)));
      feCheck(createBilinearForm(diffPhi, {mu, phi}, new feSysElm_MixedGradGrad<2>(&coeffDiffPhi)));
      feCheck(createBilinearForm(    muS,      {mu}, new feSysElm_Source(&muSource)));

      std::vector<feBilinearForm*> forms = {CHNS_momentum, uS, divU, dphidt,
        convPhi, diffMu, phiS, muMass, phiCube, phiMass, diffPhi, muS,
        // supgTracer
      };

      // CHNS_momentum->setComputeMatrixWithFD(true);

      feLinearSystem *system;
      #if defined(HAVE_MKL)
        feCheck(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
        feCheck(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
      #else
        feCheck(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
      #endif

      feNorm *errorU_L2 = nullptr, *errorP_L2 = nullptr, *errorPhi_L2 = nullptr, *errorMu_L2 = nullptr;
      feCheck(createNorm(errorU_L2, VECTOR_L2_ERROR, {u}, &sol, nullptr, &uSol));
      feCheck(createNorm(errorP_L2, L2_ERROR, {p}, &sol, &pSol));
      feCheck(createNorm(errorPhi_L2, L2_ERROR, {phi}, &sol, &phiSol));
      feCheck(createNorm(errorMu_L2 , L2_ERROR, {mu}, &sol, &muSol));
      // feCheck(createNorm(errorU_L2, VECTOR_LINF_ERROR, {u}, &sol, nullptr, &uSol));
      // feCheck(createNorm(errorP_L2, LINF_ERROR, {p}, &sol, &pSol));
      // feCheck(createNorm(errorPhi_L2, LINF_ERROR, {phi}, &sol, &phiSol));
      // feCheck(createNorm(errorMu_L2 , LINF_ERROR, {mu}, &sol, &muSol));
      std::vector<feNorm *> norms = {errorU_L2, errorP_L2, errorPhi_L2, errorMu_L2};

      // Reference test file matches those parameters and L2 error evaluated at t1:
      UNUSED(timeStepFactor);
      double t0 = 0.;
      double t1 = 0.0001;
      int nTimeSteps = 10;

      // For other convergence studies:
      // double t1 = 0.00001;
      // double t1 = 5.;
      // int nTimeSteps = 10 * pow(2, timeStepFactor);

      timesteps = nTimeSteps;
      TimeIntegrator *solver;
      feCheck(createTimeIntegrator(solver, timeIntegratorScheme::BDF2, NLoptions,
        system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));

      // const double dt = solver->getTimeStep();

      // Start BDF2 with manufactured solution
      static_cast<BDF2Integrator*>(solver)->setStartingMethod(BDF2Starter::InitialCondition);

      // Solve
      feCheck(solver->makeSteps(nTimeSteps));

      // Use L2 at t1 to compare with reference test files:
      L2ErrorU   = errorU_L2->compute();
      L2ErrorP   = errorP_L2->compute();
      L2ErrorPhi = errorPhi_L2->compute();
      L2ErrorMu  = errorMu_L2->compute();

      // Or use L1 or Linf error in time:
      // std::vector<std::vector<double>> postProc = solver->getPostProcessingData();
      // L2ErrorU = 0.;
      // L2ErrorP = 0.;
      // L2ErrorPhi = 0.;
      // L2ErrorMu = 0.;
      // for(int j = 1; j < nTimeSteps+1; ++j) {

      //   // L1 norm in time, method of rectangles
      //   L2ErrorU   += dt * postProc[1][j];
      //   L2ErrorP   += dt * postProc[2][j];
      //   L2ErrorPhi += dt * postProc[3][j];
      //   L2ErrorMu  += dt * postProc[4][j];

      //   // // Linf norm in time
      //   // L2ErrorU = fmax(L2ErrorU, postProc[1][j]);
      //   // L2ErrorP = fmax(L2ErrorP, postProc[2][j]);
      //   // L2ErrorPhi = fmax(L2ErrorPhi, postProc[3][j]);
      //   // L2ErrorMu = fmax(L2ErrorMu, postProc[4][j]);
      // }

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
      std::vector<int> timesteps(numMeshes);
      std::vector<int> nVertices(numMeshes);
      std::vector<double> errU(numMeshes, 0.);
      std::vector<double> errP(numMeshes, 0.);
      std::vector<double> errPhi(numMeshes, 0.);
      std::vector<double> errMu(numMeshes, 0.);
      // std::vector<double> errU_H1(numMeshes, 0.);
      // std::vector<double> errP_H1(numMeshes, 0.);
      // std::vector<double> errPhi_H1(numMeshes, 0.);
      // std::vector<double> errMu_H1(numMeshes, 0.);

      bool writeErrorToFile = true;
      std::string outputDirectory = "./";
      std::ofstream errorFileL2, errorFileH1, timeErrorFile;
      if(writeErrorToFile) {
        errorFileL2.open(outputDirectory + "errorL2.txt");
        errorFileH1.open(outputDirectory + "errorH1.txt");
        timeErrorFile.open(outputDirectory + "timeErrorFile.txt");
      }

      for(int i = 0; i < numMeshes; ++i)
      {
        const int timeStepFactor = i;
        std::string meshFile = "../../../data/stokes" + std::to_string(i+1) + ".msh";
        // std::string meshFile = "../../../data/stokes4.msh";
        cahnHilliardNavierStokes::MMS::solve(meshFile, timeStepFactor, orderU, orderPhiMu,
          degreeQuadrature, NLoptions, nElm[i], nVertices[i],
          errU[i], errP[i], errPhi[i], errMu[i], timesteps[i]);
          // errU_H1[i], errP_H1[i], errPhi_H1[i], errMu_H1[i]);

        if(writeErrorToFile) {
          errorFileL2 << nVertices[i] << "\t" << nElm[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errU[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errP[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errPhi[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errMu[i] << "\n";
          errorFileL2.flush();
          timeErrorFile << timesteps[i] << "\t" << nElm[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errU[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errP[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errPhi[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errMu[i] << "\n";
          timeErrorFile.flush();
        }
      }

      if(writeErrorToFile) {
        errorFileL2.close();
        errorFileH1.close();
        timeErrorFile.close();
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

      // resultBuffer << "Time convergence:" << std::endl; 
      // computeAndPrintConvergence(1, numMeshes, errU, timesteps, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      // computeAndPrintConvergence(1, numMeshes, errP, timesteps, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      // computeAndPrintConvergence(1, numMeshes, errPhi, timesteps, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      // computeAndPrintConvergence(1, numMeshes, errMu, timesteps, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

      return 0;
    }

  } // namespace MMS
}

///////////////////////////////////////////////////////////////////////////////
// Modify here for another MMS

// #define TEST_TIME_CONVERGENCE

#define CONSTANT_MOBILITY(M, phi) M
#define DCONSTANT_MOBILITY_DPHI(M, phi) 0. 

// Time dependency.
// The u dot grad(phi) term in the tracer equation is of order FT*FT,
// so keep this in mind when choosing a time integrator.
// E.g., FT = t yields t^2 time dependency, which should be captured by a BDF2.
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

  #define DUDX(x,y,t) FT(t) * M_PI *   cos(M_PI*x) * sin(M_PI*y)
  #define DUDY(x,y,t) FT(t) * M_PI *   sin(M_PI*x) * cos(M_PI*y)
  #define DVDX(x,y,t) FT(t) * M_PI * - sin(M_PI*x) * cos(M_PI*y)
  #define DVDY(x,y,t) FT(t) * M_PI * - cos(M_PI*x) * sin(M_PI*y)

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
  #define DPDYY(x,y,t)   M_PI * M_PI * sin(M_PI*x) * cos(M_PI*y)

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

//
// Parameters for volume-averaged formulation of CHNS
//
#define DECLARE_VARIABLES(args, par)    \
  const double x   = args.pos[0];       \
  const double y   = args.pos[1];       \
  const double t   = args.t;            \
  const double phi = args.u;            \
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
  UNUSED(x,y,t,phi,rhoA, rhoB,          \
    viscA, viscB, M, gamma,             \
    epsilon, lambda);                   \


struct CahnHilliardVolumeAveraged {

  static void uSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);
    res[0] = U(x,y,t);
    res[1] = V(x,y,t);
  }

  static void grad_uSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);
    const double gradu[2][2] = GRADU(x,y,t);
    res[0] = gradu[0][0];
    res[1] = gradu[0][1];
    res[2] = gradu[1][0];
    res[3] = gradu[1][1];
  }

  static double pSol_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return P(x,y,t);
  }

  static void grad_pSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);
    const double gradp[2] = GRADP(x,y,t);
    res[0] = gradp[0];
    res[1] = gradp[1];
  }

  static double phiSol_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return PHI(x,y,t);
  }

  static void grad_phiSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);
    const double gradphi[2] = GRADPHI(x,y,t);
    res[0] = gradphi[0];
    res[1] = gradphi[1];
  }

  static double muSol_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return MU(x,y,t);
  }

  static void grad_muSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);
    const double gradmu[2] = GRADMU(x,y,t);
    res[0] = gradmu[0];
    res[1] = gradmu[1];
  }

  //
  // CHNS models for density, viscosity and mobility
  //
  static double rho_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return (rhoA - rhoB)/2. * phi + (rhoA + rhoB)/2.;
  }

  static double drhodphi_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return (rhoA - rhoB)/2.;
  }

  static double viscosity_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return (viscA - viscB)/2. * phi + (viscA + viscB)/2.;
  }

  static double dviscdphi_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return (viscA - viscB)/2.;
  }

  static double constant_mobility_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return CONSTANT_MOBILITY(M, phi);
  }

  //
  // Degenerate mobility (0 in pure phases)
  // Here we take (1 - phi^2) as the dependence in phi,
  // that is:
  //
  // M(phi) = M * (1 - phi^2).
  //
  // Since phi is a product of sines, phi is in [-1,1] and M(phi) >= 0.
  //
  // static double degenerate_mobility_f(const feFunctionArguments &args, const std::vector<double> &par)
  // {
  //   DECLARE_VARIABLES(args, par);
  //   return MOBILITY(M, phi);
  // }

  //
  // Source terms
  //
  static double pSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    const double div_u = DIVU(x,y,t);
    return - (div_u);
  }

  static void uSource_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);

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

  static double phiSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
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
    const double mu      = MU(x,y,t);
    const double phiRef  = PHI(x,y,t);
    const double phiCube = PHI(x,y,t) * PHI(x,y,t) * PHI(x,y,t);
    const double lapPhi  = LAP_PHI(x,y,t);

    return - (mu - lambda / (epsilon*epsilon) * (phiCube - phiRef) + lambda * lapPhi);
  }

  // CahnHilliardVolumeAveraged
  static feStatus solve(const std::string &meshFile,
                        const int timeStepFactor,
                        const int orderU,
                        const int orderPhiMu,
                        const int degreeQuadrature,
                        const feNLSolverOptions &NLoptions,
                        int &numInteriorElements,
                        int &numVertices,
                        double &L2ErrorU,
                        double &L2ErrorP,
                        double &L2ErrorPhi,
                        double &L2ErrorMu,
                        double &integralP,
                        double &H1ErrorU,
                        double &H1ErrorP,
                        double &H1ErrorPhi,
                        double &H1ErrorMu,
                        int &timesteps)
  {
    const double rhoA    = 1000.;
    const double rhoB    = 1.;
    const double viscA   = 100.;
    const double viscB   = 1.;
    const double M       = 1e-6;
    const double gamma   = 0.01;
    const double epsilon = 1e-3;

    const double lambda = 3. / (2.*sqrt(2.)) * epsilon * gamma;

    std::vector<double> CHNSparameters = {gamma, epsilon};
    std::vector<double> params = {rhoA, rhoB, viscA, viscB, M, gamma, epsilon, lambda};

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

    feFunction   density(      rho_f, params);
    feFunction  drhodphi( drhodphi_f, params);
    feFunction viscosity(viscosity_f, params);
    feFunction dviscdphi(dviscdphi_f, params);

    feFunction mobility(constant_mobility_f, params);

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

    feBilinearForm *CHNS;
    feCheckReturn(createBilinearForm(CHNS, {u, p, phi, mu}, new CHNS_VolumeAveraged<2>(
      &density,
      &drhodphi,
      &viscosity,
      &dviscdphi,
      &mobility,
      &volumeForce,
      &pSource,
      &uSource,
      &phiSource,
      &muSource,
      CHNSparameters)));
    std::vector<feBilinearForm*> forms = {CHNS};

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      feCheck(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheck(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
    #else
      feCheck(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
    #endif

    feNorm *errorU_L2 = nullptr, *errorP_L2 = nullptr, *errorPhi_L2 = nullptr, *errorMu_L2 = nullptr;
    feNorm *errorU_H1 = nullptr, *errorP_H1 = nullptr, *errorPhi_H1 = nullptr, *errorMu_H1 = nullptr;

    feCheckReturn(createNorm(  errorU_L2, VECTOR_L2_ERROR,   {u}, &sol, nullptr, &uSol));
    feCheckReturn(createNorm(  errorP_L2,        L2_ERROR,   {p}, &sol, &pSol));
    feCheckReturn(createNorm(errorPhi_L2,        L2_ERROR, {phi}, &sol, &phiSol));
    feCheckReturn(createNorm( errorMu_L2,        L2_ERROR,  {mu}, &sol, &muSol));

    feCheckReturn(createNorm(  errorU_H1, VECTOR_SEMI_H1_ERROR,   {u}, &sol, nullptr, &grad_uSol));
    feCheckReturn(createNorm(  errorP_H1,        SEMI_H1_ERROR,   {p}, &sol, nullptr, &grad_pSol));
    feCheckReturn(createNorm(errorPhi_H1,        SEMI_H1_ERROR, {phi}, &sol, nullptr, &grad_phiSol));
    feCheckReturn(createNorm( errorMu_H1,        SEMI_H1_ERROR,  {mu}, &sol, nullptr, &grad_muSol));

    std::vector<feNorm *> norms = {};

    // Reference test file matches those parameters and L2 error evaluated at t1:
    UNUSED(timeStepFactor);
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 20;
    timesteps = nTimeSteps;

    TimeIntegrator *solver;
    timeIntegratorScheme scheme = timeIntegratorScheme::BDF2;
    feCheck(createTimeIntegrator(solver, scheme, NLoptions,
      system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));

    // Start BDF2 with manufactured solution
    if(scheme == timeIntegratorScheme::BDF2) {
      static_cast<BDF2Integrator*>(solver)->setStartingMethod(BDF2Starter::InitialCondition);
    }

    // Solve
    feCheck(solver->makeSteps(nTimeSteps));

    // Use L2 at t1 to compare with reference test files:
    L2ErrorU   = errorU_L2->compute();
    L2ErrorP   = errorP_L2->compute();
    L2ErrorPhi = errorPhi_L2->compute();
    L2ErrorMu  = errorMu_L2->compute();

    // Use L2 at t1 to compare with reference test files:
    L2ErrorU   = errorU_L2->compute();
    L2ErrorP   = errorP_L2->compute();
    L2ErrorPhi = errorPhi_L2->compute();
    L2ErrorMu  = errorMu_L2->compute();

    H1ErrorU   = errorU_H1->compute();
    H1ErrorP   = errorP_H1->compute();
    H1ErrorPhi = errorPhi_H1->compute();
    H1ErrorMu  = errorMu_H1->compute();

    // Integral of pressure (target is 0)
    feNorm *intP;
    feCheckReturn(createNorm(intP, INTEGRAL, {p}, &sol));
    integralP = fabs(intP->compute());
    delete intP;

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

    return FE_STATUS_OK;
  }

  static feStatus ezsolve(const std::string &meshFile,
                        const int timeStepFactor,
                        const int orderU,
                        const int orderPhiMu,
                        const int degreeQuadrature,
                        const feNLSolverOptions &NLoptions,
                        int &numInteriorElements,
                        int &numVertices,
                        double &L2ErrorU,
                        double &L2ErrorP,
                        double &L2ErrorPhi,
                        double &L2ErrorMu,
                        double &integralP,
                        double &H1ErrorU,
                        double &H1ErrorP,
                        double &H1ErrorPhi,
                        double &H1ErrorMu,
                        int &timesteps)
  {
    UNUSED(NLoptions);

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
    CHNSparameters.formulation = Parameters::CHNS::Formulation::volumeAveraged;
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

    VectorFEDescriptor   u("U", "Domaine", vlag, orderU, degreeQuadrature, &uSol);
    ScalarFEDescriptor   p("P", "Domaine", lag, orderU-1, degreeQuadrature, &pSol);
    ScalarFEDescriptor phi("Phi", "Domaine", lag, orderPhiMu, degreeQuadrature, &phiSol);
    ScalarFEDescriptor  mu("Mu", "Domaine", lag, orderPhiMu, degreeQuadrature, &muSol);

    std::vector<FEDescriptor> fieldDescriptors = {u, p, phi, mu};

    // The exact solutions and gradients for those fields, and the associated norms to compute
    std::vector<std::pair<int, const feFunction *>> exactScalarSol = {{1, &pSol}, {2, &phiSol}, {3, &muSol}};
    std::vector<std::pair<int, const feVectorFunction *>> exactVectorSol = {{0, &grad_uSol}, {1, &grad_pSol}, {2, &grad_phiSol}, {3, &grad_muSol}};
    std::vector<std::pair<int, normType>> normTypes = {

      {1, L2_ERROR},
      {2, L2_ERROR},
      {3, L2_ERROR},

      {0, VECTOR_SEMI_H1_ERROR},
      {1, SEMI_H1_ERROR},
      {2, SEMI_H1_ERROR},
      {3, SEMI_H1_ERROR}};

    VectorFEDescriptor uB("U", "Bord", vlag, orderU, degreeQuadrature, &uSol);
    ScalarFEDescriptor pP("P", "PointPression", lag, orderU-1, degreeQuadrature, &pSol);
    ScalarFEDescriptor phiB("Phi", "Bord", lag, orderPhiMu, degreeQuadrature, &phiSol);
    ScalarFEDescriptor muB("Mu", "Bord", lag, orderPhiMu, degreeQuadrature, &muSol);

    BoundaryConditions::VectorDirichlet uBC(uB);
    BoundaryConditions::ScalarDirichlet pBC(pP), phiBC(phiB), muBC(muB);

    std::vector<const BoundaryConditions::BoundaryCondition*> BC = {&uBC, &pBC, &phiBC, &muBC};

    Parameters::NonLinearSolver NLSolver_param;

    Parameters::TimeIntegration TimeIntegration_param;
    TimeIntegration_param.t_initial  = 0.;
    TimeIntegration_param.t_final    = 1.;

    #if defined(TEST_TIME_CONVERGENCE)
      TimeIntegration_param.nTimeSteps = 10 * pow(2, timeStepFactor);
    #else 
      TimeIntegration_param.nTimeSteps = 10 * pow(2, timeStepFactor);
    #endif

    timesteps = TimeIntegration_param.nTimeSteps;

    TimeIntegration_param.method =
      Parameters::TimeIntegration::TimeIntegrationMethod::bdf2;
    TimeIntegration_param.bdf2starter =
      // Parameters::TimeIntegration::BDF2Starter::bdf1;
      Parameters::TimeIntegration::BDF2Starter::exactSolution;


    CHNS_Solver *solver;
    feCheck(createCHNS_Solver(solver,
                              CHNSparameters,
                              NLSolver_param,
                              TimeIntegration_param,
                              fieldDescriptors,
                              BC,
                              nullptr,
                              &pSource,
                              &uSource,
                              &phiSource,
                              &muSource));

    solver->setScalarSolutions(exactScalarSol);
    solver->setVectorSolutions(exactVectorSol);
    solver->setNormTypes(normTypes);

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();
    numVertices = mesh.getNumVertices();

    std::vector<feSpace *> spaces, essentialSpaces;
    solver->createSpaces(&mesh, spaces, essentialSpaces);

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution   sol(numbering.getNbDOFs(), spaces, essentialSpaces);

    feExporter *exporter = nullptr;
    feExportData exportData = {exporter, 1, "sol"};

    feCheck(solver->solve(
      &mesh, &sol, &numbering, spaces, exportData));

    std::vector<std::vector<double>> postProcData = solver->getPostProcessingData();

    L2ErrorU = 0.;
    L2ErrorP = 0.;
    L2ErrorPhi = 0.;
    L2ErrorMu = 0.;

    H1ErrorU = 0.;
    H1ErrorP = 0.;
    H1ErrorPhi = 0.;
    H1ErrorMu = 0.;

    // First post proc vector is the time
    const double dt = postProcData[0][1] - postProcData[0][0];

    // Compute L1 norm in time of the computed spatial errors
    for(size_t i = 0; i < postProcData[0].size(); ++i)
    {
      // L2ErrorU   += dt * postProcData[1][i];
      L2ErrorP   += dt * postProcData[1][i];
      L2ErrorPhi += dt * postProcData[2][i];
      L2ErrorMu  += dt * postProcData[3][i];

      H1ErrorU   += dt * postProcData[4][i];
      H1ErrorP   += dt * postProcData[5][i];
      H1ErrorPhi += dt * postProcData[6][i];
      H1ErrorMu  += dt * postProcData[7][i];
    }

    feNorm *errorU_L2 = nullptr;
    // feNorm *errorP_L2 = nullptr, *errorPhi_L2 = nullptr, *errorMu_L2 = nullptr;
    // feNorm *errorU_H1 = nullptr, *errorP_H1 = nullptr, *errorPhi_H1 = nullptr, *errorMu_H1 = nullptr;

    feSpace *uS = spaces[0];
    // feSpace *pS = spaces[1];
    // feSpace *phiS = spaces[2];
    // feSpace *muS = spaces[3];
    feCheckReturn(createNorm(  errorU_L2, VECTOR_L2_ERROR,   {uS}, &sol, nullptr, &uSol));
    // feCheckReturn(createNorm(  errorP_L2,        L2_ERROR,   {pS}, &sol, &pSol));
    // feCheckReturn(createNorm(errorPhi_L2,        L2_ERROR, {phiS}, &sol, &phiSol));
    // feCheckReturn(createNorm( errorMu_L2,        L2_ERROR,  {muS}, &sol, &muSol));

    // feCheckReturn(createNorm(  errorU_H1, VECTOR_SEMI_H1_ERROR,   {uS}, &sol, nullptr, &grad_uSol));
    // feCheckReturn(createNorm(  errorP_H1,        SEMI_H1_ERROR,   {pS}, &sol, nullptr, &grad_pSol));
    // feCheckReturn(createNorm(errorPhi_H1,        SEMI_H1_ERROR, {phiS}, &sol, nullptr, &grad_phiSol));
    // feCheckReturn(createNorm( errorMu_H1,        SEMI_H1_ERROR,  {muS}, &sol, nullptr, &grad_muSol));

    // Use L2 at t1 to compare with reference test files:
    L2ErrorU   = errorU_L2->compute();
    // L2ErrorP   = errorP_L2->compute();
    // L2ErrorPhi = errorPhi_L2->compute();
    // L2ErrorMu  = errorMu_L2->compute();

    // H1ErrorU   = errorU_H1->compute();
    // H1ErrorP   = errorP_H1->compute();
    // H1ErrorPhi = errorPhi_H1->compute();
    // H1ErrorMu  = errorMu_H1->compute();

    // Integral of pressure (target is 0)
    // feNorm *intP;
    // feCheckReturn(createNorm(intP, INTEGRAL, {pS}, &sol));
    // integralP = fabs(intP->compute());
    // delete intP;
    integralP = 0.;

    // delete errorU_L2;
    // delete errorP_L2;
    // delete errorPhi_L2;
    // delete errorMu_L2;

    for(feSpace *s : spaces)
      delete s;
    // delete solver;

    return FE_STATUS_OK;
  }
};

#undef DECLARE_VARIABLES

#define DECLARE_VARIABLES(args, par)    \
  const double x   = args.pos[0];       \
  const double y   = args.pos[1];       \
  const double t   = args.t;            \
  const double phi = args.u;            \
                                        \
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
  UNUSED(x,y,t,phi,rhoA, rhoB,          \
    viscA, viscB, M, gamma,             \
    epsilon, beta, tau, alpha);         \

//
// Alternative mass-conservative formulation of CHNS
//
namespace alternativeCHNS_constantMobility {

  namespace MMS {

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

    double phiSol_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      DECLARE_VARIABLES(args, par);
      return PHI(x,y,t);
    }

    double muSol_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      DECLARE_VARIABLES(args, par);
      return MU(x,y,t);
    }

    double rho_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      DECLARE_VARIABLES(args, par);
      return (rhoA - rhoB)/2. * phi + (rhoA + rhoB)/2.;
    }

    double drhodphi_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      DECLARE_VARIABLES(args, par);
      return (rhoA - rhoB)/2.;
    }

    double viscosity_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      DECLARE_VARIABLES(args, par);
      return (viscA - viscB)/2. * phi + (viscA + viscB)/2.;
    }

    double dviscdphi_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      DECLARE_VARIABLES(args, par);
      return (viscA - viscB)/2.;
    }

    double mobility_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      DECLARE_VARIABLES(args, par);
      return CONSTANT_MOBILITY(M, phi);
    }

    double mobilityAlpha_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      DECLARE_VARIABLES(args, par);
      return alpha * CONSTANT_MOBILITY(M, phi);
    }

    double mobilityAlphaSquared_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      DECLARE_VARIABLES(args, par);
      return alpha * alpha * CONSTANT_MOBILITY(M, phi);
    }

    double pSource_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      DECLARE_VARIABLES(args, par);

      const double div_u = DIVU(x,y,t);
      const double lapP  = LAP_P(x,y,t);
      const double lapMu = LAP_MU(x,y,t);

      const double phiRef = PHI(x,y,t);
      UNUSED(phiRef);
      const double Mphi   = CONSTANT_MOBILITY(M, phiRef);
      const double dMdphi = DCONSTANT_MOBILITY_DPHI(M,phiRef);

      const double gradPhi[2] = GRADPHI(x,y,t);
      const double gradMu[2]  = GRADMU(x,y,t);
      const double gradP[2]   = GRADP(x,y,t);
      const double gradPhi_dot_gradMuPlusAlphaP = gradPhi[0] * (gradMu[0] + alpha * gradP[0])
                                                + gradPhi[1] * (gradMu[1] + alpha * gradP[1]);

      return - (div_u - alpha * (Mphi * (lapMu + alpha * lapP) + dMdphi * gradPhi_dot_gradMuPlusAlphaP));
    }

    void uSource_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
    {
      DECLARE_VARIABLES(args, par);

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

    double phiSource_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      DECLARE_VARIABLES(args, par);
      const double phiRef      = PHI(x,y,t);
      const double u[2]        = UVEC(x,y,t);
      const double gradphi[2]  = GRADPHI(x,y,t);
      const double dphidt      = DPHIDT(x,y,t);
      const double lapP        = LAP_P(x,y,t);
      const double lapMu       = LAP_MU(x,y,t);
      const double uDotGradPhi = u[0]*gradphi[0] + u[1]*gradphi[1];
      const double divu        = DIVU(x,y,t);
      const double div_phiu    = phiRef * divu + uDotGradPhi;

      const double Mphi   = CONSTANT_MOBILITY(M, phiRef);
      const double dMdphi = DCONSTANT_MOBILITY_DPHI(M,phiRef);

      const double gradPhi[2] = GRADPHI(x,y,t);
      const double gradMu[2]  = GRADMU(x,y,t);
      const double gradP[2]   = GRADP(x,y,t);
      const double gradPhi_dot_gradMuPlusAlphaP = gradPhi[0] * (gradMu[0] + alpha * gradP[0])
                                                + gradPhi[1] * (gradMu[1] + alpha * gradP[1]);

      return - (dphidt + div_phiu - (Mphi * (lapMu + alpha * lapP) + dMdphi * gradPhi_dot_gradMuPlusAlphaP));
    }

    double muSource_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      DECLARE_VARIABLES(args, par);
      const double mu      = MU(x,y,t);
      const double phiRef  = PHI(x,y,t);
      const double phiCube = PHI(x,y,t) * PHI(x,y,t) * PHI(x,y,t);
      const double lapPhi  = LAP_PHI(x,y,t);

      return - (mu + tau * lapPhi - beta * (phiCube - phiRef));
    }

    int solve(const std::string &meshFile,
              const int timeStepFactor,
              const int orderU,
              const int orderPhiMu,
              const int degreeQuadrature,
              const feNLSolverOptions &NLoptions,
              int &numInteriorElements,
              int &numVertices,
              double &L2ErrorU,
              double &L2ErrorP,
              double &L2ErrorPhi,
              double &L2ErrorMu,
              double &integralP,
              double &H1ErrorU,
              int &timesteps)
    {
      const double rhoA    = 1.;
      const double rhoB    = 1000.;
      const double viscA   = 1.;
      const double viscB   = 1000.;
      const double M       = 1.;
      const double gamma   = 1.;
      const double epsilon = 1.;

      const double alpha = (rhoB - rhoA) / (rhoA + rhoB);
      const double beta   = 3. / (2. * sqrt(2.)) * gamma / epsilon;
      const double tau    = 3. / (2. * sqrt(2.)) * gamma * epsilon;

      std::vector<double> params = {rhoA, rhoB, viscA, viscB, M, gamma, epsilon, beta, tau, alpha};

      feVectorFunction uSol(  uSol_f, params);
      feFunction       pSol(  pSol_f, params);
      feFunction     phiSol(phiSol_f, params);
      feFunction      muSol( muSol_f, params);

      feVectorFunction grad_uSol(grad_uSol_f, params);

      feVectorFunction   uSource(  uSource_f, params);
      feFunction         pSource(  pSource_f, params);
      feFunction       phiSource(phiSource_f, params);
      feFunction        muSource( muSource_f, params); 

      feFunction   density(      rho_f, params);
      feFunction  drhodphi( drhodphi_f, params);
      feFunction viscosity(viscosity_f, params);
      feFunction dviscdphi(dviscdphi_f, params);

      feConstantVectorFunction volumeForce({0, -9.81});
      
      // // Continuity equation
      // feConstantFunction coeffDiffMu_p(alpha * M);
      // feConstantFunction  coeffDiffP_p(alpha * alpha * M);

      // // Phi equation
      // feConstantFunction coeffDiffMu_phi(M);
      // feConstantFunction  coeffDiffP_phi(M * alpha);

      feFunction mobility(mobility_f, params);
      feFunction mobilityAlpha(mobilityAlpha_f, params);
      feFunction mobilityAlphaSquared(mobilityAlphaSquared_f, params);
      
      // Mu equation
      feConstantFunction coeffDiffPhi(- tau);
      feConstantFunction coeffPhiCube(- beta);
      feConstantFunction coeffPhi    (+ beta);

      feMesh2DP1 mesh(meshFile);
      numInteriorElements = mesh.getNumInteriorElements();
      numVertices = mesh.getNumVertices();

      // Velocity and pressure
      feSpace *u = nullptr, *uB = nullptr, *p = nullptr,*pB = nullptr; // *pP = nullptr;
      // Phase marker and chemical potential
      feSpace *phi = nullptr, *phiB = nullptr, *mu = nullptr,  *muB = nullptr;
      // feSpace *l;

      // Initialize with exact initial condition
      // Use Taylor-Hood elements for (u,p) and equal order for (phi,u)
      int dQ = degreeQuadrature;
      feCheck(createFiniteElementSpace( u, &mesh, elementType::VECTOR_LAGRANGE, orderU, "U",       "Domaine", dQ, &uSol));
      feCheck(createFiniteElementSpace(uB, &mesh, elementType::VECTOR_LAGRANGE, orderU, "U",          "Bord", dQ, &uSol));
      feCheck(createFiniteElementSpace( p, &mesh, elementType::LAGRANGE,      orderU-1, "P",       "Domaine", dQ, &pSol));
      // feCheck(createFiniteElementSpace(pP, &mesh, elementType::LAGRANGE,      orderU-1, "P", "PointPression", dQ, &pSol));
      feCheck(createFiniteElementSpace(pB, &mesh, elementType::LAGRANGE,      orderU-1, "P",          "Bord", dQ, &pSol));

      feCheck(createFiniteElementSpace( phi, &mesh, elementType::LAGRANGE, orderPhiMu, "Phi", "Domaine", dQ, &phiSol));
      feCheck(createFiniteElementSpace(phiB, &mesh, elementType::LAGRANGE, orderPhiMu, "Phi",    "Bord", dQ, &phiSol));
      feCheck(createFiniteElementSpace(  mu, &mesh, elementType::LAGRANGE, orderPhiMu,  "Mu", "Domaine", dQ, &muSol));
      feCheck(createFiniteElementSpace( muB, &mesh, elementType::LAGRANGE, orderPhiMu,  "Mu",    "Bord", dQ, &muSol));

      // feCheck(createFiniteElementSpace(   l, &mesh, elementType::LAGRANGE, orderU-1,   "L", "Domaine", dQ, &scalarConstant::zero));

      std::vector<feSpace*> spaces = {u, uB, p, pB, phi, phiB, mu, muB};
      std::vector<feSpace*> essentialSpaces = {uB, pB, phiB, muB};

      feMetaNumber numbering(&mesh, spaces, essentialSpaces);
      feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

      std::vector<feBilinearForm*> forms;

      // Continuity
      feBilinearForm *divU = nullptr, *diffMu_p = nullptr, *diffP_p = nullptr, *pS = nullptr;
      feCheck(createBilinearForm(    divU,  {p, u}, new feSysElm_MixedDivergence<2>(&scalarConstant::one)));
      // feCheck(createBilinearForm(diffMu_p, {p, mu}, new feSysElm_MixedGradGrad<2>(&mobilityAlpha)));
      // feCheck(createBilinearForm( diffP_p,     {p}, new feSysElm_Diffusion<2>(&mobilityAlphaSquared)));
      feCheck(createBilinearForm(diffMu_p, {p, mu, phi}, new feSysElm_MixedGradGradFieldDependentCoeff<2>(&mobilityAlpha)));
      feCheck(createBilinearForm( diffP_p,     {p, phi}, new feSysElm_DiffusionFieldDependentCoeff<2>(&mobilityAlphaSquared)));
      feCheck(createBilinearForm(      pS,     {p}, new feSysElm_Source(&pSource)));
      forms.push_back(divU);
      forms.push_back(diffMu_p);
      forms.push_back(diffP_p);
      forms.push_back(pS);

      // Momentum
      feBilinearForm *CHNS_momentum = nullptr, *uS = nullptr;
      feCheck(createBilinearForm(CHNS_momentum, {u, p, phi, mu}, new feSysElm_CHNS_Momentum_Alternative<2>(&density,
                                                                                                           &drhodphi,
                                                                                                           &viscosity,
                                                                                                           &dviscdphi,
                                                                                                           &volumeForce)));
      feCheck(createBilinearForm(           uS,             {u}, new feSysElm_VectorSource<2>(&uSource)));
      forms.push_back(CHNS_momentum);
      forms.push_back(uS);

      // Tracer convection
      feBilinearForm *dphidt = nullptr, *divPhiU = nullptr, *diffMu_phi = nullptr, *diffP_phi = nullptr, *phiS = nullptr;
      feCheck(createBilinearForm(    dphidt,     {phi}, new feSysElm_TransientMass(&scalarConstant::one)));
      feCheck(createBilinearForm(   divPhiU,  {phi, u}, new feSysElm_ScalarVectorProduct<2>(&scalarConstant::minusOne)));
      // feCheck(createBilinearForm(diffMu_phi, {phi, mu}, new feSysElm_MixedGradGrad<2>(&mobility)));
      // feCheck(createBilinearForm( diffP_phi,  {phi, p}, new feSysElm_MixedGradGrad<2>(&mobilityAlpha)));
      feCheck(createBilinearForm(diffMu_phi, {phi, mu, phi}, new feSysElm_MixedGradGradFieldDependentCoeff<2>(&mobility)));
      feCheck(createBilinearForm( diffP_phi,  {phi, p, phi}, new feSysElm_MixedGradGradFieldDependentCoeff<2>(&mobilityAlpha)));
      feCheck(createBilinearForm(      phiS,     {phi}, new feSysElm_Source(&phiSource)));
      forms.push_back(dphidt);
      forms.push_back(divPhiU);
      forms.push_back(diffMu_phi);
      forms.push_back(diffP_phi);
      forms.push_back(phiS);

      // Chemical potential
      feBilinearForm *muMass = nullptr, *phiCube = nullptr, *phiMass = nullptr, *diffPhi = nullptr, *muS = nullptr;
      feCheck(createBilinearForm( muMass,      {mu}, new feSysElm_Mass(&scalarConstant::one)));
      feCheck(createBilinearForm(diffPhi, {mu, phi}, new feSysElm_MixedGradGrad<2>(&coeffDiffPhi)));
      feCheck(createBilinearForm(phiCube, {mu, phi}, new feSysElm_MixedMassPower(&coeffPhiCube, 3)));
      feCheck(createBilinearForm(phiMass, {mu, phi}, new feSysElm_MixedMass(&coeffPhi)));
      feCheck(createBilinearForm(    muS,      {mu}, new feSysElm_Source(&muSource)));
      forms.push_back(muMass);
      forms.push_back(diffPhi);
      forms.push_back(phiCube);
      forms.push_back(phiMass);
      forms.push_back(muS);

      // phiCube->setComputeMatrixWithFD(true);
      // phiMass->setComputeMatrixWithFD(true);

      //////////////////////////////////////////////////////////////
      // sol.initialize(&mesh);
      // double overallError = 0.;
      // std::string highestErrorForm;
      // for(feBilinearForm *form : forms) {
      //   if(form->getWeakFormName() != "CHNS_MOMENTUM_ALTERNATIVE") {
      //     feInfo("");
      //     feInfo("");
      //     feInfo("");
      //     feInfo("Forme %s", form->getWeakFormName().data());
      //     feInfo("");
      //     feInfo("");
      //     feInfo("");
      //     double error = 0.;
      //     for(int iElm = 0; iElm < mesh.getNumInteriorElements(); ++iElm) {
      //       feInfo("elm %d", iElm);
      //       double err = form->compareAnalyticalAndFDMatrices(&sol, iElm);
      //       error = fmax(error, err);
      //       if( err > overallError) {
      //         overallError = fmax(overallError, err);
      //         highestErrorForm = form->getWeakFormName();
      //       }
      //     }
      //     feInfo("Max error = %1.6e", error);
      //   }
      // }
      // feInfo("Max overall error = %1.6e for form %s", overallError, highestErrorForm.data());
      // exit(-1);
      //////////////////////////////////////////////////////////////
      // divPhiU->setComputeMatrixWithFD(false);
      // // Zero-average pressure field
      // feBilinearForm *massLP, *massPL;
      // feCheck(createBilinearForm(massLP, {l, p}, new feSysElm_MixedMass(&scalarConstant::one)));
      // feCheck(createBilinearForm(massPL, {p, l}, new feSysElm_MixedMass(&scalarConstant::one)));
      // forms.push_back(massLP);
      // forms.push_back(massPL);

      feLinearSystem *system;
      #if defined(HAVE_MKL)
        feCheck(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
        feCheck(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
      #else
        feCheck(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
      #endif

      feNorm *errorU_L2 = nullptr, *errorP_L2 = nullptr, *errorPhi_L2 = nullptr, *errorMu_L2 = nullptr;
      feNorm *errorU_H1 = nullptr;

      feCheck(createNorm(  errorU_L2, VECTOR_L2_ERROR,   {u}, &sol, nullptr, &uSol));
      feCheck(createNorm(  errorP_L2,        L2_ERROR,   {p}, &sol, &pSol));
      feCheck(createNorm(errorPhi_L2,        L2_ERROR, {phi}, &sol, &phiSol));
      feCheck(createNorm( errorMu_L2,        L2_ERROR,  {mu}, &sol, &muSol));

      feCheck(createNorm(  errorU_H1, VECTOR_SEMI_H1_ERROR, {u}, &sol, nullptr, &grad_uSol));

      // feCheck(createNorm(  errorU_L2, VECTOR_LINF_ERROR,   {u}, &sol, nullptr, &uSol));
      // feCheck(createNorm(  errorP_L2,        LINF_ERROR,   {p}, &sol, &pSol));
      // feCheck(createNorm(errorPhi_L2,        LINF_ERROR, {phi}, &sol, &phiSol));
      // feCheck(createNorm( errorMu_L2,        LINF_ERROR,  {mu}, &sol, &muSol));

      std::vector<feNorm *> norms = {errorU_L2, errorP_L2, errorPhi_L2, errorMu_L2};

      // Reference test file matches those parameters and L2 error evaluated at t1:
      UNUSED(timeStepFactor);
      double t0 = 0.;
      double t1 = 1.;
      int nTimeSteps = 10;
      // int nTimeSteps = 10 * pow(2, timeStepFactor);

      // // For other convergence studies:
      // const double t0 = 0.;
      // const double t1 = 1e-3;
      // int nTimeSteps = 10 * pow(2, timeStepFactor);

      timesteps = nTimeSteps;
      TimeIntegrator *solver;
      timeIntegratorScheme scheme = timeIntegratorScheme::BDF1;
      // timeIntegratorScheme scheme = timeIntegratorScheme::BDF2;

      feCheck(createTimeIntegrator(solver, scheme, NLoptions,
        system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));

      // Start BDF2 with manufactured solution
      if(scheme == timeIntegratorScheme::BDF2) {
        static_cast<BDF2Integrator*>(solver)->setStartingMethod(BDF2Starter::InitialCondition);
      }

      // Solve
      feCheck(solver->makeSteps(nTimeSteps));

      // Use L2 at t1 to compare with reference test files:
      L2ErrorU   = errorU_L2->compute();
      L2ErrorP   = errorP_L2->compute();
      L2ErrorPhi = errorPhi_L2->compute();
      L2ErrorMu  = errorMu_L2->compute();

      H1ErrorU   = errorU_H1->compute();

      // Integral of pressure (target is 0)
      feNorm *intP;
      feCheckReturn(createNorm(intP, INTEGRAL, {p}, &sol));
      integralP = intP->compute();
      delete intP;

      // Or use L1 or Linf error in time:
      // std::vector<std::vector<double>> postProc = solver->getPostProcessingData();
      // L2ErrorU = 0.;
      // L2ErrorP = 0.;
      // L2ErrorPhi = 0.;
      // L2ErrorMu = 0.;
      // for(int j = 1; j < nTimeSteps+1; ++j) {

      //   // L1 norm in time, method of rectangles
      //   L2ErrorU   += dt * postProc[1][j];
      //   L2ErrorP   += dt * postProc[2][j];
      //   L2ErrorPhi += dt * postProc[3][j];
      //   L2ErrorMu  += dt * postProc[4][j];

      //   // // Linf norm in time
      //   // L2ErrorU = fmax(L2ErrorU, postProc[1][j]);
      //   // L2ErrorP = fmax(L2ErrorP, postProc[2][j]);
      //   // L2ErrorPhi = fmax(L2ErrorPhi, postProc[3][j]);
      //   // L2ErrorMu = fmax(L2ErrorMu, postProc[4][j]);
      // }

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
      feNLSolverOptions NLoptions{1e-10, 1e-10, 1e8, 50, 4, 1e-1};

      std::vector<int> nElm(numMeshes);
      std::vector<int> timesteps(numMeshes);
      std::vector<int> nVertices(numMeshes);
      std::vector<double> errU(numMeshes, 0.);
      std::vector<double> errP(numMeshes, 0.);
      std::vector<double> errPhi(numMeshes, 0.);
      std::vector<double> errMu(numMeshes, 0.);
      std::vector<double> integralP(numMeshes, 0.);
      std::vector<double> errU_H1(numMeshes, 0.);
      // std::vector<double> errP_H1(numMeshes, 0.);
      // std::vector<double> errPhi_H1(numMeshes, 0.);
      // std::vector<double> errMu_H1(numMeshes, 0.);

      bool writeErrorToFile = true;
      std::string outputDirectory = "./";
      std::ofstream errorFileL2, errorFileH1, timeErrorFile;
      if(writeErrorToFile) {
        errorFileL2.open(outputDirectory + "errorL2.txt");
        errorFileH1.open(outputDirectory + "errorH1.txt");
        timeErrorFile.open(outputDirectory + "timeErrorFile.txt");
      }

      for(int i = 0; i < numMeshes; ++i)
      {
        const int timeStepFactor = i;
        // std::string meshFile = "../../../data/stokes" + std::to_string(i+1) + ".msh";
        // std::string meshFile = "../data/stokes" + std::to_string(i+1) + ".msh";
        std::string meshFile = "../data/transfiniteSquare" + std::to_string(i+1) + ".msh";
        alternativeCHNS_constantMobility::MMS::solve(meshFile, timeStepFactor, orderU, orderPhiMu,
          degreeQuadrature, NLoptions, nElm[i], nVertices[i],
          errU[i], errP[i], errPhi[i], errMu[i], integralP[i], errU_H1[i], timesteps[i]);
          // errU_H1[i], errP_H1[i], errPhi_H1[i], errMu_H1[i]);

        if(writeErrorToFile) {
          errorFileL2 << nVertices[i] << "\t" << nElm[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errU[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errP[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errPhi[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errMu[i] << "\n";
          errorFileL2.flush();
          timeErrorFile << timesteps[i] << "\t" << nElm[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errU[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errP[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errPhi[i] << "\t"
            << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errMu[i] << "\n";
          timeErrorFile.flush();
        }
      }

      if(writeErrorToFile) {
        errorFileL2.close();
        errorFileH1.close();
        timeErrorFile.close();
      }

      resultBuffer << "   (u,p) : Taylor-Hood elements P" << orderU << "-P" << orderU-1 << std::endl; 
      resultBuffer << "(phi,mu) : Equal order elements P" << orderPhiMu << "-P" << orderPhiMu << std::endl; 
      resultBuffer << "Cahn-Hilliard Navier-Stokes MMS - L2 Error on velocity" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errU, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      resultBuffer << "Cahn-Hilliard Navier-Stokes MMS - H1 Error on velocity" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errU_H1, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      resultBuffer << "Cahn-Hilliard Navier-Stokes MMS - L2 Error on pressure" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errP, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      resultBuffer << "Cahn-Hilliard Navier-Stokes MMS - L2 Error on phase marker phi" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errPhi, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      resultBuffer << "Cahn-Hilliard Navier-Stokes MMS - L2 Error on chemical potential mu" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errMu, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

      // resultBuffer << "Time convergence:" << std::endl; 
      // computeAndPrintConvergence(1, numMeshes, errU, timesteps, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      // computeAndPrintConvergence(1, numMeshes, errP, timesteps, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      // computeAndPrintConvergence(1, numMeshes, errPhi, timesteps, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      // computeAndPrintConvergence(1, numMeshes, errMu, timesteps, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

      return 0;
    }

  } // namespace MMS
} // namespace alternativeCHNS_constantMobility

#define MOBILITY(M, phi) M * (1. - phi * phi)
#define DMOBILITY_DPHI(M, phi) -2. * M * phi

struct alternativeCHNS_degenerateMobility {

  static void uSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);
    res[0] = U(x,y,t);
    res[1] = V(x,y,t);
  }

  static void grad_uSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);
    const double gradu[2][2] = GRADU(x,y,t);
    res[0] = gradu[0][0];
    res[1] = gradu[0][1];
    res[2] = gradu[1][0];
    res[3] = gradu[1][1];
  }

  static double pSol_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return P(x,y,t);
  }

  static void grad_pSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);
    const double gradp[2] = GRADP(x,y,t);
    res[0] = gradp[0];
    res[1] = gradp[1];
  }

  static double phiSol_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return PHI(x,y,t);
  }

  static void grad_phiSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);
    const double gradphi[2] = GRADPHI(x,y,t);
    res[0] = gradphi[0];
    res[1] = gradphi[1];
  }

  static double muSol_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return MU(x,y,t);
  }

  static void grad_muSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);
    const double gradmu[2] = GRADMU(x,y,t);
    res[0] = gradmu[0];
    res[1] = gradmu[1];
  }

  //
  // CHNS models for density, viscosity and mobility
  //
  static double rho_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return (rhoA - rhoB)/2. * phi + (rhoA + rhoB)/2.;
  }

  static double drhodphi_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return (rhoA - rhoB)/2.;
  }

  static double viscosity_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return (viscA - viscB)/2. * phi + (viscA + viscB)/2.;
  }

  static double dviscdphi_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return (viscA - viscB)/2.;
  }

  //
  // Degenerate mobility (0 in pure phases)
  // Here we take (1 - phi^2) as the dependence in phi,
  // that is:
  //
  // M(phi) = M * (1 - phi^2).
  //
  // Since phi is a product of sines, phi is in [-1,1] and M(phi) >= 0.
  //
  static double mobility_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return MOBILITY(M, phi);
  }

  static double mobilityAlpha_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return alpha * MOBILITY(M, phi);
  }

  static double mobilityAlphaSquared_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
    return alpha * alpha * MOBILITY(M, phi);
  }

  static double pSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);

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

    return - (div_u - alpha * (Mphi * (lapMu + alpha * lapP) + dMdphi * gradPhi_dot_gradMuPlusAlphaP));
  }

  static void uSource_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    DECLARE_VARIABLES(args, par);

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

  static double phiSource_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    DECLARE_VARIABLES(args, par);
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
    const double mu      = MU(x,y,t);
    const double phiRef  = PHI(x,y,t);
    const double phiCube = PHI(x,y,t) * PHI(x,y,t) * PHI(x,y,t);
    const double lapPhi  = LAP_PHI(x,y,t);

    return - (mu + tau * lapPhi - beta * (phiCube - phiRef));
  }

  // alternativeCHNS_degenerateMobility
  static feStatus solve(const std::string &meshFile,
                        const int timeStepFactor,
                        const int orderU,
                        const int orderPhiMu,
                        const int degreeQuadrature,
                        const feNLSolverOptions &NLoptions,
                        int &numInteriorElements,
                        int &numVertices,
                        double &L2ErrorU,
                        double &L2ErrorP,
                        double &L2ErrorPhi,
                        double &L2ErrorMu,
                        double &integralP,
                        double &H1ErrorU,
                        double &H1ErrorP,
                        double &H1ErrorPhi,
                        double &H1ErrorMu,
                        int &timesteps)
  {
    const double rhoA    = 1.;
    const double rhoB    = 1.;
    const double viscA   = 1.;
    const double viscB   = 1.;
    const double M       = 1.;
    const double gamma   = 1.;
    const double epsilon = 1.;

    const double alpha = (rhoB - rhoA) / (rhoA + rhoB);
    const double beta   = 3. / (2. * sqrt(2.)) * gamma / epsilon;
    const double tau    = 3. / (2. * sqrt(2.)) * gamma * epsilon;

    std::vector<double> CHNSparameters = {alpha, gamma, epsilon};

    std::vector<double> params = {rhoA, rhoB, viscA, viscB, M, gamma, epsilon, beta, tau, alpha};

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

    feFunction   density(      rho_f, params);
    feFunction  drhodphi( drhodphi_f, params);
    feFunction viscosity(viscosity_f, params);
    feFunction dviscdphi(dviscdphi_f, params);

    feFunction mobility(mobility_f, params);
    feFunction mobilityAlpha(mobilityAlpha_f, params);
    feFunction mobilityAlphaSquared(mobilityAlphaSquared_f, params);

    feConstantVectorFunction volumeForce({0, -9.81});
    
    // Mu equation
    feConstantFunction coeffDiffPhi(- tau);
    feConstantFunction coeffPhiCube(- beta);
    feConstantFunction coeffPhi    (+ beta);

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();
    numVertices = mesh.getNumVertices();

    // Velocity and pressure
    feSpace *u = nullptr, *uB = nullptr, *p = nullptr,*pB = nullptr;
    // Phase marker and chemical potential
    feSpace *phi = nullptr, *phiB = nullptr, *mu = nullptr,  *muB = nullptr;

    // Initialize with exact initial condition
    // Use Taylor-Hood elements for (u,p) and equal order for (phi,u)
    int dQ = degreeQuadrature;
    feCheckReturn(createFiniteElementSpace( u, &mesh, elementType::VECTOR_LAGRANGE, orderU, "U",  "Domaine", dQ, &uSol));
    feCheckReturn(createFiniteElementSpace(uB, &mesh, elementType::VECTOR_LAGRANGE, orderU, "U",     "Bord", dQ, &uSol));
    feCheckReturn(createFiniteElementSpace( p, &mesh, elementType::LAGRANGE,      orderU-1, "P",  "Domaine", dQ, &pSol));
    feCheckReturn(createFiniteElementSpace(pB, &mesh, elementType::LAGRANGE,      orderU-1, "P",     "Bord", dQ, &pSol));

    feCheckReturn(createFiniteElementSpace( phi, &mesh, elementType::LAGRANGE, orderPhiMu, "Phi", "Domaine", dQ, &phiSol));
    feCheckReturn(createFiniteElementSpace(phiB, &mesh, elementType::LAGRANGE, orderPhiMu, "Phi",    "Bord", dQ, &phiSol));
    feCheckReturn(createFiniteElementSpace(  mu, &mesh, elementType::LAGRANGE, orderPhiMu,  "Mu", "Domaine", dQ, &muSol));
    feCheckReturn(createFiniteElementSpace( muB, &mesh, elementType::LAGRANGE, orderPhiMu,  "Mu",    "Bord", dQ, &muSol));

    std::vector<feSpace*> spaces = {u, uB, p, pB, phi, phiB, mu, muB};
    std::vector<feSpace*> essentialSpaces = {uB, pB, phiB, muB};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

    //////////////////////////////////////////////////////////////////////////
    // std::vector<feBilinearForm*> forms;

    // // Continuity
    // feBilinearForm *divU = nullptr, *diffMu_p = nullptr, *diffP_p = nullptr, *pS = nullptr;
    // feCheckReturn(createBilinearForm(    divU,  {p, u}, new feSysElm_MixedDivergence<2>(&scalarConstant::one)));
    // // feCheckReturn(createBilinearForm(diffMu_p, {p, mu}, new feSysElm_MixedGradGrad<2>(&coeffDiffMu_p)));
    // // feCheckReturn(createBilinearForm( diffP_p,     {p}, new feSysElm_Diffusion<2>(&coeffDiffP_p)));
    // feCheckReturn(createBilinearForm(diffMu_p, {p, mu, phi}, new feSysElm_MixedGradGradFieldDependentCoeff<2>(&mobilityAlpha)));
    // feCheckReturn(createBilinearForm( diffP_p,     {p, phi}, new feSysElm_DiffusionFieldDependentCoeff<2>(&mobilityAlphaSquared)));
    // feCheckReturn(createBilinearForm(      pS,     {p}, new feSysElm_Source(&pSource)));
    // forms.push_back(divU);
    // forms.push_back(diffMu_p);
    // forms.push_back(diffP_p);
    // forms.push_back(pS);

    // // Momentum
    // feBilinearForm *CHNS_momentum = nullptr, *uS = nullptr;
    // feCheckReturn(createBilinearForm(CHNS_momentum, {u, p, phi, mu}, new feSysElm_CHNS_Momentum_Alternative<2>(&density,
    //                                                                                                      &drhodphi,
    //                                                                                                      &viscosity,
    //                                                                                                      &dviscdphi,
    //                                                                                                      &volumeForce)));
    // feCheckReturn(createBilinearForm(           uS,             {u}, new feSysElm_VectorSource<2>(&uSource)));
    // forms.push_back(CHNS_momentum);
    // forms.push_back(uS);

    // // Tracer convection
    // feBilinearForm *dphidt = nullptr, *divPhiU = nullptr, *diffMu_phi = nullptr, *diffP_phi = nullptr, *phiS = nullptr;
    // feCheckReturn(createBilinearForm(    dphidt,     {phi}, new feSysElm_TransientMass(&scalarConstant::one)));
    // feCheckReturn(createBilinearForm(   divPhiU,  {phi, u}, new feSysElm_ScalarVectorProduct<2>(&scalarConstant::minusOne)));
    // feCheckReturn(createBilinearForm(diffMu_phi, {phi, mu, phi}, new feSysElm_MixedGradGradFieldDependentCoeff<2>(&mobility)));
    // feCheckReturn(createBilinearForm( diffP_phi,  {phi, p, phi}, new feSysElm_MixedGradGradFieldDependentCoeff<2>(&mobilityAlpha)));
    // feCheckReturn(createBilinearForm(      phiS,     {phi}, new feSysElm_Source(&phiSource)));
    // forms.push_back(dphidt);
    // forms.push_back(divPhiU);
    // forms.push_back(diffMu_phi);
    // forms.push_back(diffP_phi);
    // forms.push_back(phiS);

    // // Chemical potential
    // feBilinearForm *muMass = nullptr, *phiCube = nullptr, *phiMass = nullptr, *diffPhi = nullptr, *muS = nullptr;
    // feCheckReturn(createBilinearForm( muMass,      {mu}, new feSysElm_Mass(&scalarConstant::one)));
    // feCheckReturn(createBilinearForm(diffPhi, {mu, phi}, new feSysElm_MixedGradGrad<2>(&coeffDiffPhi)));
    // feCheckReturn(createBilinearForm(phiCube, {mu, phi}, new feSysElm_MixedMassPower(&coeffPhiCube, 3)));
    // feCheckReturn(createBilinearForm(phiMass, {mu, phi}, new feSysElm_MixedMass(&coeffPhi)));
    // feCheckReturn(createBilinearForm(    muS,      {mu}, new feSysElm_Source(&muSource)));
    // forms.push_back(muMass);
    // forms.push_back(diffPhi);
    // forms.push_back(phiCube);
    // forms.push_back(phiMass);
    // forms.push_back(muS);
    //////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    // Or with a single weak form:
    feBilinearForm *CHNS;
    feCheckReturn(createBilinearForm(CHNS, {u, p, phi, mu}, new feSysElm_CHNS_Alternative<2>(
      &density,
      &drhodphi,
      &viscosity,
      &dviscdphi,
      &mobility,
      &volumeForce,
      &pSource,
      &uSource,
      &phiSource,
      &muSource,
      CHNSparameters)));
    std::vector<feBilinearForm*> forms = {CHNS};
    ///////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////
    // sol.initialize(&mesh);
    // double overallError = 0.;
    // std::string highestErrorForm;
    // for(feBilinearForm *form : forms) {
    //   if(form->getWeakFormName() != "CHNS_MOMENTUM_ALTERNATIVE") {
    //     feInfo("");
    //     feInfo("");
    //     feInfo("");
    //     feInfo("Forme %s", form->getWeakFormName().data());
    //     feInfo("");
    //     feInfo("");
    //     feInfo("");
    //     double error = 0.;
    //     for(int iElm = 0; iElm < mesh.getNumInteriorElements(); ++iElm) {
    //       feInfo("elm %d", iElm);
    //       double err = form->compareAnalyticalAndFDMatrices(&sol, iElm);
    //       error = fmax(error, err);
    //       if( err > overallError) {
    //         overallError = fmax(overallError, err);
    //         highestErrorForm = form->getWeakFormName();
    //       }
    //     }
    //     feInfo("Max error = %1.6e", error);
    //   }
    // }
    // feInfo("Max overall error = %1.6e for form %s", overallError, highestErrorForm.data());
    // exit(-1);
    //////////////////////////////////////////////////////////////
    // divPhiU->setComputeMatrixWithFD(false);

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
    #else
      feCheckReturn(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
    #endif

    feNorm *errorU_L2 = nullptr, *errorP_L2 = nullptr, *errorPhi_L2 = nullptr, *errorMu_L2 = nullptr;
    feNorm *errorU_H1 = nullptr, *errorP_H1 = nullptr, *errorPhi_H1 = nullptr, *errorMu_H1 = nullptr;

    feCheckReturn(createNorm(  errorU_L2, VECTOR_L2_ERROR,   {u}, &sol, nullptr, &uSol));
    feCheckReturn(createNorm(  errorP_L2,        L2_ERROR,   {p}, &sol, &pSol));
    feCheckReturn(createNorm(errorPhi_L2,        L2_ERROR, {phi}, &sol, &phiSol));
    feCheckReturn(createNorm( errorMu_L2,        L2_ERROR,  {mu}, &sol, &muSol));

    feCheckReturn(createNorm(  errorU_H1, VECTOR_SEMI_H1_ERROR,   {u}, &sol, nullptr, &grad_uSol));
    feCheckReturn(createNorm(  errorP_H1,        SEMI_H1_ERROR,   {p}, &sol, nullptr, &grad_pSol));
    feCheckReturn(createNorm(errorPhi_H1,        SEMI_H1_ERROR, {phi}, &sol, nullptr, &grad_phiSol));
    feCheckReturn(createNorm( errorMu_H1,        SEMI_H1_ERROR,  {mu}, &sol, nullptr, &grad_muSol));

    std::vector<feNorm *> norms = {};

    // Reference test file matches those parameters and L2 error evaluated at t1:
    UNUSED(timeStepFactor);
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 20;
    // int nTimeSteps = 10 * pow(2, timeStepFactor);

    timesteps = nTimeSteps;
    TimeIntegrator *solver;
    // timeIntegratorScheme scheme = timeIntegratorScheme::BDF1;
    timeIntegratorScheme scheme = timeIntegratorScheme::BDF2;

    feCheckReturn(createTimeIntegrator(solver, scheme, NLoptions,
      system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));

    // Start BDF2 with manufactured solution
    if(scheme == timeIntegratorScheme::BDF2) {
      static_cast<BDF2Integrator*>(solver)->setStartingMethod(BDF2Starter::InitialCondition);
    }

    // Solve
    feCheckReturn(solver->makeSteps(nTimeSteps));

    // Use L2 at t1 to compare with reference test files:
    L2ErrorU   = errorU_L2->compute();
    L2ErrorP   = errorP_L2->compute();
    L2ErrorPhi = errorPhi_L2->compute();
    L2ErrorMu  = errorMu_L2->compute();

    H1ErrorU   = errorU_H1->compute();
    H1ErrorP   = errorP_H1->compute();
    H1ErrorPhi = errorPhi_H1->compute();
    H1ErrorMu  = errorMu_H1->compute();

    // Integral of pressure (target is 0)
    feNorm *intP;
    feCheckReturn(createNorm(intP, INTEGRAL, {p}, &sol));
    integralP = fabs(intP->compute());
    delete intP;

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

    return FE_STATUS_OK;
  }
}; // alternativeCHNS_degenerateMobility

template<typename problemType>
feStatus meshConvergence(std::stringstream &resultBuffer,
                         const std::string &meshFileRoot,
                         const int orderU,
                         const int orderPhiMu,
                         const int numMeshes,
                         const int degreeQuadrature)
{
  feNLSolverOptions NLoptions{1e-10, 1e-10, 1e8, 50, 4, 1e-1};

  std::vector<int> nElm(numMeshes);
  std::vector<int> timesteps(numMeshes);
  std::vector<int> nVertices(numMeshes);
  std::vector<double> errU(numMeshes, 0.);
  std::vector<double> errP(numMeshes, 0.);
  std::vector<double> errPhi(numMeshes, 0.);
  std::vector<double> errMu(numMeshes, 0.);

  std::vector<double> integralP(numMeshes, 0.);
  
  std::vector<double> errU_H1(numMeshes, 0.);
  std::vector<double> errP_H1(numMeshes, 0.);
  std::vector<double> errPhi_H1(numMeshes, 0.);
  std::vector<double> errMu_H1(numMeshes, 0.);

  bool writeErrorToFile = true;
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

  for(int i = 0; i < numMeshes; ++i)
  {
    const int timeStepFactor = i;

    #if defined(TEST_TIME_CONVERGENCE)
      std::string meshFile = meshFileRoot + "1.msh";
    #else
      std::string meshFile = meshFileRoot + std::to_string(i+1) + ".msh";
    #endif

    feCheckReturn(problemType::ezsolve(meshFile, timeStepFactor, orderU, orderPhiMu,
                       degreeQuadrature, NLoptions, nElm[i], nVertices[i],
                       errU[i],
                       errP[i],
                       errPhi[i],
                       errMu[i],
                       integralP[i],
                       errU_H1[i],
                       errP_H1[i],
                       errPhi_H1[i],
                       errMu_H1[i],
                       timesteps[i]));

    if(writeErrorToFile) {
    #if defined(TEST_TIME_CONVERGENCE)
      timeErrorFile << timesteps[i] << "\t" << nElm[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errU_H1[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errP[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errPhi_H1[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errMu_H1[i] << "\n";
      timeErrorFile.flush();
    #else
      errorFileL2 << nVertices[i] << "\t" << nElm[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errU[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errP[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errPhi[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errMu[i] << "\n";
      errorFileL2.flush();
      errorFileH1 << nVertices[i] << "\t" << nElm[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errU_H1[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errP_H1[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errPhi_H1[i] << "\t"
        << std::setprecision(DEFAULT_SIGNIFICANT_DIGITS) << errMu_H1[i] << "\n";
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
  resultBuffer << "   (u,p) : Taylor-Hood elements P" << orderU << "-P" << orderU-1 << std::endl; 
  resultBuffer << "(phi,mu) : Equal order elements P" << orderPhiMu << "-P" << orderPhiMu << std::endl;

  resultBuffer << "==================== L2 error norms ==================" << std::endl;
  resultBuffer << "U Velocity" << std::endl;
  computeAndPrintConvergence(2, numMeshes, errU, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  resultBuffer << "P Pressure" << std::endl;
  computeAndPrintConvergence(2, numMeshes, errP, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  resultBuffer << "Phi Phase marker" << std::endl;
  computeAndPrintConvergence(2, numMeshes, errPhi, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  resultBuffer << "Mu Chemical potential" << std::endl;
  computeAndPrintConvergence(2, numMeshes, errMu, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

  resultBuffer << "==================== H1 error norms ==================" << std::endl;
  resultBuffer << "U Velocity" << std::endl;
  computeAndPrintConvergence(2, numMeshes, errU_H1, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  resultBuffer << "P Pressure" << std::endl;
  computeAndPrintConvergence(2, numMeshes, errP_H1, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  resultBuffer << "Phi Phase marker" << std::endl;
  computeAndPrintConvergence(2, numMeshes, errPhi_H1, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  resultBuffer << "Mu Chemical potential" << std::endl;
  computeAndPrintConvergence(2, numMeshes, errMu_H1, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

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
  cahnHilliard::MMS::meshConvergence(resultBuffer, order, 4, degreeQuadrature);
  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

// TEST(CahnHilliard, CHNS)
// {
//   initialize(my_argc, my_argv);
//   setVerbose(1);
  
//   std::string testRoot = "../../../tests/withLinearSolver/CHNS_MMS";

//   std::stringstream resultBuffer;

//   int degreeQuadrature = 8;
//   int orderU = 2;
//   int orderPhiMu = 1;
//   cahnHilliardNavierStokes::MMS::meshConvergence(resultBuffer, orderU, orderPhiMu, 4, degreeQuadrature);
//   EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
//   finalize();
// }

TEST(CahnHilliard, CHNS_VolumeAveraged)
{
  initialize(my_argc, my_argv);
  setVerbose(1);
  
  std::string testRoot = "../../../tests/withLinearSolver/CHNS_VolumeAveraged";
  std::string meshFileRoot = "../../../data/mmsMeshes/square";

  std::stringstream resultBuffer;

  int degreeQuadrature = 8;
  int orderU = 2;
  int orderPhiMu = 1;
  ASSERT_TRUE(meshConvergence<CahnHilliardVolumeAveraged>(resultBuffer, 
    meshFileRoot, orderU, orderPhiMu, 3, degreeQuadrature) == FE_STATUS_OK);
  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

// TEST(CahnHilliard, CHNS_Alt_Constant_Mobility)
// {
//   initialize(my_argc, my_argv);
//   setVerbose(1);
  
//   std::string testRoot = "../../../tests/withLinearSolver/CHNS_Alt_MMS";

//   std::stringstream resultBuffer;

//   int degreeQuadrature = 8;
//   int orderU = 2;
//   int orderPhiMu = 1;
//   alternativeCHNS_constantMobility::MMS::meshConvergence(resultBuffer, orderU, orderPhiMu, 5, degreeQuadrature);
//   EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
//   finalize();
// }

// TEST(CahnHilliard, CHNS_Alt_Degenerate_Mobility)
// {
//   initialize(my_argc, my_argv);
//   setVerbose(1);
  
//   std::string testRoot = "../../../tests/withLinearSolver/CHNS_Alt_MMS";

//   std::stringstream resultBuffer;

//   int degreeQuadrature = 12;
//   int orderU = 2;
//   int orderPhiMu = 1;
//   // alternativeCHNS_degenerateMobility::MMS::meshConvergence(resultBuffer, orderU, orderPhiMu, 4, degreeQuadrature);
//   ASSERT_TRUE(meshConvergence<alternativeCHNS_degenerateMobility>(resultBuffer, orderU, orderPhiMu, 4, degreeQuadrature) == FE_STATUS_OK);
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