
#include "feAPI.h"

/*
  Poiseuille flow in a rectangular channel with dimensions L x H
  under constant pressure gradient dpdx = -1.
  Solves either the steady Stokes equations or the steady Navier-Stokes equations,
  depending on if the convective term is added to the vector of weak forms.
*/

double solP(const feFunctionArguments &args, const std::vector<double> &par)
{
  double x = args.pos[0];
  double L = par[0];
  double dpdx = par[1];
  return -dpdx * (L-x);
}

void solU(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
{
  double y = args.pos[1];
  double H = par[0];
  double dpdx = par[1];
  res[0] = -dpdx/2. * y * (H - y);
  res[1] = 0.;
}

int main(int argc, char **argv)
{
  initialize(argc, argv);

  // Default parameters
  const char* meshFile = "../data/poiseuille0.msh";
  int orderVelocity = 2;
  int orderPressure = 1;
  int degreeQuadrature = 10;
  const char *solverType = "PETSc_MUMPS";
  bool divergenceFormulation = true;

  feOptionsParser options(argc, argv);
  options.addOption(&meshFile, "-m", "--mesh", "Mesh file (gmsh format .msh only)");
  options.addOption(&orderVelocity, "-o", "--order_velocity", "Polynomial order for the velocity");
  options.addOption(&solverType, "-linsol", "--linear_solver_type", "Choice of linear solver (PETSc or MKL Pardiso)");
  options.addOption(&divergenceFormulation, "-div", "--enable_divergence_formulation",
                                          "-nodiv", "--disable_divergence_formulation","Divergence or Laplacian formulation");
  feCheck(options.parse());

  orderPressure = fmin(1, orderVelocity-1);

  setVerbose(1);

  feConstantFunction zero(0.);
  feConstantFunction one(1.);
  feConstantFunction minusOne(-1.);
  feConstantVectorFunction zeroVector({0., 0.});
  feConstantVectorFunction  oneVector({1., 1.});
  feConstantVectorFunction  uniform({1., 0.});

  double Re = 1.;
  feConstantFunction viscosity(1./Re);
  feConstantFunction density(1.);

  double H = 1.; // Height of the channel
  double L = 5.; // Length of the channel
  double dpdx = -1.0;
  feVectorFunction uExact(solU, {H, dpdx});
  feFunction pExact(solP, {L, dpdx});

  feMesh2DP1 mesh(meshFile);

  feSpace *uDomain = nullptr, *pDomain = nullptr;
  feSpace *uInlet = nullptr, *uOutlet = nullptr, *uNoSlip = nullptr;

  // Poiseuille
  feCheck(createFiniteElementSpace( uInlet, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",  "Inlet", degreeQuadrature, &uExact));
  feCheck(createFiniteElementSpace(uOutlet, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Outlet", degreeQuadrature, &zeroVector));
  feCheck(createFiniteElementSpace(uNoSlip, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "NoSlip", degreeQuadrature, &zeroVector));
  feCheck(createFiniteElementSpace(uDomain, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Domain", degreeQuadrature, &zeroVector));
  feCheck(createFiniteElementSpace(pDomain, &mesh, elementType::LAGRANGE,        orderPressure, "P", "Domain", degreeQuadrature, &zero));

  std::vector<feSpace*> spaces = {uInlet, uNoSlip, uDomain, pDomain};
  std::vector<feSpace*> essentialSpaces = {uNoSlip, uInlet};

  if(divergenceFormulation) {
    // Solve for u-velocity on outlet but set v-velocity to zero
    // This in turn sets p = 0 on the outlet
    spaces.push_back(uOutlet);
    uOutlet->setEssentialComponent(1, true);
  } else {
    // Free BC on the outlet :
    // effectively sets p = 0 and dvdx = 0 on the outlet
  }

  feMetaNumber numbering(&mesh, spaces, essentialSpaces);
  feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

  // Continuity
  feBilinearForm *divU = nullptr;
  feCheck(createBilinearForm(divU, {pDomain, uDomain}, new feSysElm_MixedDivergence(&one)));
  std::vector<feBilinearForm*> forms = {divU};

  // Momentum
  feBilinearForm *diffU = nullptr, *gradP = nullptr, *divSigma = nullptr;
  if(divergenceFormulation) {
    feCheck(createBilinearForm(divSigma, {uDomain, pDomain}, new feSysElm_DivergenceNewtonianStress(&one, &viscosity)));
    forms.push_back(divSigma);
   } else {
    feCheck(createBilinearForm(gradP, {uDomain, pDomain}, new feSysElm_MixedGradient(&minusOne)));
    feCheck(createBilinearForm(diffU,           {uDomain}, new feSysElm_VectorDiffusion(&minusOne, &viscosity)));
    forms.push_back(gradP);
    forms.push_back(diffU);
  }

  // Nonlinear term (optional)
  feBilinearForm *convU = nullptr;
  feCheck(createBilinearForm(convU, {uDomain}, new feSysElm_VectorConvectiveAcceleration(&one)));
  forms.push_back(convU);

  feLinearSystem *system;
  if(std::string(solverType).compare("PETSc") == 0) {
    feCheck(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
  } else if(std::string(solverType).compare("PETSc_MUMPS") == 0) {
    feCheck(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
  } else if(std::string(solverType).compare("Pardiso") == 0) {
    feCheck(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
  } else {
    feErrorMsg(FE_STATUS_ERROR, "Unrecognized linear solver type: %s", solverType);
    return 1;
  }

  // system->setDisplayMatrixInConsole(true);
  // system->setDisplayRHSInConsole(true);

  feExporter *exporter;
  feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
  // int exportEveryNSteps = 1;
  std::string vtkFileRoot = "poiseuille";
  feExportData exportData = {exporter, 1, vtkFileRoot};

  std::vector<feNorm*> norms = {};

  TimeIntegrator *solver;
  feNLSolverOptions tol{1e-14, 1e-14, 1e4, 20, 3, 1e-1};
  feCheck(createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, tol, system, &sol, &mesh, norms, exportData));
  feCheck(solver->makeSteps(1));

  feNorm normU(VECTOR_L2_ERROR, {uDomain}, &sol, nullptr, &uExact);
  feNorm normP(L2_ERROR, {pDomain}, &sol, &pExact);
  feInfo("Error on velocity is %1.14e", normU.compute());
  feInfo("Error on pressure is %1.14e", normP.compute());

  delete solver;
  delete system;
  delete exporter;
  for(feBilinearForm* f : forms)
    delete f;
  for(feSpace *s : spaces)
    delete s;

  finalize();
  return 0;
}