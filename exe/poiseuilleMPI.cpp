
#include "feAPI.h"
#if defined(HAVE_PETSC)
  #include "petscsys.h"
#endif
#include "mpi.h"

double solP(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  UNUSED(t, pos, par);
  // double x = pos[0];
  return 0.;
}

void solU(const double t, const std::vector<double> &pos, const std::vector<double> &par, std::vector<double> &res)
{
  UNUSED(t);
  double y = pos[1];
  double H = par[0];
  double dpdx = par[1];

  res[0] = -dpdx/2. * y * (H - y);
  res[1] = 0.;
}

int main(int argc, char **argv)
{
  initialize(argc, argv);

  const char* meshFile = "../data/poiseuille0.msh";
  int orderVelocity = 2;
  int orderPressure = 1;
  int degreeQuadrature = 8;
  const char *solverType = "Pardiso";
  bool divergenceFormulation = false;

#if defined(HAVE_PETSC)
  // The options parser conflicts with runtime PETSc options
  char meshFileTmp[64];
  PetscBool permute;
  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "Poiseuille example options", "");
  {
    PetscCall(PetscOptionsString("-m", "Mesh file", "", meshFileTmp, meshFileTmp, sizeof(meshFileTmp), NULL));
    permute  = PETSC_FALSE;
    PetscCall(PetscOptionsBool("-permute", "Reorder the matrix?", "", permute, &permute, NULL));
  }
  PetscOptionsEnd();
  meshFile = &meshFileTmp[0];
#else
  // Create an option parser and parse the command line arguments.
  // The parser is not compatible with the PETSc options for now.
  // Returns an error if a command line argument is ill-formed.
  // If a command line argument is provided, it will overwrite the default parameter.
  // Each parameter added with "addOption" is optional by default, to make it required,
  // set the 5th argument of addOption to "true".
  feOptionsParser options(argc, argv);
  options.addOption(&meshFile, "-m", "--mesh", "Mesh file (gmsh format .msh only)");
  options.addOption(&orderVelocity, "-o", "--order_velocity", "Polynomial order for the velocity");
  options.addOption(&solverType, "-linsol", "--linear_solver_type", "Choice of linear solver (PETSc or MKL Pardiso)");
  // options.addOption(&orderPressure, "-op", "--order_pressure", "Polynomial order for the pressure");
  feCheck(options.parse());
#endif

  orderPressure = fmin(1, orderVelocity-1);

  setVerbose(1);

  feConstantFunction zero(0.);
  feConstantFunction one(1.);
  feConstantFunction minusOne(-1.);
  feConstantVectorFunction zeroVector({0., 0.});
  feConstantVectorFunction  oneVector({1., 1.});
  feConstantVectorFunction  uniform({1., 0.});
  feConstantVectorFunction  traction({0., 0.});

  double Re = 1.;

  feConstantFunction viscosity(1./Re);
  feConstantFunction density(1.);
  feConstantVectorFunction volumeForce({0., 0.});

  double H = 1.;
  double dpdx = -1.0;
  feVectorFunction uExact(solU, {H, dpdx});

  feMesh2DP1 mesh(meshFile);

  feSpace *uDomaine = nullptr, *pDomaine = nullptr;
  feSpace *uInlet = nullptr, *uOutlet = nullptr, *uNoSlip = nullptr; //, *pPoint = nullptr;
  feSpace *phi = nullptr, *mu = nullptr;

  // Poiseuille
  feCheck(createFiniteElementSpace(uInlet, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Entree", degreeQuadrature, &uExact));
  // feCheck(createFiniteElementSpace(uInlet, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Entree", degreeQuadrature, &uniform));
  feCheck(createFiniteElementSpace(uOutlet, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Sortie", degreeQuadrature, &zeroVector));
  feCheck(createFiniteElementSpace(uNoSlip, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "NoSlip", degreeQuadrature, &zeroVector));
  feCheck(createFiniteElementSpace(uDomaine, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Domaine", degreeQuadrature, &zeroVector));
  feCheck(createFiniteElementSpace(pDomaine, &mesh, elementType::LAGRANGE, orderPressure, "P", "Domaine", degreeQuadrature, &zero));
  feCheck(createFiniteElementSpace(     phi, &mesh, elementType::LAGRANGE, orderPressure, "Phi", "Domaine", degreeQuadrature, &zero));
  feCheck(createFiniteElementSpace(      mu, &mesh, elementType::LAGRANGE, orderPressure,  "Mu", "Domaine", degreeQuadrature, &zero));

  std::vector<feSpace*> spaces = {uInlet, uNoSlip, uDomaine, pDomaine, phi, mu};
  std::vector<feSpace*> essentialSpaces = {uNoSlip, uInlet, phi, mu};

  if(divergenceFormulation) {
    // Solve for u-velocity on outlet but set v-velocity to zero
    spaces.push_back(uOutlet);
    uOutlet->setEssentialComponent(1, true);
  } else {
    // Free BC on the outlet
  }

  feMetaNumber numbering(&mesh, spaces, essentialSpaces);
  feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

  // Continuity
  feBilinearForm *divU = nullptr;
  feCheck(createBilinearForm( divU, {pDomaine, uDomaine}, new feSysElm_MixedDivergence(&one)));
  std::vector<feBilinearForm*> forms = {divU};

  // Momentum
  feBilinearForm *diffU = nullptr, *gradP = nullptr, *divSigma = nullptr;
  if(divergenceFormulation) {
    feCheck(createBilinearForm( divSigma, {uDomaine, pDomaine}, new feSysElm_DivergenceNewtonianStress(&one, &viscosity) ));
    forms.push_back(divSigma);
   } else {
    feCheck(createBilinearForm( gradP, {uDomaine, pDomaine}, new feSysElm_MixedGradient(&one)                    ));
    feCheck(createBilinearForm( diffU,           {uDomaine}, new feSysElm_VectorDiffusion(&minusOne, &viscosity) ));
    forms.push_back(gradP);
    forms.push_back(diffU);
  }

  // Nonlinear term (optional)
  feBilinearForm *convU = nullptr;
  feCheck(createBilinearForm( convU, {uDomaine}, new feSysElm_VectorConvectiveAcceleration(&one) ));
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

#if defined(HAVE_PETSC)
  system->setRelativeTol(1e-8);
  system->setReorderingStatus(permute);
#endif

  feExporter *exporter;
  feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
  // int exportEveryNSteps = 1;
  std::string vtkFileRoot = "vec";
  feExportData exportData = {exporter, 1, vtkFileRoot};
  // feExportData exportData = {nullptr, 1, ""};

  std::vector<feNorm*> norms = {};

  TimeIntegrator *solver;
  feTolerances tol{1e-14, 1e-14, 1e4, 20};
  feCheck(createTimeIntegrator(solver, STATIONARY, tol, system, &numbering, &sol, &mesh, norms, exportData));
  feCheck(solver->makeStep());

  // int rank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // if(rank == 0) {
  //   const std::vector<double> &mysol = sol.getSolutionReference();
  //   for(auto val : mysol) {
  //     feInfo("Sol = %f", val);
  //   }
  // }

  feNorm norm(VECTOR_L2_ERROR, {uDomaine}, &sol, nullptr, &uExact);
  feInfo("Error is %1.14e", norm.compute());

  // delete solver;
  delete system;
  delete exporter;
  // delete diffU;
  // delete divU;
  // delete gradP;
  // delete divSigma;
  delete uDomaine;
  delete pDomaine;
  delete uInlet;
  delete uOutlet;
  delete uNoSlip;

  finalize();
  return 0;
}