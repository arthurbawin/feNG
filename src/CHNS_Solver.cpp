
#include "CHNS_Solver.h"

#include "TransientAdapter.h"
#include "feTimeIntegration.h"

feStatus
createCHNS_Solver(CHNS_Solver                      *&solver,
                  const Parameters::CHNS            &CHNS_parameters,
                  const Parameters::NonLinearSolver &NLSolver_parameters,
                  const Parameters::TimeIntegration &TimeIntegration_parameters,
                  const std::vector<FEDescriptor>   &fieldDescriptors,
                  const std::vector<BoundaryCondition> &boundaryConditions,
                  PostProcCallback                     *postProcCallback,
                  feFunction                           *pSource,
                  feVectorFunction                     *uSource,
                  feFunction                           *phiSource,
                  feFunction                           *muSource)
{
  if (fieldDescriptors.size() != 4)
  {
    return feErrorMsg(FE_STATUS_ERROR,
                      "CHNS Solver expects 4 unknown finite element spaces: u, "
                      "p, phi, mu. Only %d field descriptor(s) were provided.",
                      fieldDescriptors.size());
  }

  std::string expectedName[4] = {"U", "P", "Phi", "Mu"};

  for (size_t i = 0; i < fieldDescriptors.size(); ++i)
  {
    if (fieldDescriptors[i]._fieldName != expectedName[i])
    {
      return feErrorMsg(FE_STATUS_ERROR,
                        "FE space %d: Expected field name: %s - Received: %s",
                        i,
                        expectedName[i].data(),
                        fieldDescriptors[i]._fieldName.data());
    }
  }

  solver = new CHNS_Solver(CHNS_parameters,
                           NLSolver_parameters,
                           TimeIntegration_parameters,
                           fieldDescriptors,
                           boundaryConditions,
                           postProcCallback,
                           pSource,
                           uSource,
                           phiSource,
                           muSource);

  return FE_STATUS_OK;
}

double density_f(const feFunctionArguments &args,
                 const std::vector<double> &par)
{
  // Assumes phi is in [-1,1]
  // Property of A when phi = 1, of B when phi = -1
  double phi  = args.u;
  double rhoA = par[0];
  double rhoB = par[1];
  return (rhoA - rhoB) / 2. * phi + (rhoA + rhoB) / 2.;
}

double densityLimiter_f(const feFunctionArguments &args,
                        const std::vector<double> &par)
{
  double phi  = fmax(-1., fmin(1., args.u));
  double rhoA = par[0];
  double rhoB = par[1];
  return (rhoA - rhoB) / 2. * phi + (rhoA + rhoB) / 2.;
}

double drhodphi_f(const feFunctionArguments & /*args*/,
                  const std::vector<double> &par)
{
  double rhoA = par[0];
  double rhoB = par[1];
  return (rhoA - rhoB) / 2.;
}

double viscosity_f(const feFunctionArguments &args,
                   const std::vector<double> &par)
{
  double phi   = args.u;
  double viscA = par[0];
  double viscB = par[1];
  return (viscA - viscB) / 2. * phi + (viscA + viscB) / 2.;
}

double viscosityLimiter_f(const feFunctionArguments &args,
                          const std::vector<double> &par)
{
  double phi   = fmax(-1., fmin(1., args.u));
  double viscA = par[0];
  double viscB = par[1];
  return (viscA - viscB) / 2. * phi + (viscA + viscB) / 2.;
}

double dviscdphi_f(const feFunctionArguments & /*args*/,
                   const std::vector<double> &par)
{
  double viscA = par[0];
  double viscB = par[1];
  return (viscA - viscB) / 2.;
}

double degenerateMobility_f(const feFunctionArguments &args,
                            const std::vector<double> &par)
{
  double       phi         = args.u;
  const double mobilityVal = par[0];
  return mobilityVal * fabs(1. - phi * phi);
}

CHNS_Solver::CHNS_Solver(
  const Parameters::CHNS               &CHNS_parameters,
  const Parameters::NonLinearSolver    &NLSolver_parameters,
  const Parameters::TimeIntegration    &TimeIntegration_parameters,
  const std::vector<FEDescriptor>      &fieldDescriptors,
  const std::vector<BoundaryCondition> &boundaryConditions,
  PostProcCallback                     *postProcCallback,
  feFunction                           *pSource,
  feVectorFunction                     *uSource,
  feFunction                           *phiSource,
  feFunction                           *muSource)
  : SolverBase(NLSolver_parameters,
               TimeIntegration_parameters,
               fieldDescriptors,
               boundaryConditions,
               postProcCallback)
  , _CHNS_parameters(CHNS_parameters)
  , _pSource(pSource)
  , _uSource(uSource)
  , _phiSource(phiSource)
  , _muSource(muSource)
{
  _numFields = 4;

  const double rhoA = CHNS_parameters.fluids[0].density;
  const double rhoB = CHNS_parameters.fluids[1].density;

  const double viscA = CHNS_parameters.fluids[0].dynamic_viscosity;
  const double viscB = CHNS_parameters.fluids[1].dynamic_viscosity;

  const double M = CHNS_parameters.mobility;

  //
  // Initialize CHNS callbacks
  //
  if (CHNS_parameters.limitPhaseMarker)
  {
    _density   = new feFunction(densityLimiter_f, {rhoA, rhoB});
    _viscosity = new feFunction(viscosityLimiter_f, {viscA, viscB});
  }
  else
  {
    _density   = new feFunction(density_f, {rhoA, rhoB});
    _viscosity = new feFunction(viscosity_f, {viscA, viscB});
  }

  _drhodphi  = new feFunction(drhodphi_f, {rhoA, rhoB});
  _dviscdphi = new feFunction(dviscdphi_f, {viscA, viscB});

  switch (CHNS_parameters.mobilityType)
  {
    case Parameters::CHNS::MobilityType::constant:
      _mobility = new feConstantFunction({M});
      break;
    case Parameters::CHNS::MobilityType::degenerate:
      _mobility = new feFunction(degenerateMobility_f, {M});
      break;
  }

  _volumeForce =
    new feConstantVectorFunction({0., CHNS_parameters.gravity, 0.});

  if (pSource == nullptr)
  {
    _pSource   = new feConstantFunction(0.);
    _uSource   = new feConstantVectorFunction({0., 0., 0.});
    _phiSource = new feConstantFunction(0.);
    _muSource  = new feConstantFunction(0.);
  }

  //
  // Initialize derived coefficients
  //
  _CHNS_parameters.mass_alpha = (rhoB - rhoA) / (rhoA + rhoB);
}

feStatus
CHNS_Solver::createBilinearForms(const std::vector<feSpace *>  &spaces,
                                     std::vector<feBilinearForm *> &forms) const
{
  forms.clear();

  feSpace *u   = spaces[0];
  feSpace *p   = spaces[1];
  feSpace *phi = spaces[2];
  feSpace *mu  = spaces[3];

  std::vector<double> scalarParameters;

  feBilinearForm *CHNS;

  switch (_CHNS_parameters.formulation)
  {
    case Parameters::CHNS::Formulation::volumeAveraged:

      scalarParameters.push_back(_CHNS_parameters.surfaceTension);
      scalarParameters.push_back(_CHNS_parameters.epsilon);

      feCheckReturn(
        createBilinearForm(CHNS,
                           {u, p, phi, mu},
                           new CHNS_VolumeAveraged<2>(_density,
                                                      _drhodphi,
                                                      _viscosity,
                                                      _dviscdphi,
                                                      _mobility,
                                                      _volumeForce,
                                                      _pSource,
                                                      _uSource,
                                                      _phiSource,
                                                      _muSource,
                                                      scalarParameters)));
      break;

    case Parameters::CHNS::Formulation::massAveraged:

      scalarParameters.push_back(_CHNS_parameters.mass_alpha);
      scalarParameters.push_back(_CHNS_parameters.surfaceTension);
      scalarParameters.push_back(_CHNS_parameters.epsilon);

      feCheckReturn(
        createBilinearForm(CHNS,
                           {u, p, phi, mu},
                           new feSysElm_CHNS_Alternative<2>(_density,
                                                            _drhodphi,
                                                            _viscosity,
                                                            _dviscdphi,
                                                            _mobility,
                                                            _volumeForce,
                                                            _pSource,
                                                            _uSource,
                                                            _phiSource,
                                                            _muSource,
                                                            scalarParameters)));
      break;
  }

  forms = {CHNS};
  
  //
  // TODO: Loop over boundaries and impose flux BC if needed
  //

  return FE_STATUS_OK;
}

void CHNS_Solver::updateBeforeFixedPointIteration(const int iter)
{
  const Parameters::CHNS::MeshAdaptation &p = _CHNS_parameters.adaptation;

  // Reduce interface thickness
  if (iter > 0 && (iter % p.epsilonDecreaseFrequency) == 0)
  {
    const double prev = _CHNS_parameters.epsilon;
    _CHNS_parameters.epsilon /= p.epsilonDecrease;
    feInfo("CHNS epsilon was updated from %+-1.6e to %+-1.6e",
      prev, _CHNS_parameters.epsilon);
  }

  // // Update functions depending on epsilon or M
  const double M = 0.1 * _CHNS_parameters.epsilon * _CHNS_parameters.epsilon;
  _mobility->resetParameters({M});
}

//
// Compute error on scalar fields
//
feStatus CHNS_Solver::computeError(const std::vector<feSpace *> &spaces,
                                   feMesh                       *mesh,
                                   feSolution                   *sol,
                                   const int         iIntervalReferenceSolution,
                                   TransientAdapter &adapter,
                                   std::vector<double> &errors) const
{
  if (adapter._adapt_parameters.readReferenceSolution)
  {
    feInfo("Computing error for CHNS based on reference finer solution");

    errors.resize(4, 0.);

    // feSpace *u   = spaces[0];
    feSpace *p   = spaces[1];
    feSpace *phi = spaces[2];
    feSpace *mu  = spaces[3];

    feNorm *errorOnScalarSpaces = nullptr;
    feCheckReturn(createNorm(errorOnScalarSpaces, NONE, {p, phi, mu}, sol));
    feCheckReturn(errorOnScalarSpaces->computeErrorNormFromExternalSolution(
      adapter.referenceTestCases[iIntervalReferenceSolution].mesh,
      adapter.referenceTestCases[iIntervalReferenceSolution].sol,
      adapter.referenceTestCases[iIntervalReferenceSolution].spaces,
      mesh,
      sol,
      errors));
    feInfo("Computed error on P  : %+-1.6e", errors[0]);
    feInfo("Computed error on Phi: %+-1.6e", errors[1]);
    feInfo("Computed error on Mu : %+-1.6e", errors[2]);

    delete errorOnScalarSpaces;
  }
  else
  {
    feInfo("Not computing CHNS error as no reference solution was provided.");
  }
  return FE_STATUS_OK;
}