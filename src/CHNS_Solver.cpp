
#include "CHNS_Solver.h"

#include "TransientAdapter.h"
#include "feTimeIntegration.h"

feStatus createCHNS_Solver(
  CHNS_Solver                      *&solver,
  const Parameters::CHNS            &CHNS_parameters,
  const Parameters::NonLinearSolver &NLSolver_parameters,
  const Parameters::TimeIntegration &TimeIntegration_parameters,
  const std::vector<FEDescriptor>   &fieldDescriptors,
  const std::vector<const BoundaryConditions::BoundaryCondition *>
                   &boundaryConditions,
  PostProcCallback *postProcCallback,
  feFunction       *pSource,
  feVectorFunction *uSource,
  feFunction       *phiSource,
  feFunction       *muSource)
{
  size_t requiredDescriptors = 4;
  if (CHNS_parameters.formulation ==
      Parameters::CHNS::Formulation::volumeAveragedGeneric)
    requiredDescriptors = 5;

  if (fieldDescriptors.size() != requiredDescriptors)
  {
    return feErrorMsg(
      FE_STATUS_ERROR,
      "CHNS Solver expects %d unknown finite element spaces: u, "
      "p, phi, mu (,v). %d field descriptor(s) were provided.",
      requiredDescriptors,
      fieldDescriptors.size());
  }

  std::vector<std::string> expectedName = {"U", "P", "Phi", "Mu"};
  if (CHNS_parameters.formulation ==
      Parameters::CHNS::Formulation::volumeAveragedGeneric)
    expectedName.push_back("V");

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

  //
  // Test validity of boundary conditions
  //
  for (auto b : boundaryConditions)
  {
    if (b->_boundaryType == BoundaryConditions::Type::contactAngle)
    {
      const BoundaryConditions::FluxContactAngle *flux =
        dynamic_cast<const BoundaryConditions::FluxContactAngle *>(b);

      if (flux->_contactAngle <= 0. || flux->_contactAngle >= 180.)
      {
        return feErrorMsg(FE_STATUS_ERROR,
                          "Contact angle must be in ]0,180[ degrees.");
      }

      // Contact angle is only defined for the phase marker,
      // and the associated test space must be the potential.
      if (b->_descriptor_trialSpace._fieldName != "Phi")
      {
        return feErrorMsg(
          FE_STATUS_ERROR,
          "Contact angle condition for the Cahn-Hilliard Navier-Stokes system "
          "is only defined for the phase marker variable \"Phi\".");
      }
      if (b->_descriptor_testSpace._fieldName != "Mu")
      {
        return feErrorMsg(FE_STATUS_ERROR,
                          "The test FE space for the contact angle condition "
                          "for the Cahn-Hilliard Navier-Stokes system must be "
                          "the chemical potential \"Mu\"");
      }
    }
  }

  //
  // Formulation-specific constraints
  //
  if (CHNS_parameters.formulation ==
        Parameters::CHNS::Formulation::massAveraged &&
      CHNS_parameters.mobilityType == Parameters::CHNS::MobilityType::constant)
  {
    return feErrorMsg(
      FE_STATUS_ERROR,
      "Mobility was set to constant, but mass-averaged CHNS formulation is "
      "consistent only with degenerate mobility.");
  }

  if (CHNS_parameters.formulation == Parameters::CHNS::Formulation::khanwale &&
      TimeIntegration_parameters.method !=
        Parameters::TimeIntegration::TimeIntegrationMethod::bdf1)
  {
    return feErrorMsg(
      FE_STATUS_ERROR,
      "CHNS formulation from Khanwale et al. requires implicit Euler (BDF1) "
      "time stepping for consistency with energy stability properties.");
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
  const double phi  = args.u;
  const double rhoA = par[0];
  const double rhoB = par[1];
  return (rhoA - rhoB) / 2. * phi + (rhoA + rhoB) / 2.;
}

double densityLimiter_f(const feFunctionArguments &args,
                        const std::vector<double> &par)
{
  const double phi  = fmax(-1., fmin(1., args.u));
  const double rhoA = par[0];
  const double rhoB = par[1];
  return (rhoA - rhoB) / 2. * phi + (rhoA + rhoB) / 2.;
}

double drhodphi_f(const feFunctionArguments & /*args*/,
                  const std::vector<double> &par)
{
  const double rhoA = par[0];
  const double rhoB = par[1];
  return (rhoA - rhoB) / 2.;
}

template <int whichFluid, bool withLimiter>
double indicatorDensity_f(const feFunctionArguments &args,
                          const std::vector<double> &par)
{
  const double phi = args.u;
  double       density;
  if constexpr (withLimiter)
  {
    density = densityLimiter_f(args, par);
  }
  else
  {
    density = density_f(args, par);
  }
  if constexpr (whichFluid == 0)
  {
    // In fluid A: phi >= 0.
    return (phi >= 0) ? density : 0;
  }
  else
  {
    // In fluid B: phi <= 0
    return (phi <= 0) ? density : 0;
  }
}

double viscosity_f(const feFunctionArguments &args,
                   const std::vector<double> &par)
{
  const double phi   = args.u;
  const double viscA = par[0];
  const double viscB = par[1];
  return (viscA - viscB) / 2. * phi + (viscA + viscB) / 2.;
}

double viscosityLimiter_f(const feFunctionArguments &args,
                          const std::vector<double> &par)
{
  const double phi   = fmax(-1., fmin(1., args.u));
  const double viscA = par[0];
  const double viscB = par[1];
  return (viscA - viscB) / 2. * phi + (viscA + viscB) / 2.;
}

double dviscdphi_f(const feFunctionArguments & /*args*/,
                   const std::vector<double> &par)
{
  const double viscA = par[0];
  const double viscB = par[1];
  return (viscA - viscB) / 2.;
}

template <int whichFluid, bool withLimiter>
double indicatorViscosity_f(const feFunctionArguments &args,
                            const std::vector<double> &par)
{
  const double phi = args.u;
  double       viscosity;
  if constexpr (withLimiter)
  {
    viscosity = viscosityLimiter_f(args, par);
  }
  else
  {
    viscosity = viscosity_f(args, par);
  }
  if constexpr (whichFluid == 0)
  {
    // In fluid A: phi >= 0.
    return (phi >= 0) ? viscosity : 0;
  }
  else
  {
    // In fluid B: phi <= 0
    return (phi <= 0) ? viscosity : 0;
  }
}

double degenerateMobility_f(const feFunctionArguments &args,
                            const std::vector<double> &par)
{
  const double phi         = args.u;
  const double mobilityVal = par[0];
  return mobilityVal * fabs(1. - phi * phi);
}

CHNS_Solver::CHNS_Solver(
  const Parameters::CHNS            &CHNS_parameters,
  const Parameters::NonLinearSolver &NLSolver_parameters,
  const Parameters::TimeIntegration &TimeIntegration_parameters,
  const std::vector<FEDescriptor>   &fieldDescriptors,
  const std::vector<const BoundaryConditions::BoundaryCondition *>
                   &boundaryConditions,
  PostProcCallback *postProcCallback,
  feFunction       *pSource,
  feVectorFunction *uSource,
  feFunction       *phiSource,
  feFunction       *muSource)
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
  bool limiter = CHNS_parameters.limitPhaseMarker;
  if (limiter)
  {
    _density   = new feFunction(densityLimiter_f, {rhoA, rhoB});
    _viscosity = new feFunction(viscosityLimiter_f, {viscA, viscB});

    _densityIndicatorFluidA =
      new feFunction(indicatorDensity_f<0, true>, {rhoA, rhoB});
    _densityIndicatorFluidB =
      new feFunction(indicatorDensity_f<1, true>, {rhoA, rhoB});
    _viscosityIndicatorFluidA =
      new feFunction(indicatorViscosity_f<0, true>, {rhoA, rhoB});
    _viscosityIndicatorFluidB =
      new feFunction(indicatorViscosity_f<1, true>, {rhoA, rhoB});
  }
  else
  {
    _density   = new feFunction(density_f, {rhoA, rhoB});
    _viscosity = new feFunction(viscosity_f, {viscA, viscB});

    _densityIndicatorFluidA =
      new feFunction(indicatorDensity_f<0, false>, {rhoA, rhoB});
    _densityIndicatorFluidB =
      new feFunction(indicatorDensity_f<1, false>, {rhoA, rhoB});
    _viscosityIndicatorFluidA =
      new feFunction(indicatorViscosity_f<0, false>, {rhoA, rhoB});
    _viscosityIndicatorFluidB =
      new feFunction(indicatorViscosity_f<1, false>, {rhoA, rhoB});
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

  const double lambda = 3. / (2. * sqrt(2.)) * _CHNS_parameters.surfaceTension *
                        _CHNS_parameters.epsilon;
  _coeffFluxContactAngle = new feConstantFunction(lambda);

  feInfo("CHNS SOLVER");
  feInfo("\tCreated solver with parameters:");
  feInfo("\t\t                   rhoA = %+-1.6e", rhoA);
  feInfo("\t\t                   rhoB = %+-1.6e", rhoB);
  feInfo("\t\t                    muA = %+-1.6e", viscA);
  feInfo("\t\t                    muB = %+-1.6e", viscB);
  feInfo("\t\t                epsilon = %+-1.6e", _CHNS_parameters.epsilon);
  feInfo("\t\t        surface tension = %+-1.6e",
         _CHNS_parameters.surfaceTension);
  feInfo("\t\t                gravity = %+-1.6e", _CHNS_parameters.gravity);
  feInfo("\t\t               mobility = %+-1.6e", _CHNS_parameters.mobility);
  feInfo("\t\t    degenerate mobility = %s",
         (CHNS_parameters.mobilityType ==
          Parameters::CHNS::MobilityType::degenerate) ?
           "yes" :
           "no");
  feInfo("\t\tlimiter on phase marker = %s",
         (CHNS_parameters.limitPhaseMarker) ? "yes" : "no");
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
    case Parameters::CHNS::Formulation::abels:

      scalarParameters.push_back(_CHNS_parameters.surfaceTension);
      scalarParameters.push_back(_CHNS_parameters.epsilon);

      feCheckReturn(createBilinearForm(CHNS,
                                       {u, p, phi, mu},
                                       new CHNS_Abels<2>(_density,
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
                           new CHNS_MassAveraged<2>(_density,
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

    case Parameters::CHNS::Formulation::khanwale:
    {
      const double rhoA = _CHNS_parameters.fluids[0].density;
      const double rhoB = _CHNS_parameters.fluids[1].density;
      
      scalarParameters.push_back(_CHNS_parameters.Re);
      scalarParameters.push_back(_CHNS_parameters.Pe);
      scalarParameters.push_back(_CHNS_parameters.Cn);
      scalarParameters.push_back(_CHNS_parameters.We);
      scalarParameters.push_back(_CHNS_parameters.Fr);
      scalarParameters.push_back(rhoA);
      scalarParameters.push_back(rhoB);

      feCheckReturn(
        createBilinearForm(CHNS,
                           {u, p, phi, mu},
                           new CHNS_Khanwale<2>(_density,
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
    case Parameters::CHNS::Formulation::volumeAveragedGeneric:
    {
      feSpace *v = spaces[4];

      scalarParameters.push_back(_CHNS_parameters.mass_alpha);
      scalarParameters.push_back(_CHNS_parameters.surfaceTension);
      scalarParameters.push_back(_CHNS_parameters.epsilon);

      feCheckReturn(createBilinearForm(
        CHNS,
        {u, p, phi, mu, v},
        new CHNS_VolumeAveragedGeneric<2>(_density,
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
  }

  forms = {CHNS};

  for (auto b : _boundaryConditions)
  {
    //
    // Contact angle boundary condition
    //
    if (b->_boundaryType == BoundaryConditions::Type::contactAngle)
    {
      const BoundaryConditions::FluxContactAngle *fluxBC =
        dynamic_cast<const BoundaryConditions::FluxContactAngle *>(b);

      // Find FE spaces for phase marker and potential
      feSpace *muContact = nullptr, *phiContact = nullptr;
      for (auto &s : spaces)
      {
        if (s->representsSameFieldAs(b->_descriptor_trialSpace))
        {
          phiContact = s;
        }
        if (s->representsSameFieldAs(b->_descriptor_testSpace))
        {
          muContact = s;
        }
      }
      if (!muContact || !phiContact)
      {
        return feErrorMsg(
          FE_STATUS_ERROR,
          "Could not create a FluxPrescribedAngle bilinear form because the FE "
          "space for Phi or Mu on the prescribed boundary could not be found "
          "in the vector of all FE spaces. This should not happen unless the "
          "flux boundary condition is ill-formed.");
      }

      feBilinearForm *flux;
      feCheck(createBilinearForm(flux,
                                 {muContact, phiContact, phi},
                                 new feSysElm_FluxPrescribedNormalAngle<2>(
                                   _coeffFluxContactAngle,
                                   fluxBC->_contactAngle * M_PI / 180.)));
      forms.push_back(flux);
    }
  }

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
           prev,
           _CHNS_parameters.epsilon);
  }

  // Update functions depending on epsilon or M
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