
#include "feAPI.h"
#include <numeric>
#include "CHNS_Solver.h"

//
// User data to handle post processing
// for the various benchmarks
//
// class PostProcData
// {
// public:
//   // A vector of post-processing results, which can be anything. Even a boat.
//   std::vector<double> res;

//   // The 0 isoline for CHNS test cases
//   std::vector<double> isoline;

// public:
//   PostProcData(){};
// };

class AdaptData
{
public:
  std::string readDir, writeDir, initialMesh;
  bool restartFromAdaptedMeshes = false;

  bool adapt = true;
  bool exportSolution = true;
  int exportFrequency = 5;

  int nVertices = 5000; // On each spatial mesh (!= space-time complexity)
  double p = 4.; // Lp norm
  double dim = 2.;

  // Adapt with the usual metric for steady solution if true
  bool steadyMetric = false;

  int nFixedPoint = 1;
  int iFixedPoint = 0;
  int nIntervals = 4;

  bool enableGradation = false;
  double gradation = 1.6;

  // Time integration
  double t0 = 0.;
  double tEnd = 0.1;
  double currentTime = t0;
  int currentStep = 0;
  int nTimeStepsPerInterval = 4;
  double dt = (tEnd - t0) / (nIntervals * nTimeStepsPerInterval);

  std::vector<feMesh2DP1*> allMeshes;
  std::vector<feMetric*> allMetrics;
  std::vector<feMetricOptions> allOptions;
  std::vector<std::map<int, MetricTensor>> allHi;
  std::vector<feSolutionContainer*> allContainers;
  std::vector<feMetaNumber*> allNumberings;
  std::vector<std::vector<feSpace*>> allSpaces;
  std::vector<std::vector<feSpace*>> allEssentialSpaces;

  struct referenceTestCaseStruct {
    feMesh2DP1 *mesh;
    feMetaNumber *numbering;
    feSolution *sol;
    std::vector<feSpace*> spaces;
  };

  // A fine mesh solution for each interval, against which the error on coarser meshes is evaluated 
  std::vector<referenceTestCaseStruct> referenceTestCase;

  // A buffer container to use when transfering solution from one mesh to the next
  feSolutionContainer bufferContainer;

public:
  AdaptData(){};
  ~AdaptData()
  {
    for(auto *ptr : allMeshes) delete ptr;
    for(auto *ptr : allMetrics) delete ptr;
    for(auto *ptr : allContainers) delete ptr;
    for(auto *ptr : allNumberings) delete ptr;
    for(auto &vec : allSpaces)
      for(auto *ptr : vec) delete ptr;

    for(auto &refCase : referenceTestCase) {
      delete refCase.mesh;
      for(auto &s : refCase.spaces)
        delete s;
      delete refCase.numbering;
      delete refCase.sol;
    }
  }

  feStatus parseAndInitialize(int argc, char **argv,
                              const std::string &writeDirectory,
                              const std::string &initialMeshName,
                              const std::string &readDirectory,
                              const int orderForAdaptation,
                              const int nIntervalsArg = -1,
                              const int nVerticesArg  = -1,
                              const int nTimeStepsPerIntervalArg = -1,
                              const int nIntervalsReferenceSolution = 1)
  {
    writeDir = writeDirectory;
    initialMesh = initialMeshName;
    readDir = readDirectory;

    bool ignoreUnrecognizedOptions = true;
    feOptionsParser options(argc, argv, ignoreUnrecognizedOptions);
    options.addOption(&adapt, "-adapt", "--adapt", "-no-adapt", "--no-adapt", "");
    options.addOption(&restartFromAdaptedMeshes, "-restart", "--restart", "-no-restart", "--no-restart", "");
    options.addOption(&exportSolution, "-export", "--export", "-no-export", "--no-export", "");
    options.addOption(&exportFrequency, "-exportFrequency", "--exportFrequency", "");
    options.addOption(&nVertices, "-nV", "--nVertices", "");
    options.addOption(&nFixedPoint, "-nfp", "--nfp", "");
    options.addOption(&nIntervals, "-nI", "--nIntervals", "");
    options.addOption(&p, "-Lp", "--Lp", "");
    options.addOption(&gradation, "-grad", "--gradation", "");
    options.addOption(&t0, "-t0", "--t0", "");
    options.addOption(&tEnd, "-tEnd", "--tEnd", "");
    options.addOption(&nTimeStepsPerInterval, "-nT", "--nT", "");
    options.addOption(&steadyMetric, "-steadyMetric", "--sM", "-no_steadyMetric", "-noSM", "");
    feCheckReturn(options.parse());

    // Override (for convergence study)
    if(nIntervalsArg > 0)
      nIntervals = nIntervalsArg;
    if(nVerticesArg > 0)
      nVertices = nVerticesArg;
    if(nTimeStepsPerIntervalArg > 0)
      nTimeStepsPerInterval = nTimeStepsPerIntervalArg;

    currentTime = t0;
    dt = (tEnd - t0) / (nIntervals * nTimeStepsPerInterval);

    // Create structures for each time interval
    allMeshes.resize(nIntervals, nullptr);
    allMetrics.resize(nIntervals, nullptr);
    allHi.resize(nIntervals);
    allContainers.resize(nIntervals, nullptr);
    allNumberings.resize(nIntervals, nullptr);
    allSpaces.resize(nIntervals);
    allEssentialSpaces.resize(nIntervals);

    referenceTestCase.resize(nIntervalsReferenceSolution);

    // Metric options common to all meshes
    feMetricOptions metricOptions(initialMesh);
    metricOptions.polynomialDegree = orderForAdaptation;

    switch(orderForAdaptation)
    {
      case 1 : metricOptions.method = adaptationMethod::ANISO_P1; break;
      case 2 : metricOptions.method = adaptationMethod::ANISO_P2; break;
      default :
        metricOptions.method = adaptationMethod::ANISO_PN;
    }
    metricOptions.nTargetVertices = nVertices;
    metricOptions.LpNorm = p;
    metricOptions.hMin = 1e-10;
    metricOptions.hMax = 10.;
    metricOptions.enableGradation = false;
    metricOptions.gradation = gradation;
    metricOptions.gradationSpace = GradationSpace::Mixed;
    metricOptions.mmgInputMeshfile = writeDir + "inputMMG.msh";
    metricOptions.mmgOutputMeshfile = writeDir + "outputMMG.msh";
    metricOptions.debug = false;
    
    allOptions.resize(nIntervals, metricOptions);

    for(int iI = 0; iI < nIntervals; ++iI)
    {
      // Different adapted mesh name for each sub-interval
      allOptions[iI].adaptedMeshName = writeDir + "adapted_ifp_0_interval_" + std::to_string(iI) + ".msh";

      if(restartFromAdaptedMeshes)
      {
        initialMesh = readDir + "adapted_ifp_2_interval_" + std::to_string(iI) + ".msh";
        allOptions[iI].backgroundMeshfile = initialMesh;
      }

      allMeshes[iI] = new feMesh2DP1(initialMesh);

      feCheckReturn(createMetricField(allMetrics[iI], allOptions[iI], allMeshes[iI]));
      allMetrics[iI]->setMetricsToZero();
      allMetrics[iI]->getMetrics(allHi[iI]);
      allMetrics[iI]->setMetricsToIdentity();
    }

    return FE_STATUS_OK;
  }

  feStatus readAdaptedMeshes()
  {
    for(int iI = 0; iI < nIntervals; ++iI)
    {
      // Read adapted mesh
      feInfo("Deleting mesh %d", iI);
      delete allMeshes[iI];
      std::string meshFile = writeDir + "adapted_ifp_" + std::to_string(iFixedPoint-1) + "_interval_" + std::to_string(iI) + ".msh";
      allMeshes[iI] = new feMesh2DP1(meshFile);

      // Create new metric field
      allOptions[iI].backgroundMeshfile = meshFile;
      allOptions[iI].adaptedMeshName = writeDir + "adapted_ifp_" + std::to_string(iFixedPoint) + "_interval_" + std::to_string(iI) + ".msh";
      delete allMetrics[iI];
      allMetrics[iI] = nullptr;
      feCheck(createMetricField(allMetrics[iI], allOptions[iI], allMeshes[iI]));
      allHi[iI].clear();
      allMetrics[iI]->setMetricsToZero();
      allMetrics[iI]->getMetrics(allHi[iI]);
      allMetrics[iI]->setMetricsToIdentity();
    }
    return FE_STATUS_OK;
  }
};

static void printInt(const char *name, int val) {
  feInfo("%35s = %d", name, val);
} 
static void printDouble(const char *name, double val) {
  feInfo("%35s = %+-1.16e", name, val);
} 
static void printString(const char *name, std::string val) {
  feInfo("%35s = %s", name, val.data());
} 
#define PRINT_INT(name) printInt(#name, (name))
#define PRINT_DBL(name) printDouble(#name, (name))
#define PRINT_STR(name) printString(#name, (name))


//
// An abstract solver: creates FE spaces, solves and compute metric tensors
//
class GenericSolver
{
public:
  bool useAnalyticalSolution = false;
  bool _isCHNS = false;

public:
  virtual ~GenericSolver(){};

  virtual int getNumFields() = 0;
  virtual std::string getReadDir() = 0;
  virtual std::string getWriteDir() = 0;
  virtual std::string getInitialMesh() = 0;
  virtual int getOrderForAdaptation() = 0;

  // Is this a solver for the CHNS system?
  bool isCHNS() { return _isCHNS; };

  virtual feStatus parseAndCreateFunctions(int argc, char **argv) = 0;

  virtual void printParameters() {}

  virtual feStatus updateParameters(const int /*fixedPointIteration*/)
  {
    return FE_STATUS_OK;
  }

  virtual feStatus createSpaces(feMesh *mesh,
                                std::vector<feSpace*> &spaces,
                                std::vector<feSpace*> &essentialSpaces) = 0;

  feStatus projectSolution(feMesh2DP1 *currentMesh,
                           feMesh2DP1 *nextMesh,
                           feSolutionContainer *container,
                           feSolutionContainer &bufferContainer,
                           feSolutionContainer *&nextContainer,
                           feMetaNumber *numbering,
                           const std::vector<feSpace*> spaces,
                           const std::vector<feSpace*> essentialSpaces)
  {
    UNUSED(currentMesh, nextMesh, container, numbering, spaces, essentialSpaces);

    // Copy current container into the buffer
    bufferContainer = *container;

    // For checks : last BDF solutions
    std::vector<feSolution> currentSolutions(bufferContainer.getNbSol());
    for(int i = 0; i < bufferContainer.getNbSol(); ++i) {
      currentSolutions[i] = feSolution(bufferContainer, i, spaces, essentialSpaces);
    }

    // Projection
    std::vector<feSpace*> nextSpaces, nextEssentialSpaces;
    feCheckReturn(this->createSpaces(nextMesh, nextSpaces, nextEssentialSpaces));
    feMetaNumber nextNumbering(nextMesh, nextSpaces, nextEssentialSpaces);
    feCheckReturn(currentMesh->transfer(nextMesh, numbering, &nextNumbering, &bufferContainer, spaces, essentialSpaces, nextSpaces));

    delete nextContainer;
    nextContainer = new feSolutionContainer(bufferContainer);

    // Check projection (visualization) and integral of each projected solution in container
    for(int i = 0; i < bufferContainer.getNbSol(); ++i)
    {
      feInfo("Checking solution %d in bufferContainer :", i);
      feSolution projectedSolution(bufferContainer, i, spaces, essentialSpaces);
      // feExporter *nextExporter;
      // feCheck(createVisualizationExporter(nextExporter, VTK, &nextNumbering, &projectedSolution, nextMesh, nextSpaces));
      // feCheck(nextExporter->writeStep("afterTransfer" + std::to_string(i) + ".vtk"));
      // delete nextExporter;

      // Integral of highest dimensionsal fields before and after transfer
      for(size_t j = 0; j < spaces.size(); ++j)
      {
        if(spaces[j]->getDim() == 2) {
          feNorm *integral, *integralNextMesh, *intDot, *intDotNextMesh;
          feCheck(createNorm(integral, INTEGRAL, {spaces[j]}, &currentSolutions[i]));
          feCheck(createNorm(integralNextMesh, INTEGRAL, {nextSpaces[j]}, &projectedSolution));
          feInfo("Integral of field %s on current mesh = %+-1.6e", spaces[j]->getFieldID().data(), integral->compute());
          feInfo("Integral of field %s on next    mesh = %+-1.6e", nextSpaces[j]->getFieldID().data(), integralNextMesh->compute());
          feCheck(createNorm(intDot, INTEGRAL_DT, {spaces[j]}, &currentSolutions[i]));
          feCheck(createNorm(intDotNextMesh, INTEGRAL_DT, {nextSpaces[j]}, &projectedSolution));
          feInfo("Integral of field %s dot on current mesh = %+-1.6e", spaces[j]->getFieldID().data(), intDot->compute());
          feInfo("Integral of field %s dot on next    mesh = %+-1.6e", nextSpaces[j]->getFieldID().data(), intDotNextMesh->compute());
          delete integral;
          delete integralNextMesh;
          delete intDot;
          delete intDotNextMesh;
        }
      }
    }      

    for(feSpace *s : nextSpaces)
      delete s;

    return FE_STATUS_OK;
  }

  virtual feStatus solve(AdaptData &adapter,
                         const int iInterval,
                         feMesh *mesh,
                         feSolution *sol,
                         feMetaNumber *numbering,
                         std::vector<feSpace*> &spaces,
                         const int /*numUnknowns*/,
                         feExportData &exportData) = 0;

  feStatus computeAndIncrementMetric(const double currentTime,
                                     const int iStep,
                                     const int iInterval,
                                     AdaptData &adapter,
                                     feMetric *metricField,
                                     feNewRecovery *recoveredField)
  {
    metricField->setRecoveredFields({recoveredField});
    metricField->setCurrentTime(currentTime);
    metricField->setMetricScaling(false);
    feCheckReturn(metricField->computeMetrics());

    // Trapeze rule
    if(iStep == 0 || iStep == adapter.nTimeStepsPerInterval-1)
      metricField->addMetricsToOther(adapter.dt/2., adapter.allHi[iInterval]);
    else
      metricField->addMetricsToOther(adapter.dt   , adapter.allHi[iInterval]);

    return FE_STATUS_OK;
  }

  virtual feStatus computeLastMetric(AdaptData &adapter,
                                     const int iInterval,
                                     feMesh *mesh,
                                     feSolution *sol,
                                     feMetaNumber *numbering,
                                     std::vector<feSpace*> &spaces) = 0;

  virtual feStatus computeError(const int iIntervalReference,
                                AdaptData &adapter,
                                const std::vector<feSpace*> &spaces,
                                feMesh *mesh,
                                feSolution *sol,
                                std::vector<double> &errors) = 0;

  virtual feStatus computePostProcessing(AdaptData &adapter,
                                         const std::vector<feSpace*> &spaces,
                                         feMesh *mesh,
                                         feSolution *sol,
                                         PostProcData &data)
  {
    UNUSED(adapter, spaces, mesh, sol, data);
    return FE_STATUS_OK;
  };

  //
  // Create the structures for a reference test case
  // used as an ersatz of the exact solution.
  //
  feStatus readReferenceTestCase(const int iInterval,
                                 AdaptData &adapter,
                                 const std::string meshName,
                                 const std::string solutionFileName)
  {
    AdaptData::referenceTestCaseStruct &refCase = adapter.referenceTestCase[iInterval];

    // // Free existing solution
    // delete refCase.mesh;
    // for(auto &s : refCase.spaces) delete s;
    // delete refCase.numbering;
    // delete refCase.sol;

    refCase.mesh = new feMesh2DP1(meshName);
    std::vector<feSpace*> refEssentialSpaces;
    this->createSpaces(refCase.mesh, refCase.spaces, refEssentialSpaces);
    refCase.numbering = new feMetaNumber(refCase.mesh, refCase.spaces, refEssentialSpaces);


    std::filebuf fb;
    if(fb.open(solutionFileName, std::ios::in)) {
      fb.close();

      // Read solution file if it exists
      refCase.sol = new feSolution(refCase.numbering->getNbDOFs(), solutionFileName);

    } else {

      return feErrorMsg(FE_STATUS_ERROR, "Could not create reference test case because file \"%s\" does not exist.",
        solutionFileName.data());
      
    }
    return FE_STATUS_OK;
  }
};

//
// Generic solver for CHNS system with slip BC
//
namespace CHNS_Generic 
{
  double rho_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    // Assumes phi is in [-1,1]
    // Property of A when phi = 1, of B when phi = -1
    double phi = fmax(-1., fmin(1., args.u));
    // double phi = args.u;
    double rhoA = par[0];
    double rhoB = par[1];
    return (rhoA - rhoB)/2. * phi + (rhoA + rhoB)/2.;
  }

  double drhodphi_f(const feFunctionArguments &/*args*/, const std::vector<double> &par)
  {
    double rhoA = par[0];
    double rhoB = par[1];
    return (rhoA - rhoB)/2.;
  }

  double viscosity_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    // Assumes phi is in [-1,1]
    // Property of A when phi = 1, of B when phi = -1
    double phi = fmax(-1., fmin(1., args.u));
    // double phi = args.u;
    double viscA = par[0];
    double viscB = par[1];
    return (viscA - viscB)/2. * phi + (viscA + viscB)/2.;
  }

  double dviscdphi_f(const feFunctionArguments &/*args*/, const std::vector<double> &par)
  {
    double viscA = par[0];
    double viscB = par[1];
    return (viscA - viscB)/2.;
  }

  //
  // Degenerate mobility (0 in pure phases)
  //
  double mobility_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    // double phi = fmax(-1., fmin(1., args.u));
    double phi = args.u;
    const double mobilityVal = par[0];
    return mobilityVal * fabs(1. - phi*phi);
  }

  double mobilityAlpha_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    // double phi = fmax(-1., fmin(1., args.u));
    double phi = args.u;
    const double mobilityVal = par[0];
    const double alpha       = par[1];
    return alpha * mobilityVal * fabs(1. - phi*phi);
  }

  double mobilityAlphaSquared_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    // double phi = fmax(-1., fmin(1., args.u));
    double phi = args.u;
    const double mobilityVal = par[0];
    const double alpha       = par[1];
    return alpha * alpha * mobilityVal * fabs(1. - phi*phi);
  }

  class GenericSolverCHNS : public GenericSolver
  {
  protected:
    std::string _readDir = "";
    std::string _writeDir = "";
    const char *_initialMesh = "";

    int _orderVelocity    = 2;
    int _orderPressure    = 1;
    int _orderPhi         = 1;
    int _orderMu          = 1;
    int _degreeQuadrature = 7;

    bool _allSlip = true;
    bool _withSUPG = false;

    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 50, 4, 1.};

    // Physical properties
    double _rhoA, _rhoB, _viscA, _viscB, _gamma, _epsilon, _mobilityVal;

    // Generic CHNS functions
    feConstantFunction *coeffKorteweg, *coeffPhiCube, *coeffPhi, *coeffDiffPhi;
    feFunction *density, *drhodphi, *viscosity, *dviscdphi, *mobility;

    // Functions initialized by the derived class
    feFunction *phiInit = nullptr;
    feVectorFunction *volumeForce = nullptr;

  public:
    GenericSolverCHNS()
    {
      _isCHNS = true; 
    };
    virtual ~GenericSolverCHNS(){};

    int getNumFields() { return 4; }
    std::string getReadDir(){ return _readDir; };
    std::string getWriteDir(){ return _writeDir; };
    std::string getInitialMesh(){ return std::string(_initialMesh); };
    int getOrderForAdaptation() { return _orderPhi; };

    // GenericSolverCHNS
    void printParameters() override
    {
      feInfo("CHNS Solver (with Zeineb's parameters) parameters :");
      PRINT_STR(_readDir);
      PRINT_STR(_writeDir);
      PRINT_STR(std::string(_initialMesh));

      PRINT_INT(_orderVelocity);
      PRINT_INT(_orderPressure);
      PRINT_INT(_orderPhi);
      PRINT_INT(_orderMu);

      PRINT_DBL(_rhoA);
      PRINT_DBL(_rhoB);
      PRINT_DBL(_viscA);
      PRINT_DBL(_viscB);
      PRINT_DBL(_gamma);
      PRINT_DBL(_epsilon);
      PRINT_DBL(_mobilityVal);
    }

    // GenericSolverCHNS
    feStatus parseGenericCHNSParameters(int argc, char **argv)
    {
      const char* writeDirectory = _writeDir.data();

      bool ignoreUnrecognizedOptions = true;
      feOptionsParser options(argc, argv, ignoreUnrecognizedOptions);
      options.addOption(&_initialMesh, "-m", "--m", "");
      options.addOption(&writeDirectory, "-dir", "--dir", "");
      options.addOption(&_orderVelocity, "-ou", "--orderU", "");
      options.addOption(&_orderPressure, "-op", "--orderP", "");
      options.addOption(&_orderPhi, "-ophi", "--orderPhi", "");
      options.addOption(&_orderMu , "-omu", "--orderMu", "");
      options.addOption(&_rhoA, "-rhoA", "--rhoA", "");
      options.addOption(&_rhoB, "-rhoB", "--rhoB", "");
      options.addOption(&_viscA, "-viscA", "--viscosityA", "");
      options.addOption(&_viscB, "-viscB", "--viscosityB", "");
      options.addOption(&_gamma, "-g", "--gamma", "Surface tension");
      options.addOption(&_epsilon, "-e", "--epsilon", "Interface thickness");
      options.addOption(&_mobilityVal, "-mob", "--mobility", "");
      options.addOption(&_allSlip, "-allSlip", "--allSlip", "-no-allSlip", "--no-allSlip", "");
      options.addOption(&_withSUPG, "-supg", "--supg", "-no-supg", "--no-supg", "");
      feCheckReturn(options.parse());

      _writeDir = std::string(writeDirectory);
      feInfo("Write directory is : %s", _writeDir.data());

      return FE_STATUS_OK;
    }

    // GenericSolverCHNS
    // Create functions once all parameters have been parsed
    feStatus createCHNSFunctions()
    {
      // Mémoire de Zeineb:
      coeffKorteweg = new feConstantFunction(-_gamma/_epsilon);
      coeffPhiCube  = new feConstantFunction(-1.);
      coeffPhi      = new feConstantFunction(+1.);
      coeffDiffPhi  = new feConstantFunction(-_epsilon*_epsilon);

      // Mémoire de Pierre:
      // const double lambda = 3. * _epsilon * _gamma / (2.*sqrt(2.));
      // coeffKorteweg = new feConstantFunction(-1.);
      // coeffPhiCube  = new feConstantFunction(- (_epsilon * _epsilon) / lambda);
      // coeffPhi      = new feConstantFunction(+ (_epsilon * _epsilon) / lambda);
      // coeffDiffPhi  = new feConstantFunction(- lambda);

      mobility      = new feFunction(mobility_f, {_mobilityVal});
      density       = new feFunction(     rho_f, {_rhoA, _rhoB});
      drhodphi      = new feFunction(drhodphi_f, {_rhoA, _rhoB});
      viscosity     = new feFunction(viscosity_f, {_viscA, _viscB});
      dviscdphi     = new feFunction(dviscdphi_f, {_viscA, _viscB});

      return FE_STATUS_OK;
    }

    // GenericSolverCHNS
    feStatus createSpaces(feMesh *mesh,
                          std::vector<feSpace*> &spaces,
                          std::vector<feSpace*> &essentialSpaces)
    {
      for(auto *s : spaces) delete s;
      spaces.clear();
      essentialSpaces.clear();

      feSpace *u = nullptr;
      feSpace *uNoSlip = nullptr;
      feSpace *uSlipX = nullptr, *uSlipY = nullptr;
      feSpace *p = nullptr, *pPoint = nullptr;
      feSpace *phi = nullptr, *mu = nullptr;

      feCheckReturn(createFiniteElementSpace(      u, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "Domaine", _degreeQuadrature, &vectorConstant::zero));

      if(_allSlip) {
        // Slip on all boundaries
        feCheckReturn(createFiniteElementSpace(uSlipX, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "SlipX", _degreeQuadrature, &vectorConstant::zero));
        feCheckReturn(createFiniteElementSpace(uSlipY, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "SlipY", _degreeQuadrature, &vectorConstant::zero));
      } else {
        // Slip on sides and no-slip top/bottom
        feCheckReturn(createFiniteElementSpace(uNoSlip, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "NoSlip", _degreeQuadrature, &vectorConstant::zero));
        feCheckReturn(createFiniteElementSpace( uSlipY, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",  "SlipY", _degreeQuadrature, &vectorConstant::zero));
      }

      feCheckReturn(createFiniteElementSpace(      p, mesh, elementType::LAGRANGE, _orderPressure,   "P",       "Domaine", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace( pPoint, mesh, elementType::LAGRANGE, _orderPressure,   "P", "PressurePoint", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(    phi, mesh, elementType::LAGRANGE, _orderPhi,        "Phi",     "Domaine", _degreeQuadrature, phiInit));
      feCheckReturn(createFiniteElementSpace(     mu, mesh, elementType::LAGRANGE, _orderMu,         "Mu",      "Domaine", _degreeQuadrature, &scalarConstant::zero));

      // uSlip->setEssentialComponent(0, true);
      // spaces = {u, uNoSlip, uSlip, p, pPoint, phi, mu};
      // essentialSpaces = {uNoSlip, pPoint};

      if(_allSlip) {
        uSlipY->setEssentialComponent(0, true);
        uSlipX->setEssentialComponent(1, true);
        spaces = {u, p, phi, mu, uSlipX, uSlipY, pPoint};
        essentialSpaces = {pPoint};
      } else {
        uSlipY->setEssentialComponent(0, true);
        spaces = {u, p, phi, mu, uNoSlip, uSlipY, pPoint};
        essentialSpaces = {uNoSlip, pPoint};
      }

      return FE_STATUS_OK;
    }

    // GenericSolverCHNS
    feStatus solve(AdaptData &adapter,
                   const int iInterval,
                   feMesh *mesh,
                   feSolution *sol,
                   feMetaNumber *numbering,
                   std::vector<feSpace*> &spaces,
                   const int /*numUnknowns*/,
                   feExportData &exportData)
    {
      std::vector<feBilinearForm*> forms = {};

      feSpace *u   = spaces[0];
      feSpace *p   = spaces[1];
      feSpace *phi = spaces[2];
      feSpace *mu  = spaces[3];

      // Momentum
      feBilinearForm *CHNS_momentum = nullptr;
      feCheckReturn(createBilinearForm(CHNS_momentum, {u, p, phi, mu}, new feSysElm_CHNS_Momentum<2>(density,
                                                                                                  drhodphi,
                                                                                                  viscosity,
                                                                                                  dviscdphi,
                                                                                                  mobility,
                                                                                                  coeffKorteweg,
                                                                                                  volumeForce)));
      forms.push_back(CHNS_momentum);

      // Continuity
      feBilinearForm *divU = nullptr;
      feCheckReturn(createBilinearForm(divU, {p, u}, new feSysElm_MixedDivergence<2>(&scalarConstant::one)));
      forms.push_back(divU);

      // Tracer convection
      feBilinearForm *dphidt = nullptr, *convPhi = nullptr, *diffMu = nullptr;
      feCheckReturn(createBilinearForm( dphidt,     {phi}, new feSysElm_TransientMass(&scalarConstant::one)));
      feCheckReturn(createBilinearForm(convPhi,  {phi, u}, new feSysElm_TracerConvection<2>(&scalarConstant::one)));
      feCheckReturn(createBilinearForm( diffMu, {phi, mu}, new feSysElm_MixedGradGrad<2>(mobility)));
      forms.push_back(dphidt);
      forms.push_back(convPhi);
      forms.push_back(diffMu);

      // SUPG for tracer
      if(_withSUPG) {
        feBilinearForm *supgTracer = nullptr;
        feCheckReturn(createBilinearForm(supgTracer, {phi, u, mu}, new feSysElm_CHNS_Tracer_SUPG<2>(&scalarConstant::one, mobility)));
        forms.push_back(supgTracer);
      }

      // Chemical potential
      feBilinearForm *muMass = nullptr, *phiCube = nullptr, *phiMass = nullptr, *diffPhi = nullptr;
      feCheckReturn(createBilinearForm( muMass,      {mu}, new feSysElm_Mass(&scalarConstant::one)));
      feCheckReturn(createBilinearForm(phiCube, {mu, phi}, new feSysElm_MixedMassPower(coeffPhiCube, 3)));
      feCheckReturn(createBilinearForm(phiMass, {mu, phi}, new feSysElm_MixedMass(coeffPhi)));
      feCheckReturn(createBilinearForm(diffPhi, {mu, phi}, new feSysElm_MixedGradGrad<2>(coeffDiffPhi)));
      forms.push_back(muMass);
      forms.push_back(phiCube);
      forms.push_back(phiMass);
      forms.push_back(diffPhi);

      for(feBilinearForm *form : forms) {
        if(form->getWeakFormName() != "CHNS_TRACER_SUPG")
        form->setComputeMatrixWithFD(false);
      }
      CHNS_momentum->setComputeMatrixWithFD(true);
      diffMu->setComputeMatrixWithFD(true);

      feLinearSystem *system = nullptr;
      #if defined(HAVE_MKL)
        feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
        feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering));
      #else
        feCheckReturn(createLinearSystem(system, PETSC, forms, numbering));
      #endif

      std::vector<feNorm*> norms = {};

      TimeIntegrator *timeIntegrator = nullptr;
      const double t0 = adapter.currentTime;
      const double t1 = adapter.currentTime + adapter.dt * adapter.nTimeStepsPerInterval;
      const int nTimeSteps = adapter.nTimeStepsPerInterval;
      feCheckReturn(createTimeIntegrator(timeIntegrator, timeIntegratorScheme::BDF2, NLoptions,
        system, sol, mesh, norms, exportData, t0, t1, nTimeSteps));

      // Restart from solution container if not the first sub-interval
      // (restart from last solutions of previous sub-interval)
      if(iInterval > 0) {
        // BDF1Integrator *ptr = dynamic_cast<BDF1Integrator*> (timeIntegrator);
        BDF2Integrator *ptr = dynamic_cast<BDF2Integrator*> (timeIntegrator);
        
        ptr->_sC->NaNify();
        feCheckReturn(timeIntegrator->restartFromContainer(*adapter.allContainers[iInterval]));
        timeIntegrator->setCurrentStep(adapter.currentStep);

        feInfo("Time           of timeIntegrator : %f", timeIntegrator->getCurrentTime());
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[0]);
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[1]);
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[2]);
        feInfo("Step           of timeIntegrator : %d", timeIntegrator->getCurrentStep());
      }

      // Solve and compute metrics at each time step
      for(int iStep = 0; iStep < adapter.nTimeStepsPerInterval; ++iStep)
      {
        tic();
        feCheckReturn(timeIntegrator->makeSteps(1));
        feInfo("Computed time step in %f s", toc());
        // sol->setCurrentTime(t0 + adapter.dt);
        // sol->initialize(mesh);

        // Reconstruct phi
        std::string meshFile = adapter.allOptions[iInterval].backgroundMeshfile;
        std::string recoveryFile = adapter.writeDir + "phi.msh";
        feNewRecovery recPhi(phi, 0, mesh, sol, meshFile, recoveryFile, false, false, false, nullptr, numbering);

        // Compute metric field and increment metric integral Hi
        feMetric *metricField = adapter.allMetrics[iInterval];
        computeAndIncrementMetric(timeIntegrator->getCurrentTime(), iStep, iInterval, adapter, metricField, &recPhi);

        adapter.currentStep++;
      }

      adapter.currentTime = timeIntegrator->getCurrentTime();
      // currentStep = timeIntegrator->getCurrentStep();

      delete adapter.allContainers[iInterval];
      adapter.allContainers[iInterval] = new feSolutionContainer();
      adapter.allContainers[iInterval]->NaNify();
      *(adapter.allContainers[iInterval]) = timeIntegrator->getSolutionContainer();

      delete timeIntegrator;
      delete system;
      for(feBilinearForm *f : forms)
        delete f;

      return FE_STATUS_OK;
    }

    // GenericSolverCHNS
    feStatus computeLastMetric(AdaptData &adapter,
                               const int iInterval,
                               feMesh *mesh,
                               feSolution *sol,
                               feMetaNumber *numbering,
                               std::vector<feSpace*> &spaces)
    {
      feSpace *phi = spaces[2];

      // Reconstruct velocity
      std::string meshFile = adapter.allOptions[iInterval].backgroundMeshfile;
      std::string recoveryFile = adapter.writeDir + "phi.msh";
      feNewRecovery rec(phi, 0, mesh, sol, meshFile, recoveryFile, false, false, false, nullptr, numbering);

      // Compute "steady state" metric field
      feMetric *metricField = adapter.allMetrics[iInterval];
      metricField->setRecoveredFields({&rec});
      metricField->_options.debug = true;
      metricField->setMetricScaling(true);
      feCheckReturn(metricField->computeMetrics());

      return FE_STATUS_OK;
    }

    //
    // Compute error on scalar fields
    //
    feStatus computeError(const int iIntervalReference,
                          AdaptData &adapter,
                          const std::vector<feSpace*> &spaces,
                          feMesh *mesh,
                          feSolution *sol,
                          std::vector<double> &errors)
    {
      feInfo("Computing error for CHNS based on reference finer solution");

      errors.resize(this->getNumFields(), 0.);

      // feSpace *u   = spaces[0];
      feSpace *p   = spaces[1];
      feSpace *phi = spaces[2];
      feSpace *mu  = spaces[3];

      feNorm *errorOnScalarSpaces = nullptr;
      feCheckReturn(createNorm(errorOnScalarSpaces, NONE, {p, phi, mu}, sol));
      feCheckReturn(errorOnScalarSpaces->computeErrorNormFromExternalSolution(adapter.referenceTestCase[iIntervalReference].mesh,
                                                                              adapter.referenceTestCase[iIntervalReference].sol,
                                                                              adapter.referenceTestCase[iIntervalReference].spaces,
                                                                              mesh, sol, errors));
      feInfo("Computed error on P  : %+-1.6e", errors[0]);
      feInfo("Computed error on Phi: %+-1.6e", errors[1]);
      feInfo("Computed error on Mu : %+-1.6e", errors[2]);

      delete errorOnScalarSpaces;
      return FE_STATUS_OK;
    }
  }; // GenericSolverCHNS

  class SolverCHNS_VolumeAveraged : public GenericSolver
  {
  protected:
    std::string _readDir = "";
    std::string _writeDir = "";
    const char *_initialMesh = "";

    int _orderVelocity    = 2;
    int _orderPressure    = 1;
    int _orderPhi         = 1;
    int _orderMu          = 1;
    int _degreeQuadrature = 12;

    int _numSlipBoundaries = 0;

    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e8, 50, 4, 1.};

    // Physical properties
    double _rhoA, _rhoB, _viscA, _viscB;
    double _gamma, _epsilon, _M, _lambda;

    double _gravityVal;

    // The initial epsilon on uniform mesh
    double _epsilon0;
    // The minimum epsilon that can be reached with mesh adaptation
    double _epsilonMin = 2e-3;

    // Generic CHNS functions
    feFunction *_density, *_drhodphi, *_viscosity, *_dviscdphi, *_mobility;

    // Dummy source terms
    feConstantVectorFunction *_uSource;
    feConstantFunction *_pSource, *_phiSource, *_muSource;

    // Functions initialized by the derived class (specific test cases)
    feFunction *_phiInit = nullptr;
    feVectorFunction *_volumeForce = nullptr;

  public:
    SolverCHNS_VolumeAveraged()
    {
      _isCHNS = true; 
    };
    virtual ~SolverCHNS_VolumeAveraged(){};

    int getNumFields() { return 4; }
    std::string getReadDir(){ return _readDir; };
    std::string getWriteDir(){ return _writeDir; };
    std::string getInitialMesh(){ return std::string(_initialMesh); };
    int getOrderForAdaptation() { return _orderPhi; };

    // SolverCHNS_VolumeAveraged
    void printParameters() override
    {
      feInfo("Volume-averaged CHNS Solver parameters :");
      PRINT_STR(_readDir);
      PRINT_STR(_writeDir);
      PRINT_STR(std::string(_initialMesh));

      PRINT_INT(_orderVelocity);
      PRINT_INT(_orderPressure);
      PRINT_INT(_orderPhi);
      PRINT_INT(_orderMu);

      PRINT_DBL(_rhoA);
      PRINT_DBL(_rhoB);
      PRINT_DBL(_viscA);
      PRINT_DBL(_viscB);

      PRINT_DBL(_gamma);
      PRINT_DBL(_epsilon);
      PRINT_DBL(_M);
      PRINT_DBL(_lambda);
    }

    // SolverCHNS_VolumeAveraged
    feStatus parseGenericCHNSParameters(int argc, char **argv)
    {
      const char* writeDirectory = _writeDir.data();

      bool ignoreUnrecognizedOptions = true;
      feOptionsParser options(argc, argv, ignoreUnrecognizedOptions);
      options.addOption(&_initialMesh, "-m", "--m", "");
      options.addOption(&writeDirectory, "-dir", "--dir", "");
      options.addOption(&_orderVelocity, "-ou", "--orderU", "");
      options.addOption(&_orderPressure, "-op", "--orderP", "");
      options.addOption(&_orderPhi, "-ophi", "--orderPhi", "");
      options.addOption(&_orderMu , "-omu", "--orderMu", "");
      options.addOption(&_numSlipBoundaries, "-nSlipBoundaries", "--nSlipBoundaries", "");
      feCheckReturn(options.parse());

      _writeDir = std::string(writeDirectory);
      feInfo("Write directory is : %s", _writeDir.data());

      return FE_STATUS_OK;
    }

    // SolverCHNS_VolumeAveraged
    // Create functions once all parameters have been parsed
    feStatus createCHNSFunctions()
    {
      _lambda = 3. / (2.*sqrt(2.)) * _epsilon * _gamma;

      _density         = new feFunction(      rho_f, {_rhoA, _rhoB});
      _drhodphi        = new feFunction( drhodphi_f, {_rhoA, _rhoB});
      _viscosity       = new feFunction(viscosity_f, {_viscA, _viscB});
      _dviscdphi       = new feFunction(dviscdphi_f, {_viscA, _viscB});
      _mobility        = new feFunction( mobility_f, {_M});

      _uSource         = new feConstantVectorFunction({0., 0.});
      _pSource         = new feConstantFunction(0.);
      _phiSource       = new feConstantFunction(0.);
      _muSource        = new feConstantFunction(0.);

      return FE_STATUS_OK;
    }

    //
    // Update parameters as the mesh is adapted
    //
    // Decrease interface thickness (epsilon)
    //
    feStatus updateParameters(const int fixedPointIteration) override
    {
      UNUSED(fixedPointIteration);
      if(fixedPointIteration > 0 && (fixedPointIteration % 2) == 0) {
        _epsilon = fmax(_epsilonMin, _epsilon/3.);
      }

      // Update parameters depending on epsilon
      _lambda = 3. / (2.*sqrt(2.)) * _epsilon * _gamma;

      // Update functions depending on epsilon or M
      _mobility->resetParameters({_M});

      return FE_STATUS_OK;
    }

    // SolverCHNS_VolumeAveraged
    feStatus createSpaces(feMesh *mesh,
                          std::vector<feSpace*> &spaces,
                          std::vector<feSpace*> &essentialSpaces)
    {
      for(auto *s : spaces) delete s;
      spaces.clear();
      essentialSpaces.clear();

      feSpace *u = nullptr;
      feSpace *uNoSlip = nullptr;
      feSpace *uSlipX = nullptr, *uSlipY = nullptr;
      feSpace *p = nullptr, *pPoint = nullptr;
      feSpace *phi = nullptr, *mu = nullptr;

      feCheckReturn(createFiniteElementSpace(      u, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "Domaine", _degreeQuadrature, &vectorConstant::zero));

      if(_numSlipBoundaries == 4) {
        // Slip on all boundaries
        feCheckReturn(createFiniteElementSpace(uSlipX, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "SlipX", _degreeQuadrature, &vectorConstant::zero));
        feCheckReturn(createFiniteElementSpace(uSlipY, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "SlipY", _degreeQuadrature, &vectorConstant::zero));
      } else if(_numSlipBoundaries == 2)
      {
        // Slip on sides and no-slip on top/bottom
        feCheckReturn(createFiniteElementSpace(uNoSlip, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "NoSlip", _degreeQuadrature, &vectorConstant::zero));
        feCheckReturn(createFiniteElementSpace( uSlipY, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",  "SlipY", _degreeQuadrature, &vectorConstant::zero));
      } else if(_numSlipBoundaries == 0)
      {
        // No slip everywhere
        feCheckReturn(createFiniteElementSpace(uNoSlip, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "NoSlip", _degreeQuadrature, &vectorConstant::zero));
      } else {
        feErrorMsg(FE_STATUS_ERROR, "Num of slip boundaries not handled: %d", _numSlipBoundaries);
        exit(-1);
      }

      feCheckReturn(createFiniteElementSpace(      p, mesh, elementType::LAGRANGE, _orderPressure,   "P",       "Domaine", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace( pPoint, mesh, elementType::LAGRANGE, _orderPressure,   "P", "PressurePoint", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(    phi, mesh, elementType::LAGRANGE, _orderPhi,        "Phi",     "Domaine", _degreeQuadrature, _phiInit));
      feCheckReturn(createFiniteElementSpace(     mu, mesh, elementType::LAGRANGE, _orderMu,         "Mu",      "Domaine", _degreeQuadrature, &scalarConstant::zero));

      if(_numSlipBoundaries == 4) {
        uSlipY->setEssentialComponent(0, true);
        uSlipX->setEssentialComponent(1, true);
        spaces = {u, p, phi, mu, uSlipX, uSlipY, pPoint};
        essentialSpaces = {pPoint};
      }
      else if(_numSlipBoundaries == 2)
      {
        uSlipY->setEssentialComponent(0, true);
        spaces = {u, p, phi, mu, uNoSlip, uSlipY, pPoint};
        essentialSpaces = {uNoSlip, pPoint};
      }
      else if(_numSlipBoundaries == 0)
      {
        spaces = {u, p, phi, mu, uNoSlip, pPoint};
        essentialSpaces = {uNoSlip, pPoint};
      }

      return FE_STATUS_OK;
    }

    // SolverCHNS_VolumeAveraged
    feStatus solve(AdaptData &adapter,
                   const int iInterval,
                   feMesh *mesh,
                   feSolution *sol,
                   feMetaNumber *numbering,
                   std::vector<feSpace*> &spaces,
                   const int /*numUnknowns*/,
                   feExportData &exportData)
    {
      std::vector<feBilinearForm*> forms = {};

      feSpace *u   = spaces[0];
      feSpace *p   = spaces[1];
      feSpace *phi = spaces[2];
      feSpace *mu  = spaces[3];

      std::vector<double> CHNSparameters = {_gamma, _epsilon};

      feBilinearForm *CHNS;
      feCheck(createBilinearForm(CHNS, {u, p, phi, mu}, new CHNS_VolumeAveraged<2>(
        _density,
        _drhodphi,
        _viscosity,
        _dviscdphi,
        _mobility,
        _volumeForce,
        _pSource,
        _uSource,
        _phiSource,
        _muSource,
        CHNSparameters)));
      forms.push_back(CHNS);

      feLinearSystem *system = nullptr;
      #if defined(HAVE_MKL)
        feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
        feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering));
      #else
        feCheckReturn(createLinearSystem(system, PETSC, forms, numbering));
      #endif

      std::vector<feNorm*> norms = {};

      TimeIntegrator *timeIntegrator = nullptr;
      const double t0      = adapter.currentTime;
      const double t1      = adapter.currentTime + adapter.dt * adapter.nTimeStepsPerInterval;
      const int nTimeSteps = adapter.nTimeStepsPerInterval;
      feCheckReturn(createTimeIntegrator(timeIntegrator, timeIntegratorScheme::BDF2, NLoptions,
        system, sol, mesh, norms, exportData, t0, t1, nTimeSteps));

      // Restart from solution container if not the first sub-interval
      // (restart from last solutions of previous sub-interval)
      if(iInterval > 0) {
        BDF2Integrator *ptr = dynamic_cast<BDF2Integrator*> (timeIntegrator);
        
        ptr->_sC->NaNify();
        feCheckReturn(timeIntegrator->restartFromContainer(*adapter.allContainers[iInterval]));
        timeIntegrator->setCurrentStep(adapter.currentStep);

        feInfo("Time           of timeIntegrator : %f", timeIntegrator->getCurrentTime());
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[0]);
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[1]);
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[2]);
        feInfo("Step           of timeIntegrator : %d", timeIntegrator->getCurrentStep());
      }

      // Solve and compute metrics at each time step
      for(int iStep = 0; iStep < adapter.nTimeStepsPerInterval; ++iStep)
      {
        tic();
        feCheckReturn(timeIntegrator->makeSteps(1));
        feInfo("Computed time step in %f s", toc());

        // Reconstruct phi
        std::string meshFile = adapter.allOptions[iInterval].backgroundMeshfile;
        std::string recoveryFile = adapter.writeDir + "phi.msh";
        feNewRecovery recPhi(phi, 0, mesh, sol, meshFile, recoveryFile, false, false, false, nullptr, numbering);

        // Compute metric field and increment metric integral Hi
        feMetric *metricField = adapter.allMetrics[iInterval];
        computeAndIncrementMetric(timeIntegrator->getCurrentTime(), iStep, iInterval, adapter, metricField, &recPhi);

        adapter.currentStep++;
      }

      adapter.currentTime = timeIntegrator->getCurrentTime();

      delete adapter.allContainers[iInterval];
      adapter.allContainers[iInterval] = new feSolutionContainer();
      adapter.allContainers[iInterval]->NaNify();
      *(adapter.allContainers[iInterval]) = timeIntegrator->getSolutionContainer();

      delete timeIntegrator;
      delete system;
      for(feBilinearForm *f : forms)
        delete f;

      return FE_STATUS_OK;
    }

    // SolverCHNS_VolumeAveraged
    feStatus computeLastMetric(AdaptData &adapter,
                               const int iInterval,
                               feMesh *mesh,
                               feSolution *sol,
                               feMetaNumber *numbering,
                               std::vector<feSpace*> &spaces)
    {
      feSpace *phi = spaces[2];

      // Reconstruct velocity
      std::string meshFile = adapter.allOptions[iInterval].backgroundMeshfile;
      std::string recoveryFile = adapter.writeDir + "phi.msh";
      feNewRecovery rec(phi, 0, mesh, sol, meshFile, recoveryFile, false, false, false, nullptr, numbering);

      // Compute "steady state" metric field
      feMetric *metricField = adapter.allMetrics[iInterval];
      metricField->setRecoveredFields({&rec});
      metricField->_options.debug = true;
      metricField->setMetricScaling(true);
      feCheckReturn(metricField->computeMetrics());

      return FE_STATUS_OK;
    }

    //
    // Compute error on scalar fields
    //
    feStatus computeError(const int iIntervalReference,
                          AdaptData &adapter,
                          const std::vector<feSpace*> &spaces,
                          feMesh *mesh,
                          feSolution *sol,
                          std::vector<double> &errors)
    {
      feInfo("Computing error for CHNS based on reference finer solution");

      errors.resize(this->getNumFields(), 0.);

      // feSpace *u   = spaces[0];
      feSpace *p   = spaces[1];
      feSpace *phi = spaces[2];
      feSpace *mu  = spaces[3];

      feNorm *errorOnScalarSpaces = nullptr;
      feCheckReturn(createNorm(errorOnScalarSpaces, NONE, {p, phi, mu}, sol));
      feCheckReturn(errorOnScalarSpaces->computeErrorNormFromExternalSolution(adapter.referenceTestCase[iIntervalReference].mesh,
                                                                              adapter.referenceTestCase[iIntervalReference].sol,
                                                                              adapter.referenceTestCase[iIntervalReference].spaces,
                                                                              mesh, sol, errors));
      feInfo("Computed error on P  : %+-1.6e", errors[0]);
      feInfo("Computed error on Phi: %+-1.6e", errors[1]);
      feInfo("Computed error on Mu : %+-1.6e", errors[2]);

      delete errorOnScalarSpaces;
      return FE_STATUS_OK;
    }
  }; // SolverCHNS_VolumeAveraged

  class GenericSolverCHNS_Alternative : public GenericSolver
  {
  protected:
    std::string _readDir = "";
    std::string _writeDir = "";
    const char *_initialMesh = "";

    int _orderVelocity    = 2;
    int _orderPressure    = 1;
    int _orderPhi         = 1;
    int _orderMu          = 1;
    int _degreeQuadrature = 12;

    int _numSlipBoundaries = 0;

    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 50, 4, 1.};

    // Physical properties
    double _rhoA, _rhoB, _viscA, _viscB;
    double _gamma, _epsilon, _M, _alpha, _beta, _tau;

    double _gravityVal;

    // The initial epsilon on uniform mesh
    double _epsilon0;
    // The minimum epsilon that can be reached with mesh adaptation
    double _epsilonMin = 2e-3;

    // Functions initialized by the derived class (specific test cases)
    feFunction *_phiInit = nullptr;
    feVectorFunction *_volumeForce = nullptr;

    // The following functions are specific to this CHNS model
    // and should be initialized in this base class :

    // Generic CHNS functions
    feFunction *_density, *_drhodphi, *_viscosity, *_dviscdphi, *_mobility;

    // Coefficients of the continuity equation
    feFunction *_coeffDiffMu_p, *_coeffDiffP_p;
    
    // Coefficients of the tracer equation
    feFunction *_coeffDiffMu_phi, *_coeffDiffP_phi;

    // Coefficients of the chemical potential equation
    feConstantFunction *_coeffDiffPhi, *_coeffPhiCube, *_coeffPhi;

    // Dummy source terms
    feConstantVectorFunction *_uSource;
    feConstantFunction *_pSource, *_phiSource, *_muSource;

  public:
    GenericSolverCHNS_Alternative()
    {
      _isCHNS = true; 
    };
    virtual ~GenericSolverCHNS_Alternative(){};

    int getNumFields() { return 4; }
    std::string getReadDir(){ return _readDir; };
    std::string getWriteDir(){ return _writeDir; };
    std::string getInitialMesh(){ return std::string(_initialMesh); };
    int getOrderForAdaptation() { return _orderPhi; };

    // GenericSolverCHNS_Alternative
    void printParameters() override
    {
      feInfo("Alternative CHNS Solver parameters :");
      PRINT_STR(_readDir);
      PRINT_STR(_writeDir);
      PRINT_STR(std::string(_initialMesh));

      PRINT_INT(_orderVelocity);
      PRINT_INT(_orderPressure);
      PRINT_INT(_orderPhi);
      PRINT_INT(_orderMu);

      PRINT_DBL(_rhoA);
      PRINT_DBL(_rhoB);
      PRINT_DBL(_viscA);
      PRINT_DBL(_viscB);

      PRINT_DBL(_gamma);
      PRINT_DBL(_epsilon);
      PRINT_DBL(_M);
      PRINT_DBL(_alpha);
      PRINT_DBL(_beta);
      PRINT_DBL(_tau);
    }

    // GenericSolverCHNS_Alternative
    feStatus parseGenericCHNSParameters(int argc, char **argv)
    {
      const char* writeDirectory = _writeDir.data();

      bool ignoreUnrecognizedOptions = true;
      feOptionsParser options(argc, argv, ignoreUnrecognizedOptions);
      options.addOption(&_initialMesh, "-m", "--m", "");
      options.addOption(&writeDirectory, "-dir", "--dir", "");
      options.addOption(&_orderVelocity, "-ou", "--orderU", "");
      options.addOption(&_orderPressure, "-op", "--orderP", "");
      options.addOption(&_orderPhi, "-ophi", "--orderPhi", "");
      options.addOption(&_orderMu , "-omu", "--orderMu", "");
      // options.addOption(&_allSlip, "-allSlip", "--allSlip", "-no-allSlip", "--no-allSlip", "");
      options.addOption(&_numSlipBoundaries, "-nSlipBoundaries", "--nSlipBoundaries", "");
      feCheckReturn(options.parse());

      _writeDir = std::string(writeDirectory);

      return FE_STATUS_OK;
    }

    // GenericSolverCHNS_Alternative
    // Create functions once all parameters have been parsed
    feStatus createCHNSFunctions()
    {
      _alpha = (_rhoB - _rhoA) / (_rhoA + _rhoB);
      _tau   = 3. * _gamma / (2. * sqrt(2.)) * _epsilon;
      _beta  = 3. * _gamma / (2. * sqrt(2.)) / _epsilon;

      _density         = new feFunction(      rho_f, {_rhoA, _rhoB});
      _drhodphi        = new feFunction( drhodphi_f, {_rhoA, _rhoB});
      _viscosity       = new feFunction(viscosity_f, {_viscA, _viscB});
      _dviscdphi       = new feFunction(dviscdphi_f, {_viscA, _viscB});
      _mobility        = new feFunction( mobility_f, {_M});
      // _mobility = nullptr;

      _coeffDiffMu_p   = new feFunction(mobilityAlpha_f, {_M, _alpha});
      _coeffDiffP_p    = new feFunction(mobilityAlphaSquared_f, {_M, _alpha}); 

      _coeffDiffMu_phi = new feFunction(mobility_f, {_M, _alpha});
      _coeffDiffP_phi  = new feFunction(mobilityAlpha_f, {_M, _alpha});

      _coeffDiffPhi    = new feConstantFunction(- _tau);
      _coeffPhiCube    = new feConstantFunction(- _beta);
      _coeffPhi        = new feConstantFunction(+ _beta);

      _uSource         = new feConstantVectorFunction({0., 0.});
      _pSource         = new feConstantFunction(0.);
      _phiSource       = new feConstantFunction(0.);
      _muSource        = new feConstantFunction(0.);

      return FE_STATUS_OK;
    }

    //
    // Update parameters as the mesh is adapted
    //
    // Decrease interface thickness (epsilon)
    //
    feStatus updateParameters(const int fixedPointIteration) override
    {
      if(fixedPointIteration > 0) {
        _epsilon = fmax(_epsilonMin, _epsilon/2.);
      }

      // Update parameters depending on epsilon
      _M = 0.1 * _epsilon * _epsilon;
      _tau   = 3. * _gamma / (2. * sqrt(2.)) * _epsilon;
      _beta  = 3. * _gamma / (2. * sqrt(2.)) / _epsilon;

      // Update functions depending on epsilon, M, tau and/or beta
           _mobility->resetParameters({_M});

      _coeffDiffMu_p->resetParameters({_M, _alpha});
       _coeffDiffP_p->resetParameters({_M, _alpha});

      _coeffDiffMu_phi->resetParameters({_M, _alpha});
       _coeffDiffP_phi->resetParameters({_M, _alpha});

      // Coefficients for mu equation
      _coeffDiffPhi->resetParameters({- _tau});
      _coeffPhiCube->resetParameters({- _beta});
          _coeffPhi->resetParameters({+ _beta});

      return FE_STATUS_OK;
    }

    // GenericSolverCHNS_Alternative
    feStatus createSpaces(feMesh *mesh,
                          std::vector<feSpace*> &spaces,
                          std::vector<feSpace*> &essentialSpaces)
    {
      for(auto *s : spaces) delete s;
      spaces.clear();
      essentialSpaces.clear();

      feSpace *u = nullptr;
      feSpace *uNoSlip = nullptr;
      feSpace *uSlipX = nullptr, *uSlipY = nullptr;
      feSpace *p = nullptr, *pPoint = nullptr;
      feSpace *phi = nullptr, *mu = nullptr;

      feCheckReturn(createFiniteElementSpace(      u, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "Domaine", _degreeQuadrature, &vectorConstant::zero));

      if(_numSlipBoundaries == 4) {
        // Slip on all boundaries
        feCheckReturn(createFiniteElementSpace(uSlipX, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "SlipX", _degreeQuadrature, &vectorConstant::zero));
        feCheckReturn(createFiniteElementSpace(uSlipY, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "SlipY", _degreeQuadrature, &vectorConstant::zero));
      } else if(_numSlipBoundaries == 2)
      {
        // Slip on sides and no-slip on top/bottom
        feCheckReturn(createFiniteElementSpace(uNoSlip, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "NoSlip", _degreeQuadrature, &vectorConstant::zero));
        feCheckReturn(createFiniteElementSpace( uSlipY, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",  "SlipY", _degreeQuadrature, &vectorConstant::zero));
      } else if(_numSlipBoundaries == 0)
      {
        // No slip everywhere
        feCheckReturn(createFiniteElementSpace(uNoSlip, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U", "NoSlip", _degreeQuadrature, &vectorConstant::zero));
      } else {
        feErrorMsg(FE_STATUS_ERROR, "Num of slip boundaries not handled: %d", _numSlipBoundaries);
        exit(-1);
      }

      feCheckReturn(createFiniteElementSpace(      p, mesh, elementType::LAGRANGE, _orderPressure,   "P",       "Domaine", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace( pPoint, mesh, elementType::LAGRANGE, _orderPressure,   "P", "PressurePoint", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(    phi, mesh, elementType::LAGRANGE, _orderPhi,        "Phi",     "Domaine", _degreeQuadrature, _phiInit));
      feCheckReturn(createFiniteElementSpace(     mu, mesh, elementType::LAGRANGE, _orderMu,         "Mu",      "Domaine", _degreeQuadrature, &scalarConstant::zero));

      if(_numSlipBoundaries == 4) {
        uSlipY->setEssentialComponent(0, true);
        uSlipX->setEssentialComponent(1, true);
        spaces = {u, p, phi, mu, uSlipX, uSlipY, pPoint};
        essentialSpaces = {pPoint};
      }
      else if(_numSlipBoundaries == 2)
      {
        uSlipY->setEssentialComponent(0, true);
        spaces = {u, p, phi, mu, uNoSlip, uSlipY, pPoint};
        essentialSpaces = {uNoSlip, pPoint};
      }
      else if(_numSlipBoundaries == 0)
      {
        spaces = {u, p, phi, mu, uNoSlip, pPoint};
        essentialSpaces = {uNoSlip, pPoint};
      }

      return FE_STATUS_OK;
    }

    // GenericSolverCHNS_Alternative
    feStatus solve(AdaptData &adapter,
                   const int iInterval,
                   feMesh *mesh,
                   feSolution *sol,
                   feMetaNumber *numbering,
                   std::vector<feSpace*> &spaces,
                   const int /*numUnknowns*/,
                   feExportData &exportData)
    {
      std::vector<feBilinearForm*> forms = {};

      feSpace *u   = spaces[0];
      feSpace *p   = spaces[1];
      feSpace *phi = spaces[2];
      feSpace *mu  = spaces[3];

      /////////////////////////////////////////////////////////////
      // // Continuity
      // feBilinearForm *divU = nullptr, *diffMu_p = nullptr, *diffp_p = nullptr;
      // feCheck(createBilinearForm(    divU,  {p, u}, new feSysElm_MixedDivergence<2>(&scalarConstant::one)));
      // // feCheck(createBilinearForm(diffMu_p, {p, mu}, new feSysElm_MixedGradGrad<2>(_coeffDiffMu_p)));
      // // feCheck(createBilinearForm( diffp_p,     {p}, new feSysElm_Diffusion<2>(_coeffDiffP_p)));
      // feCheck(createBilinearForm(diffMu_p, {p, mu, phi}, new feSysElm_MixedGradGradFieldDependentCoeff<2>(_coeffDiffMu_p)));
      // feCheck(createBilinearForm( diffp_p,     {p, phi}, new feSysElm_DiffusionFieldDependentCoeff<2>(_coeffDiffP_p)));
      // forms.push_back(divU);
      // forms.push_back(diffMu_p);
      // forms.push_back(diffp_p);

      // // Momentum
      // feBilinearForm *CHNS_momentum = nullptr;
      // feCheckReturn(createBilinearForm(CHNS_momentum, {u, p, phi, mu}, new feSysElm_CHNS_Momentum_Alternative<2>(_density,
      //                                                                                                            _drhodphi,
      //                                                                                                            _viscosity,
      //                                                                                                            _dviscdphi,
      //                                                                                                            _volumeForce)));
      // forms.push_back(CHNS_momentum);

      // // Tracer convection
      // feBilinearForm *dphidt = nullptr, *divPhiU = nullptr, *diffMu_phi = nullptr, *diffP_phi = nullptr;
      // feCheck(createBilinearForm(    dphidt,     {phi}, new feSysElm_TransientMass(&scalarConstant::one)));
      // feCheck(createBilinearForm(   divPhiU,  {phi, u}, new feSysElm_ScalarVectorProduct<2>(&scalarConstant::minusOne)));
      // // feCheck(createBilinearForm(diffMu_phi, {phi, mu}, new feSysElm_MixedGradGrad<2>(_coeffDiffMu_phi)));
      // // feCheck(createBilinearForm( diffP_phi,  {phi, p}, new feSysElm_MixedGradGrad<2>(_coeffDiffP_phi)));
      // feCheck(createBilinearForm(diffMu_phi, {phi, mu, phi}, new feSysElm_MixedGradGradFieldDependentCoeff<2>(_coeffDiffMu_phi)));
      // feCheck(createBilinearForm( diffP_phi,  {phi, p, phi}, new feSysElm_MixedGradGradFieldDependentCoeff<2>(_coeffDiffP_phi)));
      // forms.push_back(dphidt);
      // forms.push_back(divPhiU);
      // forms.push_back(diffMu_phi);
      // forms.push_back(diffP_phi);

      // // Chemical potential
      // feBilinearForm *muMass = nullptr, *phiCube = nullptr, *phiMass = nullptr, *diffPhi = nullptr;
      // feCheck(createBilinearForm( muMass,      {mu}, new feSysElm_Mass(&scalarConstant::one)));
      // feCheck(createBilinearForm(diffPhi, {mu, phi}, new feSysElm_MixedGradGrad<2>(_coeffDiffPhi)));
      // feCheck(createBilinearForm(phiCube, {mu, phi}, new feSysElm_MixedMassPower(_coeffPhiCube, 3)));
      // feCheck(createBilinearForm(phiMass, {mu, phi}, new feSysElm_MixedMass(_coeffPhi)));
      // forms.push_back(muMass);
      // forms.push_back(diffPhi);
      // forms.push_back(phiCube);
      // forms.push_back(phiMass);
      /////////////////////////////////////////////////////////////
      // Or with a single weak form:
      std::vector<double> CHNSparameters = {_alpha, _gamma, _epsilon};

      feBilinearForm *CHNS;
      feCheck(createBilinearForm(CHNS, {u, p, phi, mu}, new feSysElm_CHNS_Alternative<2>(
        _density,
        _drhodphi,
        _viscosity,
        _dviscdphi,
        _mobility,
        _volumeForce,
        _pSource,
        _uSource,
        _phiSource,
        _muSource,
        CHNSparameters)));
      forms.push_back(CHNS);
      /////////////////////////////////////////////////////////////

      // for(auto &f : forms) {
      //   f->setComputeMatrixWithFD(true);
      // }

      feLinearSystem *system = nullptr;
      #if defined(HAVE_MKL)
        feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
        feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering));
      #else
        feCheckReturn(createLinearSystem(system, PETSC, forms, numbering));
      #endif

      std::vector<feNorm*> norms = {};

      TimeIntegrator *timeIntegrator = nullptr;
      const double t0 = adapter.currentTime;
      const double t1 = adapter.currentTime + adapter.dt * adapter.nTimeStepsPerInterval;
      const int nTimeSteps = adapter.nTimeStepsPerInterval;
      feCheckReturn(createTimeIntegrator(timeIntegrator, timeIntegratorScheme::BDF2, NLoptions,
        system, sol, mesh, norms, exportData, t0, t1, nTimeSteps));

      // Restart from solution container if not the first sub-interval
      // (restart from last solutions of previous sub-interval)
      if(iInterval > 0) {
        // BDF1Integrator *ptr = dynamic_cast<BDF1Integrator*> (timeIntegrator);
        BDF2Integrator *ptr = dynamic_cast<BDF2Integrator*> (timeIntegrator);
        
        ptr->_sC->NaNify();
        feCheckReturn(timeIntegrator->restartFromContainer(*adapter.allContainers[iInterval]));
        timeIntegrator->setCurrentStep(adapter.currentStep);

        feInfo("Time           of timeIntegrator : %f", timeIntegrator->getCurrentTime());
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[0]);
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[1]);
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[2]);
        feInfo("Step           of timeIntegrator : %d", timeIntegrator->getCurrentStep());
      }

      // Solve and compute metrics at each time step
      for(int iStep = 0; iStep < adapter.nTimeStepsPerInterval; ++iStep)
      {
        tic();
        feCheckReturn(timeIntegrator->makeSteps(1));
        feInfo("Computed time step in %f s", toc());

        // Reconstruct phi
        std::string meshFile = adapter.allOptions[iInterval].backgroundMeshfile;
        std::string recoveryFile = adapter.writeDir + "phi.msh";
        feNewRecovery recPhi(phi, 0, mesh, sol, meshFile, recoveryFile, false, false, false, nullptr, numbering);

        // Compute metric field and increment metric integral Hi
        feMetric *metricField = adapter.allMetrics[iInterval];
        computeAndIncrementMetric(timeIntegrator->getCurrentTime(), iStep, iInterval, adapter, metricField, &recPhi);

        adapter.currentStep++;
      }

      adapter.currentTime = timeIntegrator->getCurrentTime();

      delete adapter.allContainers[iInterval];
      adapter.allContainers[iInterval] = new feSolutionContainer();
      adapter.allContainers[iInterval]->NaNify();
      *(adapter.allContainers[iInterval]) = timeIntegrator->getSolutionContainer();

      delete timeIntegrator;
      delete system;
      for(feBilinearForm *f : forms)
        delete f;

      return FE_STATUS_OK;
    }

    // GenericSolverCHNS_Alternative
    feStatus computeLastMetric(AdaptData &adapter,
                               const int iInterval,
                               feMesh *mesh,
                               feSolution *sol,
                               feMetaNumber *numbering,
                               std::vector<feSpace*> &spaces)
    {
      feSpace *phi = spaces[2];

      // Reconstruct velocity
      std::string meshFile = adapter.allOptions[iInterval].backgroundMeshfile;
      std::string recoveryFile = adapter.writeDir + "phi.msh";
      feNewRecovery rec(phi, 0, mesh, sol, meshFile, recoveryFile, false, false, false, nullptr, numbering);

      // Compute "steady state" metric field
      feMetric *metricField = adapter.allMetrics[iInterval];
      metricField->setRecoveredFields({&rec});
      metricField->_options.debug = true;
      metricField->setMetricScaling(true);
      feCheckReturn(metricField->computeMetrics());

      return FE_STATUS_OK;
    }

    //
    // Compute error on scalar fields
    //
    feStatus computeError(const int iIntervalReference,
                          AdaptData &adapter,
                          const std::vector<feSpace*> &spaces,
                          feMesh *mesh,
                          feSolution *sol,
                          std::vector<double> &errors)
    {
      feInfo("Computing error for CHNS based on reference finer solution");

      errors.resize(this->getNumFields(), 0.);

      // feSpace *u   = spaces[0];
      feSpace *p   = spaces[1];
      feSpace *phi = spaces[2];
      feSpace *mu  = spaces[3];

      feNorm *errorOnScalarSpaces = nullptr;
      feCheckReturn(createNorm(errorOnScalarSpaces, NONE, {p, phi, mu}, sol));
      feCheckReturn(errorOnScalarSpaces->computeErrorNormFromExternalSolution(adapter.referenceTestCase[iIntervalReference].mesh,
                                                                              adapter.referenceTestCase[iIntervalReference].sol,
                                                                              adapter.referenceTestCase[iIntervalReference].spaces,
                                                                              mesh, sol, errors));
      feInfo("Computed error on P  : %+-1.6e", errors[0]);
      feInfo("Computed error on Phi: %+-1.6e", errors[1]);
      feInfo("Computed error on Mu : %+-1.6e", errors[2]);

      delete errorOnScalarSpaces;
      return FE_STATUS_OK;
    }
  }; // GenericSolverCHNS_Alternative
} // namespace CHNS_Generic

namespace CHNS_RT
{
  double phiInit_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    double x = args.pos[0];
    double y = args.pos[1];
    double h = par[0];
    double L = par[1];
    double epsilon = par[2];
    double y0 = h/2. + 0.1 * L * cos(2.*M_PI*x/L);
    return tanh((y - y0)/(sqrt(2.)*epsilon));
  }

  class Solver : public CHNS_Generic::GenericSolverCHNS_Alternative
  {
  public:
    // Geometry
    double h = 1.;
    double L = h/4.;

  public:
    Solver() : CHNS_Generic::GenericSolverCHNS_Alternative()
    {
      _readDir = "";
      _writeDir = "../data/CHNS/RayleighTaylor/";
      _initialMesh = "../data/CHNS/RayleighTaylor/mediumAllSlip.msh";

      // // Symmetric test case
      // _initialMesh = "../data/CHNS/RT_symmetric/mesh.msh";
      // _writeDir = "../data/CHNS/RT_symmetric/";
    };

    // Parse for CHNS parameters
    feStatus parseAndCreateFunctions(int argc, char **argv)
    {
      feCheckReturn(parseGenericCHNSParameters(argc, argv));

      // Parse specific test case parameters here, if needed
      //
      // ...

      _rhoA          = 300.;
      _rhoB          = 100.;
      _viscA         = 1.53e-1;
      _viscB         = 1.53e-1;
      _gamma         = 1e-2;
      _gravityVal    = -9.81;

      // _epsilon       = 2.766e-3;
      // _mobilityVal   = 7.67e-6;

      // const double mesh_size = 1./32;
      const double mesh_size = 1./75; // mediumAllSlip.msh
      _epsilon0 = 0.64 * mesh_size;
      _epsilon = _epsilon0;
      _M = 0.1 * _epsilon * _epsilon;

      _phiInit       = new feFunction(phiInit_f, {h, L, _epsilon});
      _volumeForce   = new feConstantVectorFunction({0, _gravityVal});
      createCHNSFunctions();

      return FE_STATUS_OK;
    }

    //
    // Reconstruct isoline phi = 0 (interface)
    //
    feStatus computePostProcessing(AdaptData &/*adapter*/,
                                   const std::vector<feSpace*> &spaces,
                                   feMesh */*mesh*/,
                                   feSolution *sol,
                                   PostProcData &data) override
    {
      const feSpace *phi = spaces[2];

      feNorm *interface;
      feCheckReturn(createNorm(interface, NONE, {phi}, sol));
      const double targetValue = 0.;
      feCheckReturn(interface->reconstructScalarIsoline(targetValue, data.isoline));
      delete interface;

      data.res = {};

      return FE_STATUS_OK;
    }
  };
} // namespace CHNS_RT

namespace CHNS_Sloshing
{
  double phiInit_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    double y = args.pos[1];
    double epsilon = par[0];
    double h0 = par[1];
    return tanh((y - h0)/(sqrt(2.)*epsilon));
  }

  void gravity_f(const feFunctionArguments &args, const std::vector<double> &/*par*/, std::vector<double> &res)
  {
    double t = args.t;
    // double x = args.pos[0];
    // double y = args.pos[1];
    res[0] = sin(5*M_PI*t);
    res[1] = -9.81;
  }

  class Solver : public CHNS_Generic::GenericSolverCHNS
  {
  public:
    // Geometry
    double H  = 2.; // Height of the tank
    double h0 = 1.; // Initial water height
    double L  = 1.; // Width of the tank

  public:
    Solver() : CHNS_Generic::GenericSolverCHNS()
    {
      _readDir = "";
      _writeDir = "../data/CHNS/sloshing/";
      _initialMesh = "../data/CHNS/sloshing/mesh.msh";

      _rhoA          = 1.;
      _rhoB          = 100.;
      _viscA         = 1.53e-1;
      _viscB         = 1.53e-1;
      _gamma         = 1e-2;
      _epsilon       = 2.766e-3;
      _mobilityVal   = 7.67e-6;
    };

    // Parse for CHNS parameters
    feStatus parseAndCreateFunctions(int argc, char **argv)
    {
      feCheckReturn(parseGenericCHNSParameters(argc, argv));
      phiInit       = new feFunction(phiInit_f, {_epsilon, h0});
      volumeForce   = new feVectorFunction(gravity_f);
      createCHNSFunctions();

      return FE_STATUS_OK;
    }
  };
} // namespace CHNS_Sloshing

namespace CHNS_Bubble
{
  double phiInit_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    double x = args.pos[0];
    double y = args.pos[1];
    double x0 = par[0];
    double y0 = par[1];
    double R  = par[2];
    double epsilon = par[3];
    double d = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
    return tanh((d - R)/(sqrt(2.)*epsilon));
  }

  //
  // Indicator function used to determine whether we are in the bubble or not.
  // Bubble is made of fluid 2, so checks if phi <= 0.
  //
  double indicator_f(const feFunctionArguments &args, const std::vector<double> &)
  {
    const double phi = args.u;
    return (double) phi <= 0.;
  }

  // Same but returns phi if in bubble
  double indicatorPhi_f(const feFunctionArguments &args, const std::vector<double> &)
  {
    const double phi = args.u;
    return phi <= 0. ? phi : 0.;
  }

  // Same but returns x coordinate
  double indicatorX_f(const feFunctionArguments &args, const std::vector<double> &)
  {
    const double phi = args.u;
    const double x = args.pos[0];
    return phi <= 0. ? x : 0.;
  }

  // Same but returns y coordinate
  double indicatorY_f(const feFunctionArguments &args, const std::vector<double> &)
  {
    const double phi = args.u;
    const double y = args.pos[1];
    return phi <= 0. ? y : 0.;
  }

  // Same but returns velocity
  void indicatorVelocity_f(const feFunctionArguments &args,
                           const std::vector<double> &,
                           std::vector<double> &res)
  {
    const double phi = args.u;
    if(phi <= 0.) {
      res[0] = args.v[0];
      res[1] = args.v[1];
      res[2] = args.v[2];
    } else {
      res[0] = res[1] = res[2] = 0.;
    }
  }

  class Solver  : public CHNS_Generic::GenericSolverCHNS_Alternative
  {
  public:
    // Geometry
    double h  = 2.;
    double L  = 1.;
    double x0 = 0.5;
    double y0 = 0.5;
    double R  = 0.25;

    int _testCase = 1;

    // Postprocessing functions
    feFunction *_indicatorArea, *_indicatorX, *_indicatorY, *_indicatorPhi;
    feVectorFunction *_indicatorVelocity;

  public:
    Solver() : CHNS_Generic::GenericSolverCHNS_Alternative()
    {
      _readDir = "../data/CHNS/bulleMontante/";
      _writeDir = "../data/CHNS/bulleMontante/";
      _initialMesh = "../data/CHNS/bulleMontante/coarseSlipNoSlip.msh";
      // _initialMesh = "../data/CHNS/bulleMontante/mediumSlipNoSlip.msh";
    };

    // Parse for CHNS parameters
    feStatus parseAndCreateFunctions(int argc, char **argv)
    {
      feCheckReturn(parseGenericCHNSParameters(argc, argv));

      // Parse specific test case parameters
      feOptionsParser options(argc, argv, true);
      options.addOption(&_testCase, "-testCase", "--testCase", "");
      feCheckReturn(options.parse());

      // First test case
      if(_testCase == 1) {
        _rhoA          = 1000.;
        _rhoB          = 100.;
        _viscA         = 10.;
        _viscB         = 1.;
        _gamma         = 24.5;
      } else {
        // Second test case
        _rhoA          = 1000.;
        _rhoB          = 1.;
        _viscA         = 10.;
        _viscB         = 0.1;
        _gamma         = 1.96;
      }

      _gravityVal    = -0.98;

      const double mesh_size = 1./32;
      _epsilon0 = 0.64 * mesh_size;
      _epsilon = _epsilon0;
      _M = 0.1 * _epsilon * _epsilon;

      createCHNSFunctions();

      _phiInit            = new feFunction(phiInit_f, {x0, y0, R, _epsilon});
      _volumeForce        = new feConstantVectorFunction({0., _gravityVal});
      _indicatorArea     = new feFunction(indicator_f);
      _indicatorPhi      = new feFunction(indicatorPhi_f);
      _indicatorX        = new feFunction(indicatorX_f);
      _indicatorY        = new feFunction(indicatorY_f);
      _indicatorVelocity = new feVectorFunction(indicatorVelocity_f);

      return FE_STATUS_OK;
    }

    // Update specific test case parameters as the mesh is adapted
    feStatus updateParameters(const int fixedPointIteration) override
    {
      GenericSolverCHNS_Alternative::updateParameters(fixedPointIteration);

      // Update derived function
      _phiInit->resetParameters({x0, y0, R, _epsilon});

      return FE_STATUS_OK;
    }

    //
    // Compute values from the bubble benchmark:
    //  - center of mass
    //  - rising velocity
    //  - circularity
    //
    feStatus computePostProcessing(AdaptData &/*adapter*/,
                                   const std::vector<feSpace*> &spaces,
                                   feMesh */*mesh*/,
                                   feSolution *sol,
                                   PostProcData &data) override
    {
      const double time = sol->getCurrentTime();
      const feSpace *u   = spaces[0];
      const feSpace *phi = spaces[2];

      // Compute center of mass
      feNorm *bubbleArea, *intX, *intY;
      feCheckReturn(createNorm(bubbleArea, INTEGRAL_SOL_DEPENDENT_F, {phi}, sol, _indicatorArea));
      feCheckReturn(createNorm(      intX, INTEGRAL_SOL_DEPENDENT_F, {phi}, sol, _indicatorX));
      feCheckReturn(createNorm(      intY, INTEGRAL_SOL_DEPENDENT_F, {phi}, sol, _indicatorY));

      const double area = bubbleArea->compute();
      const double intXval = intX->compute();
      const double intYval = intY->compute();

      const double xc = intXval / area;
      const double yc = intYval / area;

      // Compute rising velocity
      feNorm *intVelocity;
      feCheckReturn(createNorm(intVelocity, VECTOR_INTEGRAL_SOL_DEPENDENT_F, {phi, u}, sol, nullptr, _indicatorVelocity));
      
      std::vector<double> intVelocityVal(3, 0.);
      intVelocity->compute(intVelocityVal);

      const double avgRisingVelocity = intVelocityVal[1] / area;

      // Integral of phi
      feNorm *intPhi;
      feCheck(createNorm(intPhi, INTEGRAL, {phi}, sol));
      const double integral_phi = intPhi->compute();
      delete intPhi;

      // Integral of phi in bubble
      feNorm *intPhiBubble;
      feCheckReturn(createNorm(intPhiBubble, INTEGRAL_SOL_DEPENDENT_F, {phi}, sol, _indicatorPhi));
      const double integral_phi_bubble = intPhiBubble->compute();
      delete intPhiBubble;

      // Reconstruct the interface (isoline phi = 0)
      feNorm *interface;
      feCheckReturn(createNorm(interface, NONE, {phi}, sol));
      const double targetValue = 0.;
      feCheckReturn(interface->reconstructScalarIsoline(targetValue, data.isoline));
      delete interface;

      // Circularity
      // Compute length of each segment of the isoline (assumes a single curve)
      // Isoline has an even number of points
      double perimeter = 0.;
      for(size_t iPt = 0; iPt < data.isoline.size() / 6; ++iPt) {
        double edge[3] = {data.isoline[6 * iPt + 3] - data.isoline[6 * iPt + 0],
                          data.isoline[6 * iPt + 4] - data.isoline[6 * iPt + 1],
                          data.isoline[6 * iPt + 5] - data.isoline[6 * iPt + 2]};
        perimeter += sqrt(edge[0] * edge[0] + edge[1] * edge[1]);
      }
      const double circularity = 2.*M_PI*R / perimeter;

      data.res = {time, area, xc, yc, avgRisingVelocity, integral_phi, integral_phi_bubble, circularity};

      return FE_STATUS_OK;
    }
  };
} // namespace CHNS_Bubble

namespace CHNS_YoungLaplace
{
  double phiInit_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    double x = args.pos[0];
    double y = args.pos[1];
    double x0 = par[0];
    double y0 = par[1];
    double R  = par[2];
    double epsilon = par[3];
    double d = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
    return tanh((d - R)/(sqrt(2.)*epsilon));
  }

  //
  // Indicator function used to determine whether we are in the bubble or not.
  // Bubble is made of fluid 2, so checks if phi <= 0.
  //
  double indicator_f(const feFunctionArguments &args, const std::vector<double> &)
  {
    const double phi = args.u;
    return (double) phi <= 0.;
  }

  // Same but returns p if in bubble
  double indicatorPressureIn_f(const feFunctionArguments &args, const std::vector<double> &)
  {
    const double phi = args.u;
    const double p   = args.userScalarFields[1];
    return phi <= 0. ? p : 0.;
  }

  // Same but returns p if outside bubble
  double indicatorPressureOut_f(const feFunctionArguments &args, const std::vector<double> &)
  {
    const double phi = args.u;
    const double p   = args.userScalarFields[1];
    return phi >= 0. ? p : 0.;
  }

  // class Solver  : public CHNS_Generic::GenericSolverCHNS_Alternative
  class VolumeAveragedSolver  : public CHNS_Generic::SolverCHNS_VolumeAveraged
  {
  public:
    // Geometry
    double x0 = 0.;
    double y0 = 0.;
    double R  = 0.5;

    // Postprocessing functions
    feFunction *_indicatorArea, *_indicatorPin, *_indicatorPout;

  public:
    // Solver() : CHNS_Generic::GenericSolverCHNS_Alternative()
    VolumeAveragedSolver() : CHNS_Generic::SolverCHNS_VolumeAveraged()
    {
      _readDir     = "../data/CHNS/YoungLaplace/";
      _writeDir    = "../data/CHNS/YoungLaplace/";
      // _initialMesh = "../data/CHNS/YoungLaplace/coarse.msh";
      _initialMesh = "../data/CHNS/YoungLaplace/mesh.msh";
    };

    // Parse for CHNS parameters
    feStatus parseAndCreateFunctions(int argc, char **argv)
    {
      feCheckReturn(parseGenericCHNSParameters(argc, argv));

      // Parse specific test case parameters if needed
      feOptionsParser options(argc, argv, true);
      options.addOption(&R, "-radius", "--radius", "");
      feCheckReturn(options.parse());

      _rhoA          = 10.;
      _rhoB          = 10.;
      _viscA         = 1.;
      _viscB         = 1.;
      _gamma         = 1.;

      _gravityVal    = 0.;

      // const double mesh_size = 1./32;
      // _epsilon0 = 0.64 * mesh_size;
      // _epsilon = _epsilon0;
      // _M = 0.1 * _epsilon * _epsilon;

      _epsilon = 1.38e-2;
      _M       = 1.53e-5;

      createCHNSFunctions();

      _phiInit           = new feFunction(phiInit_f, {x0, y0, R, _epsilon});
      _volumeForce       = new feConstantVectorFunction({0., _gravityVal});
      _indicatorArea     = new feFunction(indicator_f);
      _indicatorPin      = new feFunction(indicatorPressureIn_f);
      _indicatorPout     = new feFunction(indicatorPressureOut_f);

      return FE_STATUS_OK;
    }

    // Update specific test case parameters as the mesh is adapted
    feStatus updateParameters(const int fixedPointIteration) override
    {
      // GenericSolverCHNS_Alternative::updateParameters(fixedPointIteration);
      CHNS_Generic::SolverCHNS_VolumeAveraged::updateParameters(fixedPointIteration);

      // Update derived function
      _phiInit->resetParameters({x0, y0, R, _epsilon});

      return FE_STATUS_OK;
    }

    //
    //
    //
    feStatus computePostProcessing(AdaptData &/*adapter*/,
                                   const std::vector<feSpace*> &spaces,
                                   feMesh */*mesh*/,
                                   feSolution *sol,
                                   PostProcData &data) override
    {
      const double time = sol->getCurrentTime();
      const feSpace *u   = spaces[0];
      const feSpace *p   = spaces[1];
      const feSpace *phi = spaces[2];

      // Velocity magnitude
      feNorm *L2_U;
      feCheckReturn(createNorm(L2_U, L2, {u}, sol));
      const double uNorm = L2_U->compute();
      delete L2_U;

      // Pressure at (0,0) and (-2., 0.)
      feNorm *probe = nullptr;
      feCheck(createNorm(probe, NONE, {p}, sol));
      std::vector<double> xIn  = { 0., 0., 0.};
      std::vector<double> xOut = {-2., 0., 0.};
      double pIn, pOut;
      feCheckReturn(probe->probeScalarField(xIn, pIn));
      feCheckReturn(probe->probeScalarField(xOut, pOut));
      delete probe;

      // Average of inner and outer pressures
      feNorm *bubbleArea;
      feCheckReturn(createNorm(bubbleArea, INTEGRAL_SOL_DEPENDENT_F, {phi}, sol, _indicatorArea));
      const double areaIn  = bubbleArea->compute();
      const double areaOut = 25. - areaIn; // Box is 5 x 5
      delete bubbleArea;

      feNorm *avgPin, *avgPout;
      feCheckReturn(createNorm( avgPin, INTEGRAL_SOL_DEPENDENT_F, {phi, p}, sol, _indicatorPin));
      feCheckReturn(createNorm(avgPout, INTEGRAL_SOL_DEPENDENT_F, {phi, p}, sol, _indicatorPout));
      const double avgPin_val  =  avgPin->compute() / areaIn;
      const double avgPout_val = avgPout->compute() / areaOut;
      delete avgPin;
      delete avgPout;

      data.res = {time, uNorm, pIn, pOut, avgPin_val, avgPout_val};

      // Reconstruct the interface (isoline phi = 0)
      feNorm *interface;
      feCheckReturn(createNorm(interface, NONE, {phi}, sol));
      const double targetValue = 0.;
      feCheckReturn(interface->reconstructScalarIsoline(targetValue, data.isoline));
      delete interface;

      return FE_STATUS_OK;
    }
  };
} // namespace CHNS_YoungLaplace

namespace CHNS_Jurin
{
  double phiInit_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    double y = args.pos[1];
    double epsilon = par[0];
    return tanh((y - 4.)/(sqrt(2.) * epsilon));
  }

  class VolumeAveragedSolver  : public CHNS_Generic::SolverCHNS_VolumeAveraged
  {
  public:
    double _contactAngle = 50.;

  public:
    VolumeAveragedSolver() : CHNS_Generic::SolverCHNS_VolumeAveraged()
    {
      _readDir     = "../data/CHNS/jurin/";
      _writeDir    = "../data/CHNS/jurin/";
      // _initialMesh = "../data/CHNS/jurin/coarse.msh";
      _initialMesh = "../data/CHNS/jurin/medium.msh";
      // _initialMesh = "../data/CHNS/jurin/fine.msh";
    };

    // Parse for CHNS parameters
    feStatus parseAndCreateFunctions(int argc, char **argv)
    {
      feCheckReturn(parseGenericCHNSParameters(argc, argv));

      // Parse specific test case parameters if needed
      feOptionsParser options(argc, argv, true);
      options.addOption(&_contactAngle, "-angle", "--angle", "");
      feCheckReturn(options.parse());

      const double m_to_mm  = 1e3;
      const double m_to_mm2 = 1e6;
      const double m_to_mm3 = 1e9;

      _rhoA          =    1. / m_to_mm3;
      _rhoB          = 2000. / m_to_mm3;

      const double nuA = 8e-5  * m_to_mm2;
      const double nuB = 1e-4  * m_to_mm2;

      _viscA         = nuA * _rhoA;
      _viscB         = nuB * _rhoB;

      _gamma         = 0.073;

      _gravityVal    = -9.81 * m_to_mm;

      // const double mesh_size = 1./32;
      // _epsilon0 = 0.64 * mesh_size;
      // _epsilon = _epsilon0;
      // _M = 0.1 * _epsilon * _epsilon;

      _epsilon = 0.044;
      _M       = 1e-7;

      createCHNSFunctions();

      _phiInit           = new feFunction(phiInit_f, {_epsilon});
      _volumeForce       = new feConstantVectorFunction({0., _gravityVal});

      return FE_STATUS_OK;
    }

    // Update specific test case parameters as the mesh is adapted
    feStatus updateParameters(const int fixedPointIteration) override
    {
      SolverCHNS_VolumeAveraged::updateParameters(fixedPointIteration);

      // Update derived function
      _phiInit->resetParameters({_epsilon});

      return FE_STATUS_OK;
    }

    feStatus createSpaces(feMesh *mesh,
                          std::vector<feSpace*> &spaces,
                          std::vector<feSpace*> &essentialSpaces)
    {
      for(auto *s : spaces) delete s;
      spaces.clear();
      essentialSpaces.clear();

      feSpace *u, *uSlipX, *uSlipY, *uSlipYcontact, *uSlipYMeasure, *uSlipSymmetry, *uTop,
        *p, *pPoint, *phi, *phiContact, *phiMeasure, *phiSymmetry, *mu, *muContact;

      feCheckReturn(createFiniteElementSpace(            u, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",       "Domaine", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(       uSlipX, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",         "SlipX", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(       uSlipY, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",         "SlipY", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(uSlipYcontact, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",  "ContactAngle", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(uSlipYMeasure, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",   "MeasureSide", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(uSlipSymmetry, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",      "Symmetry", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(         uTop, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",           "Top", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(            p, mesh, elementType::LAGRANGE,        _orderPressure,   "P",       "Domaine", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(       pPoint, mesh, elementType::LAGRANGE,        _orderPressure,   "P", "PressurePoint", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(          phi, mesh, elementType::LAGRANGE,        _orderPhi,        "Phi",     "Domaine", _degreeQuadrature, _phiInit));
      feCheckReturn(createFiniteElementSpace(   phiContact, mesh, elementType::LAGRANGE,        _orderPhi,        "Phi","ContactAngle", _degreeQuadrature, _phiInit));
      feCheckReturn(createFiniteElementSpace(   phiMeasure, mesh, elementType::LAGRANGE,        _orderPhi,        "Phi", "MeasureSide", _degreeQuadrature, _phiInit));
      feCheckReturn(createFiniteElementSpace(  phiSymmetry, mesh, elementType::LAGRANGE,        _orderPhi,        "Phi",    "Symmetry", _degreeQuadrature, _phiInit));
      feCheckReturn(createFiniteElementSpace(           mu, mesh, elementType::LAGRANGE,        _orderMu,         "Mu",      "Domaine", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(    muContact, mesh, elementType::LAGRANGE,        _orderMu,         "Mu", "ContactAngle", _degreeQuadrature, &scalarConstant::zero));

             uSlipY->setEssentialComponent(0, true);
      uSlipYcontact->setEssentialComponent(0, true);
      uSlipYMeasure->setEssentialComponent(0, true);
      uSlipSymmetry->setEssentialComponent(0, true);

      uSlipX->setEssentialComponent(1, true);
        uTop->setEssentialComponent(1, true);
      
      spaces = {u, p, phi, mu, uSlipX, uSlipY, uSlipYcontact, phiContact, muContact, uSlipYMeasure, phiMeasure, uSlipSymmetry, phiSymmetry, uTop, pPoint};
      essentialSpaces = {pPoint};

      return FE_STATUS_OK;
    }

    feStatus solve(AdaptData &adapter,
                   const int iInterval,
                   feMesh *mesh,
                   feSolution *sol,
                   feMetaNumber *numbering,
                   std::vector<feSpace*> &spaces,
                   const int /*numUnknowns*/,
                   feExportData &exportData)
    {
      std::vector<feBilinearForm*> forms = {};

      feSpace *u   = spaces[0];
      feSpace *p   = spaces[1];
      feSpace *phi = spaces[2];
      feSpace *mu  = spaces[3];
      feSpace *phiContact = spaces[7];
      feSpace *muContact  = spaces[8];

      std::vector<double> CHNSparameters = {_gamma, _epsilon};

      // CHNS
      feBilinearForm *CHNS;
      feCheck(createBilinearForm(CHNS, {u, p, phi, mu}, new CHNS_VolumeAveraged<2>(
        _density,
        _drhodphi,
        _viscosity,
        _dviscdphi,
        _mobility,
        _volumeForce,
        _pSource,
        _uSource,
        _phiSource,
        _muSource,
        CHNSparameters)));
      forms.push_back(CHNS);

      // Contact angle
      feBilinearForm *flux;
      feInfo("Contact angle = %f (deg)", _contactAngle);
      feCheck(createBilinearForm(flux, {muContact, phiContact, phi}, new feSysElm_FluxPrescribedNormalAngle<2>(&scalarConstant::one, _contactAngle * M_PI / 180.)));
      forms.push_back(flux);

      feLinearSystem *system = nullptr;
      #if defined(HAVE_MKL)
        feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
        feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering));
      #else
        feCheckReturn(createLinearSystem(system, PETSC, forms, numbering));
      #endif

      std::vector<feNorm*> norms = {};

      // timeIntegratorScheme scheme = timeIntegratorScheme::BDF1;
      timeIntegratorScheme scheme = timeIntegratorScheme::BDF2;

      TimeIntegrator *timeIntegrator = nullptr;
      const double t0 = adapter.currentTime;
      const double t1 = adapter.currentTime + adapter.dt * adapter.nTimeStepsPerInterval;
      const int nTimeSteps = adapter.nTimeStepsPerInterval;
      feCheckReturn(createTimeIntegrator(timeIntegrator, scheme, NLoptions,
        system, sol, mesh, norms, exportData, t0, t1, nTimeSteps));

      // Restart from solution container if not the first sub-interval
      // (restart from last solutions of previous sub-interval)
      if(iInterval > 0) {
        // BDF1Integrator *ptr = dynamic_cast<BDF1Integrator*> (timeIntegrator);
        // BDF2Integrator *ptr = dynamic_cast<BDF2Integrator*> (timeIntegrator);
        
        // ptr->_sC->NaNify();
        feCheckReturn(timeIntegrator->restartFromContainer(*adapter.allContainers[iInterval]));
        timeIntegrator->setCurrentStep(adapter.currentStep);

        feInfo("Time           of timeIntegrator : %f", timeIntegrator->getCurrentTime());
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[0]);
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[1]);
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[2]);
        feInfo("Step           of timeIntegrator : %d", timeIntegrator->getCurrentStep());
      }

      // Solve and compute metrics at each time step
      for(int iStep = 0; iStep < adapter.nTimeStepsPerInterval; ++iStep)
      {
        tic();
        feCheckReturn(timeIntegrator->makeSteps(1));
        feInfo("Computed time step in %f s", toc());

        // Reconstruct phi
        std::string meshFile = adapter.allOptions[iInterval].backgroundMeshfile;
        std::string recoveryFile = adapter.writeDir + "phi.msh";
        feNewRecovery recPhi(phi, 0, mesh, sol, meshFile, recoveryFile, false, false, false, nullptr, numbering);

        // Compute metric field and increment metric integral Hi
        feMetric *metricField = adapter.allMetrics[iInterval];
        computeAndIncrementMetric(timeIntegrator->getCurrentTime(), iStep, iInterval, adapter, metricField, &recPhi);

        adapter.currentStep++;
      }

      adapter.currentTime = timeIntegrator->getCurrentTime();

      delete adapter.allContainers[iInterval];
      adapter.allContainers[iInterval] = new feSolutionContainer();
      adapter.allContainers[iInterval]->NaNify();
      *(adapter.allContainers[iInterval]) = timeIntegrator->getSolutionContainer();

      delete timeIntegrator;
      delete system;
      for(feBilinearForm *f : forms)
        delete f;

      return FE_STATUS_OK;
    }

    //
    // Compute height between tip of meniscus and outer fluid
    // Export all quantities in millimeters
    //
    feStatus computePostProcessing(AdaptData &/*adapter*/,
                                   const std::vector<feSpace*> &spaces,
                                   feMesh */*mesh*/,
                                   feSolution *sol,
                                   PostProcData &data) override
    {
      const double time = sol->getCurrentTime();
      const feSpace *phi         = spaces[2];
      const feSpace *phiContact  = spaces[7];
      const feSpace *phiMeasure  = spaces[10];
      const feSpace *phiSymmetry = spaces[12];

      // Detect point where phi = 0 on contact angle boundary and on opposite boundary
      feNorm *interfaceContact, *interfaceOpposite, *interfaceSymmetry;
      feCheckReturn(createNorm(interfaceContact, NONE, {phiContact}, sol));
      feCheckReturn(createNorm(interfaceOpposite, NONE, {phiMeasure}, sol));
      feCheckReturn(createNorm(interfaceSymmetry, NONE, {phiSymmetry}, sol));
      const double targetValue = 0.;

      double measuredContactAngle = -1.;

      // Height at r = R
      // Extract height from isoline (y-component of only point in isoline)
      data.isoline.clear();
      feCheckReturn(interfaceContact->reconstructScalarIsoline(targetValue, data.isoline));
      const double hTipR = data.isoline[1];

      // Height at r = 0
      data.isoline.clear();
      feCheckReturn(interfaceSymmetry->reconstructScalarIsoline(targetValue, data.isoline));
      const double hTip0 = data.isoline[1];

      data.isoline.clear();
      feCheckReturn(interfaceOpposite->reconstructScalarIsoline(targetValue, data.isoline));
      const double hLow = data.isoline[1];

      const double deltaH    = 0.5 * (hTip0 + hTipR) - hLow;
      const double R         = 1e-3; // Half width of column
      const double rho       = _rhoB * 1e9;
      const double deltaHref = _gamma * cos(_contactAngle * M_PI / 180.) / (rho * 9.81 * R);

      // Reconstruct the interface (isoline phi = 0)
      feNorm *interface;
      feCheckReturn(createNorm(interface, NONE, {phi}, sol));
      feCheckReturn(interface->reconstructScalarIsoline(targetValue, data.isoline, &measuredContactAngle));

      const double measured_deltaHref = _gamma * cos(measuredContactAngle * M_PI / 180.) / (rho * 9.81 * R);

      data.res = {time, hTipR, hTip0, hLow, deltaH, deltaHref * 1e3, measuredContactAngle, measured_deltaHref * 1e3};

      delete interface;
      delete interfaceContact;
      delete interfaceOpposite;
      delete interfaceSymmetry;

      return FE_STATUS_OK;
    }
  }; // VolumeAveragedSolver

  class MassAveragedSolver  : public CHNS_Generic::GenericSolverCHNS_Alternative
  {
  public:
    double _contactAngle = 50.;

  public:
    MassAveragedSolver() : CHNS_Generic::GenericSolverCHNS_Alternative()
    {
      _readDir     = "../data/CHNS/jurin/";
      _writeDir    = "../data/CHNS/jurin/";
      // _initialMesh = "../data/CHNS/jurin/coarse.msh";
      _initialMesh = "../data/CHNS/jurin/medium.msh";
      // _initialMesh = "../data/CHNS/jurin/fine.msh";
    };

    // Parse for CHNS parameters
    feStatus parseAndCreateFunctions(int argc, char **argv)
    {
      feCheckReturn(parseGenericCHNSParameters(argc, argv));

      // Parse specific test case parameters if needed
      feOptionsParser options(argc, argv, true);
      options.addOption(&_contactAngle, "-angle", "--angle", "");
      feCheckReturn(options.parse());

      const double m_to_mm  = 1e3;
      const double m_to_mm2 = 1e6;
      const double m_to_mm3 = 1e9;

      _rhoA          =    1. / m_to_mm3;
      _rhoB          = 2000. / m_to_mm3;

      const double nuA = 8e-5  * m_to_mm2;
      const double nuB = 1e-4  * m_to_mm2;

      _viscA         = nuA * _rhoA;
      _viscB         = nuB * _rhoB;

      _gamma         = 0.073;

      _gravityVal    = -9.81 * m_to_mm;

      // const double mesh_size = 1./32;
      // _epsilon0 = 0.64 * mesh_size;
      // _epsilon = _epsilon0;
      // _M = 0.1 * _epsilon * _epsilon;

      _epsilon = 0.044;
      _M       = 1e-7;

      createCHNSFunctions();

      _phiInit           = new feFunction(phiInit_f, {_epsilon});
      _volumeForce       = new feConstantVectorFunction({0., _gravityVal});

      return FE_STATUS_OK;
    }

    // Update specific test case parameters as the mesh is adapted
    feStatus updateParameters(const int fixedPointIteration) override
    {
      GenericSolverCHNS_Alternative::updateParameters(fixedPointIteration);

      // Update derived function
      _phiInit->resetParameters({_epsilon});

      return FE_STATUS_OK;
    }

    feStatus createSpaces(feMesh *mesh,
                          std::vector<feSpace*> &spaces,
                          std::vector<feSpace*> &essentialSpaces)
    {
      for(auto *s : spaces) delete s;
      spaces.clear();
      essentialSpaces.clear();

      feSpace *u, *uSlipX, *uSlipY, *uSlipYcontact, *uSlipYMeasure,
        *p, *pPoint, *phi, *phiContact, *phiMeasure, *mu, *muContact;

      feCheckReturn(createFiniteElementSpace(            u, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",       "Domaine", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(       uSlipX, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",         "SlipX", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(       uSlipY, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",         "SlipY", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(uSlipYcontact, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",  "ContactAngle", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(uSlipYMeasure, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",   "MeasureSide", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(            p, mesh, elementType::LAGRANGE,        _orderPressure,   "P",       "Domaine", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(       pPoint, mesh, elementType::LAGRANGE,        _orderPressure,   "P", "PressurePoint", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(          phi, mesh, elementType::LAGRANGE,        _orderPhi,        "Phi",     "Domaine", _degreeQuadrature, _phiInit));
      feCheckReturn(createFiniteElementSpace(   phiContact, mesh, elementType::LAGRANGE,        _orderPhi,        "Phi","ContactAngle", _degreeQuadrature, _phiInit));
      feCheckReturn(createFiniteElementSpace(   phiMeasure, mesh, elementType::LAGRANGE,        _orderPhi,        "Phi", "MeasureSide", _degreeQuadrature, _phiInit));
      feCheckReturn(createFiniteElementSpace(           mu, mesh, elementType::LAGRANGE,        _orderMu,         "Mu",      "Domaine", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(    muContact, mesh, elementType::LAGRANGE,        _orderMu,         "Mu", "ContactAngle", _degreeQuadrature, &scalarConstant::zero));

      uSlipY->setEssentialComponent(0, true);
      uSlipYcontact->setEssentialComponent(0, true);
      uSlipYMeasure->setEssentialComponent(0, true);
      uSlipX->setEssentialComponent(1, true);
      spaces = {u, p, phi, mu, uSlipX, uSlipY, uSlipYcontact, phiContact, muContact, uSlipYMeasure, phiMeasure};
      essentialSpaces = {};

      return FE_STATUS_OK;
    }

    feStatus solve(AdaptData &adapter,
                   const int iInterval,
                   feMesh *mesh,
                   feSolution *sol,
                   feMetaNumber *numbering,
                   std::vector<feSpace*> &spaces,
                   const int /*numUnknowns*/,
                   feExportData &exportData)
    {
      std::vector<feBilinearForm*> forms = {};

      feSpace *u   = spaces[0];
      feSpace *p   = spaces[1];
      feSpace *phi = spaces[2];
      feSpace *mu  = spaces[3];
      feSpace *phiContact = spaces[7];
      feSpace *muContact  = spaces[8];

      std::vector<double> CHNSparameters = {_alpha, _gamma, _epsilon};

      // CHNS
      feBilinearForm *CHNS;
      feCheck(createBilinearForm(CHNS, {u, p, phi, mu}, new feSysElm_CHNS_Alternative<2>(
        _density,
        _drhodphi,
        _viscosity,
        _dviscdphi,
        _mobility,
        _volumeForce,
        _pSource,
        _uSource,
        _phiSource,
        _muSource,
        CHNSparameters)));
      forms.push_back(CHNS);

      // Contact angle
      feBilinearForm *flux;
      feInfo("Contact angle = %f (deg)", _contactAngle);
      feCheck(createBilinearForm(flux, {muContact, phiContact, phi}, new feSysElm_FluxPrescribedNormalAngle<2>(&scalarConstant::one, _contactAngle * M_PI / 180.)));
      forms.push_back(flux);

      feLinearSystem *system = nullptr;
      #if defined(HAVE_MKL)
        feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
        feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering));
      #else
        feCheckReturn(createLinearSystem(system, PETSC, forms, numbering));
      #endif

      std::vector<feNorm*> norms = {};

      // timeIntegratorScheme scheme = timeIntegratorScheme::BDF1;
      timeIntegratorScheme scheme = timeIntegratorScheme::BDF2;

      TimeIntegrator *timeIntegrator = nullptr;
      const double t0 = adapter.currentTime;
      const double t1 = adapter.currentTime + adapter.dt * adapter.nTimeStepsPerInterval;
      const int nTimeSteps = adapter.nTimeStepsPerInterval;
      feCheckReturn(createTimeIntegrator(timeIntegrator, scheme, NLoptions,
        system, sol, mesh, norms, exportData, t0, t1, nTimeSteps));

      // Restart from solution container if not the first sub-interval
      // (restart from last solutions of previous sub-interval)
      if(iInterval > 0) {
        // BDF1Integrator *ptr = dynamic_cast<BDF1Integrator*> (timeIntegrator);
        // BDF2Integrator *ptr = dynamic_cast<BDF2Integrator*> (timeIntegrator);
        
        // ptr->_sC->NaNify();
        feCheckReturn(timeIntegrator->restartFromContainer(*adapter.allContainers[iInterval]));
        timeIntegrator->setCurrentStep(adapter.currentStep);

        feInfo("Time           of timeIntegrator : %f", timeIntegrator->getCurrentTime());
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[0]);
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[1]);
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[2]);
        feInfo("Step           of timeIntegrator : %d", timeIntegrator->getCurrentStep());
      }

      // Solve and compute metrics at each time step
      for(int iStep = 0; iStep < adapter.nTimeStepsPerInterval; ++iStep)
      {
        tic();
        feCheckReturn(timeIntegrator->makeSteps(1));
        feInfo("Computed time step in %f s", toc());

        // Reconstruct phi
        std::string meshFile = adapter.allOptions[iInterval].backgroundMeshfile;
        std::string recoveryFile = adapter.writeDir + "phi.msh";
        feNewRecovery recPhi(phi, 0, mesh, sol, meshFile, recoveryFile, false, false, false, nullptr, numbering);

        // Compute metric field and increment metric integral Hi
        feMetric *metricField = adapter.allMetrics[iInterval];
        computeAndIncrementMetric(timeIntegrator->getCurrentTime(), iStep, iInterval, adapter, metricField, &recPhi);

        adapter.currentStep++;
      }

      adapter.currentTime = timeIntegrator->getCurrentTime();

      delete adapter.allContainers[iInterval];
      adapter.allContainers[iInterval] = new feSolutionContainer();
      adapter.allContainers[iInterval]->NaNify();
      *(adapter.allContainers[iInterval]) = timeIntegrator->getSolutionContainer();

      delete timeIntegrator;
      delete system;
      for(feBilinearForm *f : forms)
        delete f;

      return FE_STATUS_OK;
    }

    //
    // Compute height between tip of meniscus and outer fluid
    // Export all quantities in millimeters
    //
    feStatus computePostProcessing(AdaptData &/*adapter*/,
                                   const std::vector<feSpace*> &spaces,
                                   feMesh */*mesh*/,
                                   feSolution *sol,
                                   PostProcData &data) override
    {
      const double time = sol->getCurrentTime();
      const feSpace *phiContact = spaces[7];
      const feSpace *phiMeasure = spaces[10];

      // Detect point where phi = 0 on contact angle boundary and on opposite boundary
      feNorm *interfaceContact, *interfaceOpposite;
      feCheckReturn(createNorm(interfaceContact, NONE, {phiContact}, sol));
      feCheckReturn(createNorm(interfaceOpposite, NONE, {phiMeasure}, sol));
      const double targetValue = 0.;

      data.isoline.clear();
      feCheckReturn(interfaceContact->reconstructScalarIsoline(targetValue, data.isoline));

      // Extract height from isoline (y-component of only point in isoline)
      const double hTip = data.isoline[1];

      data.isoline.clear();
      feCheckReturn(interfaceOpposite->reconstructScalarIsoline(targetValue, data.isoline));
      const double hLow = data.isoline[1];

      const double deltaH    = hTip - hLow;
      const double R   = 1e-3; // Half width of column
      const double rho = _rhoB * 1e9;
      const double deltaHref = _gamma * cos(_contactAngle * M_PI / 180.) / (rho * 9.81 * R);

      data.res = {time, hTip, hLow, deltaH, deltaHref * 1e3};

      delete interfaceContact;
      delete interfaceOpposite;

      return FE_STATUS_OK;
    }
  }; // Mass-averaged solver

} // namespace CHNS_Jurin

namespace CHNS_BubbleContact
{
  double phiInit_f(const feFunctionArguments &args, const std::vector<double> &par)
  {
    double x       = args.pos[0];
    double y       = args.pos[1];
    double x0      = par[0];
    double y0      = par[1];
    double R       = par[2];
    double epsilon = par[3];
    double d = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
    return tanh((d - R)/(sqrt(2.)*epsilon));
  }

  class VolumeAveragedSolver  : public CHNS_Generic::SolverCHNS_VolumeAveraged
  {
  public:
    double x0 = 7.5;
    double y0 = 0.0;
    double R  = 2.5;
    double _contactAngle = 50.;

    feConstantFunction *_coeffFlux;

  public:
    VolumeAveragedSolver() : CHNS_Generic::SolverCHNS_VolumeAveraged()
    {
      _readDir     = "../data/CHNS/BubbleContact/";
      _writeDir    = "../data/CHNS/BubbleContact/";
      // _initialMesh = "../data/CHNS/BubbleContact/coarse.msh";
      _initialMesh = "../data/CHNS/BubbleContact/mesh.msh";
      // _initialMesh = "../data/CHNS/BubbleContact/fine.msh";
    };

    // Parse for CHNS parameters
    feStatus parseAndCreateFunctions(int argc, char **argv)
    {
      feCheckReturn(parseGenericCHNSParameters(argc, argv));

      // Parse specific test case parameters if needed
      feOptionsParser options(argc, argv, true);
      options.addOption(&_contactAngle, "-angle", "--angle", "");
      feCheckReturn(options.parse());

      const double m_to_mm  = 1e3;
      const double m_to_mm3 = 1e9;

      _rhoA        =   1.2 / m_to_mm3;
      _rhoB        = 1000. / m_to_mm3;

      _viscA       = 1.82e-5 / m_to_mm * 100.;
      _viscB       = 1.01e-3 / m_to_mm * 100.;

      _gamma       = 0.073;

      _gravityVal  = -9.81 * m_to_mm;

      _epsilon     = 0.4;
      _M           = 0.1 * _epsilon * _epsilon;

      createCHNSFunctions();

      _phiInit           = new feFunction(phiInit_f, {x0, y0, R, _epsilon});
      _volumeForce       = new feConstantVectorFunction({0., _gravityVal});


      const double lambda  = 3. / (2. * sqrt(2.)) * _gamma * _epsilon;
      _coeffFlux         = new feConstantFunction(lambda);

      return FE_STATUS_OK;
    }

    // Update specific test case parameters as the mesh is adapted
    feStatus updateParameters(const int fixedPointIteration) override
    {
      SolverCHNS_VolumeAveraged::updateParameters(fixedPointIteration);

      // Update derived function
      _phiInit->resetParameters({x0, y0, R, _epsilon});

      const double lambda  = 3. / (2. * sqrt(2.)) * _gamma * _epsilon;
      _coeffFlux->resetParameters({lambda});

      return FE_STATUS_OK;
    }

    feStatus createSpaces(feMesh *mesh,
                          std::vector<feSpace*> &spaces,
                          std::vector<feSpace*> &essentialSpaces)
    {
      for(auto *s : spaces) delete s;
      spaces.clear();
      essentialSpaces.clear();

      feSpace *u, *uT, *uB, *uL, *uR, *p, *pPoint, *phi, *phiB, *mu, *muB;

      feCheckReturn(createFiniteElementSpace(     u, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",       "Domaine", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(    uT, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",           "Top", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(    uB, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",        "Bottom", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(    uL, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",          "Left", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(    uR, mesh, elementType::VECTOR_LAGRANGE, _orderVelocity,   "U",         "Right", _degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace(     p, mesh, elementType::LAGRANGE,        _orderPressure,   "P",       "Domaine", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(pPoint, mesh, elementType::LAGRANGE,        _orderPressure,   "P", "PressurePoint", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(   phi, mesh, elementType::LAGRANGE,        _orderPhi,        "Phi",     "Domaine", _degreeQuadrature, _phiInit));
      feCheckReturn(createFiniteElementSpace(  phiB, mesh, elementType::LAGRANGE,        _orderPhi,        "Phi",      "Bottom", _degreeQuadrature, _phiInit));
      feCheckReturn(createFiniteElementSpace(    mu, mesh, elementType::LAGRANGE,        _orderMu,         "Mu",      "Domaine", _degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(   muB, mesh, elementType::LAGRANGE,        _orderMu,         "Mu",       "Bottom", _degreeQuadrature, &scalarConstant::zero));

      uL->setEssentialComponent(0, true);
      uR->setEssentialComponent(0, true);

      uT->setEssentialComponent(1, true);
      uB->setEssentialComponent(1, true);
      
      spaces = {u, p, phi, mu, uT, uB, uL, uR, pPoint, phiB, muB};
      essentialSpaces = {pPoint};

      return FE_STATUS_OK;
    }

    feStatus solve(AdaptData &adapter,
                   const int iInterval,
                   feMesh *mesh,
                   feSolution *sol,
                   feMetaNumber *numbering,
                   std::vector<feSpace*> &spaces,
                   const int /*numUnknowns*/,
                   feExportData &exportData)
    {
      std::vector<feBilinearForm*> forms = {};

      feSpace *u   = spaces[0];
      feSpace *p   = spaces[1];
      feSpace *phi = spaces[2];
      feSpace *mu  = spaces[3];
      feSpace *phiContact = spaces[9];
      feSpace *muContact  = spaces[10];

      std::vector<double> CHNSparameters = {_gamma, _epsilon};

      // CHNS
      feBilinearForm *CHNS;
      feCheck(createBilinearForm(CHNS, {u, p, phi, mu}, new CHNS_VolumeAveraged<2>(
        _density,
        _drhodphi,
        _viscosity,
        _dviscdphi,
        _mobility,
        _volumeForce,
        _pSource,
        _uSource,
        _phiSource,
        _muSource,
        CHNSparameters)));
      forms.push_back(CHNS);

      // Contact angle
      feBilinearForm *flux;
      feInfo("Contact angle = %f (deg)", _contactAngle);
      feCheck(createBilinearForm(flux, {muContact, phiContact, phi}, new feSysElm_FluxPrescribedNormalAngle<2>(_coeffFlux, _contactAngle * M_PI / 180.)));
      forms.push_back(flux);

      feLinearSystem *system = nullptr;
      #if defined(HAVE_MKL)
        feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
        feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering));
      #else
        feCheckReturn(createLinearSystem(system, PETSC, forms, numbering));
      #endif

      std::vector<feNorm*> norms = {};

      // timeIntegratorScheme scheme = timeIntegratorScheme::BDF1;
      timeIntegratorScheme scheme = timeIntegratorScheme::BDF2;

      TimeIntegrator *timeIntegrator = nullptr;
      const double t0 = adapter.currentTime;
      const double t1 = adapter.currentTime + adapter.dt * adapter.nTimeStepsPerInterval;
      const int nTimeSteps = adapter.nTimeStepsPerInterval;
      feCheckReturn(createTimeIntegrator(timeIntegrator, scheme, NLoptions,
        system, sol, mesh, norms, exportData, t0, t1, nTimeSteps));

      // Restart from solution container if not the first sub-interval
      // (restart from last solutions of previous sub-interval)
      if(iInterval > 0) {
        // BDF1Integrator *ptr = dynamic_cast<BDF1Integrator*> (timeIntegrator);
        // BDF2Integrator *ptr = dynamic_cast<BDF2Integrator*> (timeIntegrator);
        
        // ptr->_sC->NaNify();
        feCheckReturn(timeIntegrator->restartFromContainer(*adapter.allContainers[iInterval]));
        timeIntegrator->setCurrentStep(adapter.currentStep);

        feInfo("Time           of timeIntegrator : %f", timeIntegrator->getCurrentTime());
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[0]);
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[1]);
        feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[2]);
        feInfo("Step           of timeIntegrator : %d", timeIntegrator->getCurrentStep());
      }

      // Solve and compute metrics at each time step
      for(int iStep = 0; iStep < adapter.nTimeStepsPerInterval; ++iStep)
      {
        tic();
        feCheckReturn(timeIntegrator->makeSteps(1));
        feInfo("Computed time step in %f s", toc());

        // Reconstruct phi
        std::string meshFile = adapter.allOptions[iInterval].backgroundMeshfile;
        std::string recoveryFile = adapter.writeDir + "phi.msh";
        feNewRecovery recPhi(phi, 0, mesh, sol, meshFile, recoveryFile, false, false, false, nullptr, numbering);

        // Compute metric field and increment metric integral Hi
        feMetric *metricField = adapter.allMetrics[iInterval];
        computeAndIncrementMetric(timeIntegrator->getCurrentTime(), iStep, iInterval, adapter, metricField, &recPhi);

        adapter.currentStep++;
      }

      adapter.currentTime = timeIntegrator->getCurrentTime();

      delete adapter.allContainers[iInterval];
      adapter.allContainers[iInterval] = new feSolutionContainer();
      adapter.allContainers[iInterval]->NaNify();
      *(adapter.allContainers[iInterval]) = timeIntegrator->getSolutionContainer();

      delete timeIntegrator;
      delete system;
      for(feBilinearForm *f : forms)
        delete f;

      return FE_STATUS_OK;
    }

    //
    // Compute height between tip of meniscus and outer fluid
    // Export all quantities in millimeters
    //
    feStatus computePostProcessing(AdaptData &/*adapter*/,
                                   const std::vector<feSpace*> &spaces,
                                   feMesh */*mesh*/,
                                   feSolution *sol,
                                   PostProcData &data) override
    {
      const double time = sol->getCurrentTime();
      const feSpace *phi   = spaces[2];
      // const feSpace *phiB  = spaces[9];

      double measuredContactAngle = -1.;

      // Reconstruct the interface (isoline phi = 0)
      feNorm *interface;
      feCheckReturn(createNorm(interface, NONE, {phi}, sol));
      const double targetValue = 0.;
      feCheckReturn(interface->reconstructScalarIsoline(targetValue, data.isoline, &measuredContactAngle));

      data.res = {time, measuredContactAngle};

      delete interface;

      return FE_STATUS_OK;
    }
  }; // VolumeAveragedSolver
} // namespace CHNS_BubbleContact

namespace Cylinder {

  void uInlet_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
  {
    const double t = args.t;
    const double y = args.pos[1];
    const double h = par[0];
    const double t0 = par[1];
    const double ramp_t = fmin(t/t0, 1.);
    // const double epsilon = (y > h/2.)  ? 0.1 * (1. - ramp_t) : 0.;
    const double epsilon = 0.1 * y/h * (1. - ramp_t);
    res[0] = 1. + fabs(epsilon);
    // res[0] = 1.;
    res[1] = 0.;
  }

  class Solver  : public GenericSolver
  {
  public:
    std::string readDir = "../data/cylinderUnsteadyMetric/";
    std::string writeDir = "../data/cylinderUnsteadyMetric/";
    const char *initialMesh = "../data/cylinderUnsteadyMetric/cylinderMedium.msh";

    bool computeTemperature = true;

    int orderV = 2;
    int orderP = 1;
    int orderW = 1;
    int orderT = 1;
    int degreeQuadrature = 10;
    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 50, 4, 1.};

    double D = 1.;   // Diameter
    double h = 16.;   // Geometry height
    double tRamp = 0.5; // Time for the ramp on the initial condition
    double Re = 250.;

    feVectorFunction *fInlet = nullptr;
    feConstantFunction *density = nullptr, *viscosity = nullptr;
    feConstantFunction *Tzero = nullptr, *Tcylindre = nullptr, *diffusivity = nullptr;

  public:
    Solver() {};
    int getNumFields() { return 4; }
    std::string getReadDir(){ return readDir; };
    std::string getWriteDir(){ return writeDir; };
    std::string getInitialMesh(){ return std::string(initialMesh); };
    int getOrderForAdaptation() { return orderW; };

    feStatus parseAndCreateFunctions(int argc, char **argv)
    {
      bool ignoreUnrecognizedOptions = true;
      feOptionsParser options(argc, argv, ignoreUnrecognizedOptions);
      // options.addOption(&writeDir, "-dir", "--dir", "");
      // options.addOption(&initialMesh, "-m", "--m", "");
      options.addOption(&orderV, "-ou", "--orderU", "");
      options.addOption(&orderP, "-op", "--orderP", "");
      options.addOption(&h, "-h", "--h", "");
      options.addOption(&D, "-D", "--D", "");
      options.addOption(&tRamp, "-tRamp", "--tRamp", "");
      feCheckReturn(options.parse());

      fInlet = new feVectorFunction(uInlet_f, {h, tRamp});
      density = new feConstantFunction(1.);
      viscosity = new feConstantFunction(D/Re);
      Tzero     = new feConstantFunction(293.);
      Tcylindre = new feConstantFunction(373.);
      diffusivity = new feConstantFunction(0.01);

      return FE_STATUS_OK;
    }

    // Cylinder
    feStatus createSpaces(feMesh *mesh,
                          std::vector<feSpace*> &spaces,
                          std::vector<feSpace*> &essentialSpaces)
    {
      spaces.clear();
      essentialSpaces.clear();
      feSpace *u, *uSlip, *uNoSlip, *uInlet, *p, *w, *T, *Tcyl;
      feCheckReturn(createFiniteElementSpace(      u, mesh, elementType::VECTOR_LAGRANGE, orderV, "U", "Domaine", degreeQuadrature, &vectorConstant::one_x));
      feCheckReturn(createFiniteElementSpace(  uSlip, mesh, elementType::VECTOR_LAGRANGE, orderV, "U",    "Slip", degreeQuadrature, &vectorConstant::one_x));
      feCheckReturn(createFiniteElementSpace(uNoSlip, mesh, elementType::VECTOR_LAGRANGE, orderV, "U",  "NoSlip", degreeQuadrature, &vectorConstant::zero));
      feCheckReturn(createFiniteElementSpace( uInlet, mesh, elementType::VECTOR_LAGRANGE, orderV, "U",   "Inlet", degreeQuadrature, fInlet));
      feCheckReturn(createFiniteElementSpace(      p, mesh, elementType::LAGRANGE,        orderP, "P", "Domaine", degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(      w, mesh, elementType::LAGRANGE,        orderW, "W", "Domaine", degreeQuadrature, &scalarConstant::zero));
      feCheckReturn(createFiniteElementSpace(      T, mesh, elementType::LAGRANGE,        orderT, "T", "Domaine", degreeQuadrature, Tzero));
      feCheckReturn(createFiniteElementSpace(   Tcyl, mesh, elementType::LAGRANGE,        orderT, "T",  "NoSlip", degreeQuadrature, Tcylindre));
      uSlip->setEssentialComponent(1, true);
      spaces = {u, uSlip, uNoSlip, uInlet, p, w, T, Tcyl};
      essentialSpaces = {uNoSlip, uInlet, Tcyl};

      if(computeTemperature) {
        spaces.push_back(T);
        spaces.push_back(Tcyl);
        essentialSpaces.push_back(Tcyl);
      }

      return FE_STATUS_OK;
    }

    // Cylinder
    feStatus solve(AdaptData &adapter,
                   const int iInterval,
                   feMesh *mesh,
                   feSolution *sol,
                   feMetaNumber *numbering,
                   std::vector<feSpace*> &spaces,
                   const int /*numUnknowns*/,
                   feExportData &exportData)
    {
      std::vector<feBilinearForm*> forms = {};

      feSpace *u = spaces[0];
      feSpace *p = spaces[4];
      feSpace *w = spaces[5];
      feSpace *T = spaces[6];

      // Momentum
      feBilinearForm *dudt, *conv, *gradp, *diffu;
      feCheckReturn(createBilinearForm( dudt,    {u}, new feSysElm_TransientVectorMass<2>(&scalarConstant::one)));
      feCheckReturn(createBilinearForm( conv,    {u}, new feSysElm_VectorConvectiveAcceleration<2>(&scalarConstant::one)));
      feCheckReturn(createBilinearForm(gradp, {u, p}, new feSysElm_MixedGradient<2>(&scalarConstant::one)));
      feCheckReturn(createBilinearForm(diffu,    {u}, new feSysElm_VectorDiffusion<2>(&scalarConstant::one, viscosity)));
      forms.push_back(dudt);
      forms.push_back(conv);
      forms.push_back(gradp);
      forms.push_back(diffu);

      // Continuity
      feBilinearForm *divU;
      feCheckReturn(createBilinearForm(divU, {p, u}, new feSysElm_MixedDivergence<2>(&scalarConstant::one)));
      forms.push_back(divU);

      // Vorticity : w = curl(u)
      feBilinearForm *curl_u, *mass_w;
      feCheckReturn(createBilinearForm(curl_u, {w, u}, new feSysElm_MixedCurl<2>(&scalarConstant::minusOne)));
      feCheckReturn(createBilinearForm(mass_w,    {w}, new feSysElm_Mass(&scalarConstant::one)));
      forms.push_back(curl_u);
      forms.push_back(mass_w);

      // Energy
      feBilinearForm *dTdt, *convT, *diffT;
      if(computeTemperature) {
        feCheckReturn(createBilinearForm( dTdt,    {T}, new feSysElm_TransientMass(&scalarConstant::one)));
        feCheckReturn(createBilinearForm(convT, {T, u}, new feSysElm_TracerConvection<2>(&scalarConstant::one)));
        feCheckReturn(createBilinearForm(diffT,    {T}, new feSysElm_Diffusion<2>(diffusivity)));
        forms.push_back(dTdt);
        forms.push_back(convT);
        forms.push_back(diffT);
      }

      for(feBilinearForm *form : forms)
        form->setComputeMatrixWithFD(false);

      feLinearSystem *system;
      #if defined(HAVE_MKL)
        feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
        feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering));
      #else
        feCheckReturn(createLinearSystem(system, PETSC, forms, numbering));
      #endif

      std::vector<feNorm*> norms = {};

      TimeIntegrator *timeIntegrator;
      const double t0 = adapter.currentTime;
      const double t1 = adapter.currentTime + adapter.dt * adapter.nTimeStepsPerInterval;
      const int nTimeSteps = adapter.nTimeStepsPerInterval;
      feCheckReturn(createTimeIntegrator(timeIntegrator, timeIntegratorScheme::BDF2, NLoptions,
        system, sol, mesh, norms, exportData, t0, t1, nTimeSteps));

      // Restart from solution container if not the first sub-interval
      // (restart from last solutions of previous sub-interval)
      if(iInterval > 0) {
        feCheckReturn(timeIntegrator->restartFromContainer(*adapter.allContainers[iInterval]));
        timeIntegrator->setCurrentStep(adapter.currentStep);
      }

      // Solve and compute metrics at each time step
      for(int iStep = 0; iStep < adapter.nTimeStepsPerInterval; ++iStep)
      {
        feCheckReturn(timeIntegrator->makeSteps(1));

        // Reconstruct velocity
        std::string meshFile = adapter.allOptions[iInterval].backgroundMeshfile;
        std::string recoveryFile = adapter.writeDir + "w.msh";
        feNewRecovery rec(w, 0, mesh, sol, meshFile, recoveryFile, false, false, false, nullptr, numbering);

        // Compute metric field and increment metric integral Hi
        feMetric *metricField = adapter.allMetrics[iInterval];
        computeAndIncrementMetric(timeIntegrator->getCurrentTime(), iStep, iInterval, adapter, metricField, &rec);

        adapter.currentStep++;
      }

      adapter.currentTime = timeIntegrator->getCurrentTime();
      // currentStep = timeIntegrator->getCurrentStep();

      feInfo("Assigning container");

      delete adapter.allContainers[iInterval];
      adapter.allContainers[iInterval] = new feSolutionContainer();
      *(adapter.allContainers[iInterval]) = timeIntegrator->getSolutionContainer();

      delete timeIntegrator;
      delete system;
      for(feBilinearForm *f : forms)
        delete f;

      feInfo("Solve OK");

      return FE_STATUS_OK;
    }

    // Cylinder
    feStatus computeLastMetric(AdaptData &adapter,
                               const int iInterval,
                               feMesh *mesh,
                               feSolution *sol,
                               feMetaNumber *numbering,
                               std::vector<feSpace*> &spaces)
    {
      // feSpace *u = spaces[0];
      feSpace *w = spaces[5];

      // Reconstruct velocity
      std::string meshFile = adapter.allOptions[iInterval].backgroundMeshfile;
      std::string recoveryFile = adapter.writeDir + "w.msh";
      feNewRecovery rec(w, 0, mesh, sol, meshFile, recoveryFile, false, false, false, nullptr, numbering);

      // Compute "steady state" metric field
      feMetric *metricField = adapter.allMetrics[iInterval];
      metricField->setRecoveredFields({&rec});
      metricField->_options.debug = true;
      metricField->setMetricScaling(true);
      feCheckReturn(metricField->computeMetrics());

      return FE_STATUS_OK;
    }

    feStatus computeError(const int iIntervalReference,
                          AdaptData &adapter,
                          const std::vector<feSpace*> &spaces,
                          feMesh *mesh,
                          feSolution *sol,
                          std::vector<double> &errors)
    {
      UNUSED(iIntervalReference, adapter, spaces, mesh, sol);
      errors.resize(this->getNumFields(), 0.);
      return FE_STATUS_OK;
    }
  };
} // namespace Cylinder

namespace Analytic {

  #define D_PARAM 0.1

  double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double t = args.t;
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double c = 1.;
    const double d = D_PARAM;
    return tanh((2.*(x - c*t) - sin(5.*y))/d);
  }

  double dudt_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double t = args.t;
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double c = 1.;
    const double d = D_PARAM;
    const double cosh_val = cosh(- (2.*(x-c*t) - sin(5.*y))/d);
    const double sech = 1./cosh_val;
    // Negative sign as it appears as source term
    return - (- 2.*c/d * sech * sech);
  }

  double uSource_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double t = args.t;
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double c = 1.;
    const double d = D_PARAM;
    const double cosh_val = cosh(- (2.*(x-c*t) - sin(5.*y))/d);
    const double sech = 1./cosh_val;
    const double dudt = - 2.*c/d * sech * sech;

    const double valx = (2.0 * (x-c*t) -sin(y*5.0));
    const double uxx = 1.0/(d*d)*tanh(valx/d)*(pow(tanh(valx/d),2.0)-1.0)*8.0;
    const double uyy = -(sin(y*5.0)*(pow(tanh(valx/d),2.0)-1.0)*2.5E+1)/d+1.0/(d*d)*pow(cos(y*5.0),2.0)*tanh(valx/d)*(pow(tanh(valx/d),2.0)-1.0)*5.0E+1;
    const double lap_u = uxx + uyy;
    return - (dudt - lap_u);
  } 

  void uHessian_f(const feFunctionArguments &args,
                  const std::vector<double> &/*par*/,
                  std::vector<double> &res)
  {
    const double t = args.t;
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double c = 1.;
    const double d = D_PARAM;
    const double valx = (2.0 * (x-c*t) -sin(y*5.0));
    const double uxx = 1.0/(d*d)*tanh(valx/d)*(pow(tanh(valx/d),2.0)-1.0)*8.0;
    const double uxy = -1.0/(d*d)*cos(y*5.0)*tanh(valx/d)*(pow(tanh(valx/d),2.0)-1.0)*2.0E+1;
    const double uyy = -(sin(y*5.0)*(pow(tanh(valx/d),2.0)-1.0)*2.5E+1)/d+1.0/(d*d)*pow(cos(y*5.0),2.0)*tanh(valx/d)*(pow(tanh(valx/d),2.0)-1.0)*5.0E+1;
    res[0] = uxx;
    res[1] = uxy;
    res[2] = uxy;
    res[3] = uyy;
    return;
  }

  void uThird_f(const feFunctionArguments &args,
                const std::vector<double> &/*par*/,
                std::vector<double> &res)
  {
    const double t = args.t;
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double c = 1.;
    const double d = D_PARAM;
    const double valx = (2.0 * (x-c*t) -sin(y*5.0));
    const double valx2 = pow(tanh(valx/d),2.0);
    const double uxxx = 1.0/(d*d*d)*pow(valx2-1.0,2.0)*-1.6E+1-1.0/(d*d*d)*valx2*(valx2-1.0)*3.2E+1 + 1. * 6.0;
    const double uxxy = 1.0/(d*d*d)*cos(y*5.0)*pow(valx2-1.0,2.0)*4.0E+1+1.0/(d*d*d)*cos(y*5.0)*valx2*(valx2-1.0)*8.0E+1;
    const double uxyy = 1.0/(d*d*d)*pow(cos(y*5.0),2.0)*pow(valx2-1.0,2.0)*-1.0E+2+1.0/(d*d)*sin(y*5.0)*tanh(valx/d)*(valx2-1.0)*1.0E+2-1.0/(d*d*d)*pow(cos(y*5.0),2.0)*valx2*(valx2-1.0)*2.0E+2+  + 1. * 2.0;
    const double uyyy = 1.0/(d*d*d)*pow(cos(y*5.0),3.0)*pow(valx2-1.0,2.0)*2.5E+2-(cos(y*5.0)*(valx2-1.0)*1.25E+2)/d+1.0/(d*d*d)*pow(cos(y*5.0),3.0)*valx2*(valx2-1.0)*5.0E+2-1.0/(d*d)*cos(y*5.0)*sin(y*5.0)*tanh(valx/d)*(valx2-1.0)*7.5E+2;
    res[0] = uxxx;
    res[1] = uxxy;
    res[2] = uxxy;
    res[3] = uxyy;
    res[4] = uxxy;
    res[5] = uxyy;
    res[6] = uxyy;
    res[7] = uyyy;
    return;
  }

  class Solver  : public GenericSolver
  {
  public:
    int order = 1;
    int degreeQuadrature = 8;

    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 50, 4, 1.};

    std::string readDir = "../data/unsteadyMetricAnalytic/verificationP1_referenceSolution/";
    std::string writeDir = "../data/unsteadyMetricAnalytic/verificationP1_referenceSolution/";
    // const char *initialMesh = "../data/square4.msh";
    // const char *initialMesh = "../data/unsteadyMetricAnalytic/start.msh";
    const char *initialMesh = "../data/unsteadyMetricAnalytic/startFine.msh";

    feFunction uSol, duSoldt, uSource;
    feVectorFunction hessienne, troisieme;
  public:
    Solver() : 
      uSol(feFunction(uSol_f)),
      duSoldt(feFunction(dudt_f)),
      uSource(feFunction(uSource_f)),
      hessienne(feVectorFunction(uHessian_f)),
      troisieme(feVectorFunction(uThird_f))
    {

    };

    int getNumFields() { return 1; }
    std::string getReadDir(){ return readDir; };
    std::string getWriteDir(){ return writeDir; };
    std::string getInitialMesh(){ return std::string(initialMesh); };
    int getOrderForAdaptation() { return order; };

    feStatus parseAndCreateFunctions(int /*argc*/, char **/*argv*/)
    {
      return FE_STATUS_OK;
    };

    feStatus createSpaces(feMesh *mesh,
                          std::vector<feSpace*> &spaces,
                          std::vector<feSpace*> &essentialSpaces)
    {
      for(auto *s : spaces)
        delete s;
      spaces.clear();
      essentialSpaces.clear();
      feSpace *u = nullptr, *uB = nullptr;
      feCheckReturn(createFiniteElementSpace( u, mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &uSol));
      feCheckReturn(createFiniteElementSpace(uB, mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol));
      spaces = {u, uB};
      if(useAnalyticalSolution) {
        essentialSpaces = {u, uB};
      } else {
        essentialSpaces = {uB};
      }
      return FE_STATUS_OK;
    }

    // Analytic
    feStatus solve(AdaptData &adapter,
                   const int iInterval,
                   feMesh *mesh,
                   feSolution *sol,
                   feMetaNumber *numbering,
                   std::vector<feSpace*> &spaces,
                   const int /*numUnknowns*/,
                   feExportData &exportData)
    {
      feSpace *u = spaces[0];
      std::vector<feBilinearForm*> forms = {};
      feBilinearForm *dudt = nullptr, *diffu = nullptr, *src = nullptr;
      feLinearSystem *system = nullptr;
      std::vector<feNorm*> norms = {};
      TimeIntegrator *timeIntegrator = nullptr;

      if(!useAnalyticalSolution)
      {
        // // Solve dudt = duSoldt(x,t)
        // feCheckReturn(createBilinearForm(dudt, {u}, new feSysElm_TransientMass(&scalarConstant::one)));
        // feCheckReturn(createBilinearForm( src, {u}, new feSysElm_Source(&duSoldt)));
        // forms.push_back(dudt);
        // forms.push_back(src);

        // Solve dudt - div(grad(u)) + uSource = 0
        feCheckReturn(createBilinearForm( dudt, {u}, new feSysElm_TransientMass(&scalarConstant::one)));
        feCheckReturn(createBilinearForm(diffu, {u}, new feSysElm_Diffusion<2>(&scalarConstant::one)));
        feCheckReturn(createBilinearForm(  src, {u}, new feSysElm_Source(&uSource)));
        forms.push_back(dudt);
        forms.push_back(diffu);
        forms.push_back(src);

        #if defined(HAVE_MKL)
          feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering));
        #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
          feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering));
        #else
          feCheckReturn(createLinearSystem(system, PETSC, forms, numbering));
        #endif

        const double t0 = adapter.currentTime;
        const double t1 = adapter.currentTime + adapter.dt * adapter.nTimeStepsPerInterval;
        const int nTimeSteps = adapter.nTimeStepsPerInterval;

        feCheckReturn(createTimeIntegrator(timeIntegrator, timeIntegratorScheme::BDF2, NLoptions,
          system, sol, mesh, norms, exportData, t0, t1, nTimeSteps));

        // Start BDF2 with manufactured solution
        static_cast<BDF2Integrator*>(timeIntegrator)->setStartingMethod(BDF2Starter::InitialCondition);

        // Restart from solution container if not the first sub-interval
        // (restart from last solutions of previous sub-interval)
        if(iInterval > 0) {
          BDF2Integrator *ptr = static_cast<BDF2Integrator*> (timeIntegrator);
          ptr->_sC->NaNify();
          feCheckReturn(timeIntegrator->restartFromContainer(*adapter.allContainers[iInterval]));
          timeIntegrator->setCurrentStep(adapter.currentStep);

          feInfo("Time           of timeIntegrator : %f", timeIntegrator->getCurrentTime());
          feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[0]);
          feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[1]);
          feInfo("Previous times of timeIntegrator : %f", timeIntegrator->getTime()[2]);
          feInfo("Step           of timeIntegrator : %d", timeIntegrator->getCurrentStep());
        }
      }

      // Solve and compute metrics at each time step
      for(int iStep = 0; iStep < adapter.nTimeStepsPerInterval; ++iStep)
      {
        if(useAnalyticalSolution)
        {
          // Simply initialize the essential space
          sol->setCurrentTime(adapter.currentTime);
          sol->initialize(mesh);

          // if(adapter.exportSolution) {
          //   std::string vtkfile = exportData.fileNameRoot + std::to_string(adapter.currentStep) + ".vtk";
          //   feCheckReturn(exportData.exporter->writeStep(vtkfile));
          // }

        } else {
          feCheckReturn(timeIntegrator->makeSteps(1));
        }

        if(adapter.adapt)
        {
          // Compute metric field and increment metric integral Hi
          feMetric *metricField = adapter.allMetrics[iInterval];
          metricField->_options.useAnalyticDerivatives = true;
          metricField->_options.secondDerivatives = &hessienne;
          metricField->_options.thirdDerivatives = &troisieme;
          metricField->setCurrentTime(adapter.currentTime);
          metricField->setMetricScaling(false);
          metricField->computeMetrics();

          // Trapeze rule
          if(iStep == 0 || iStep == adapter.nTimeStepsPerInterval-1)
            metricField->addMetricsToOther(adapter.dt/2., adapter.allHi[iInterval]);
          else
            metricField->addMetricsToOther(adapter.dt   , adapter.allHi[iInterval]);
        }

        adapter.currentTime += adapter.dt;
        adapter.currentStep++;
      }

      delete adapter.allContainers[iInterval];
      if(useAnalyticalSolution)
      {
        // A container with a single solution
        adapter.allContainers[iInterval] = new feSolutionContainer(1, adapter.currentTime, numbering->getNbDOFs());
        adapter.allContainers[iInterval]->setCurrentSolution(sol);
        adapter.allContainers[iInterval]->setCurrentSolutionDot(sol);
      } else {
        // Get the container from the integrator
        adapter.allContainers[iInterval] = new feSolutionContainer();
        adapter.allContainers[iInterval]->NaNify();
        *(adapter.allContainers[iInterval]) = timeIntegrator->getSolutionContainer();

        delete timeIntegrator;
        delete system;
        for(feBilinearForm *f : forms)
          delete f;
      }

      return FE_STATUS_OK;

      // // Steady solve at t = 0
      // sol->setCurrentTime(0.);
      // sol->initialize(mesh);

      // if(adapter.exportSolution) {
      //   std::string vtkfile = exportData.fileNameRoot + std::to_string(adapter.currentStep) + ".vtk";
      //   feCheckReturn(exportData.exporter->writeStep(vtkfile));
      // }

      // // Compute metric field and increment metric integral Hi
      // feMetric *metricField = adapter.allMetrics[iInterval];
      // metricField->_options.useAnalyticDerivatives = true;
      // metricField->_options.secondDerivatives = &hessienne;
      // metricField->_options.thirdDerivatives = &troisieme;
      // metricField->setCurrentTime(0.);

      // delete adapter.allContainers[iInterval];
      // adapter.allContainers[iInterval] = new feSolutionContainer(1, 0., numDOF);
      // adapter.allContainers[iInterval]->setCurrentSolution(sol);
      // adapter.allContainers[iInterval]->setCurrentSolutionDot(sol);

      // return FE_STATUS_OK;
    }

    feStatus computeLastMetric(AdaptData &adapter,
                               const int iInterval,
                               feMesh */*mesh*/,
                               feSolution */*sol*/,
                               feMetaNumber */*numbering*/,
                               std::vector<feSpace*> &/*spaces*/)
    {
      // Compute "steady state" metric field
      feMetric *metricField = adapter.allMetrics[iInterval];
      metricField->_options.debug = true;
      metricField->setMetricScaling(true);
      feCheckReturn(metricField->computeMetrics());
      return FE_STATUS_OK;
    }

    // Compute error at the end of the sub-interval (at t_i+1)
    feStatus computeError(const int iIntervalReference,
                          AdaptData &adapter,
                          const std::vector<feSpace*> &spaces,
                          feMesh *mesh,
                          feSolution *sol,
                          std::vector<double> &errors)
    {
      // if(useAnalyticalSolution)
      // {
        UNUSED(iIntervalReference, adapter, mesh);

        feInfo("Computing error for analytical test case based on EXACT SOLUTION");
        const feSpace *u = spaces[0];
        feNorm *L2error;
        feCheckReturn(createNorm(L2error, L2_ERROR, {u}, sol, &uSol));
        errors.resize(1);
        errors[0] = L2error->compute();
        return FE_STATUS_OK;
      // }
      // else
      // {
      //   feInfo("Computing error for CHNS based on REFERENCE FINER SOLUTION");

      //   errors.resize(1, 0.);

      //   const feSpace *u   = spaces[0];

      //   feNorm *errorOnScalarSpaces = nullptr;
      //   feCheckReturn(createNorm(errorOnScalarSpaces, NONE, {u}, sol));
      //   feCheckReturn(errorOnScalarSpaces->computeErrorNormFromExternalSolution(adapter.referenceTestCase[iIntervalReference].mesh,
      //                                                                           adapter.referenceTestCase[iIntervalReference].sol,
      //                                                                           adapter.referenceTestCase[iIntervalReference].spaces,
      //                                                                           mesh, sol, errors));
      //   feInfo("Computed error on u  : %+-1.6e", errors[0]);
      //   delete errorOnScalarSpaces;
      //   return FE_STATUS_OK;
      // }
    }
  };
} // namespace Analytic

int main(int argc, char **argv)
{
  initialize(argc, argv);

  const int verbosity = 2;
  setVerbose(verbosity);

  const char *solverType = "";

  feOptionsParser options(argc, argv, true);
  options.addOption(&solverType, "-solver", "--solver", "");
  feCheckReturn(options.parse());

  GenericSolver *solver = nullptr;

  std::string solverTypeStr(solverType);
  if(solverTypeStr == "CHNS_RT")
    solver = new CHNS_RT::Solver();
  else if(solverTypeStr == "CHNS_Sloshing")
    solver = new CHNS_Sloshing::Solver();
  else if(solverTypeStr == "CHNS_Bubble")
    solver = new CHNS_Bubble::Solver();
  else if(solverTypeStr == "CHNS_YoungLaplace_VolumeAveraged")
    solver = new CHNS_YoungLaplace::VolumeAveragedSolver();
  else if(solverTypeStr == "CHNS_Jurin_MassAveraged")
    solver = new CHNS_Jurin::MassAveragedSolver();
  else if(solverTypeStr == "CHNS_Jurin_VolumeAveraged")
    solver = new CHNS_Jurin::VolumeAveragedSolver();
  else if(solverTypeStr == "CHNS_BubbleContact_VolumeAveraged")
    solver = new CHNS_BubbleContact::VolumeAveragedSolver();
  else if(solverTypeStr == "Cylinder")
    solver = new Cylinder::Solver();
  else if(solverTypeStr == "Analytic")
    solver = new Analytic::Solver();
  else {
    feErrorMsg(FE_STATUS_ERROR, "Unexpected solver type : %s", solverTypeStr.data());
    return -1;
  }

  solver->useAnalyticalSolution = false;

  feCheck(solver->parseAndCreateFunctions(argc, argv));
  solver->printParameters();

  // Use this block when creating a new reference solution
  // or when running a single simulation.
  const int nConv = 1;

  bool generateNewReferenceSolution = true;
  bool exportReferenceSolution = true;
  bool readReferenceSolution = false;
  const int nIntervalsReferenceSolution = 10;

  int nIntervals            = 10;
  int nVertices             = 1000;
  int nTimeStepsPerInterval = 50;

  const int intervalFactor  = 1;
  const int verticesFactor  = 2;
  const int timestepFactor  = 1;

  // // Use this block when running convergence tests
  // const int nConv = 6;

  // bool generateNewReferenceSolution = false;
  // bool exportReferenceSolution = false;

  // bool readReferenceSolution = true;
  // const int nIntervalsReferenceSolution = 15;

  // int nIntervals            = 15;
  // int nVertices             = 32000;
  // int nTimeStepsPerInterval = 50;

  // const int intervalFactor  = 1;
  // const int verticesFactor  = 2;
  // const int timestepFactor  = 1;

  std::ofstream errorFile;
  errorFile.open(solver->getWriteDir() + "st_error.txt");

  PostProcData postProcData;

  //
  // Convergence loop
  //
  for(int iConv = 0; iConv < nConv; ++iConv, nIntervals *= intervalFactor,
                                              nVertices *= verticesFactor,
                                  nTimeStepsPerInterval *= timestepFactor)
  {
    std::ofstream postProcFile;
    postProcFile.open(solver->getWriteDir() + "postProcessing_iConv" + std::to_string(iConv) + ".txt");

    AdaptData adapter;
    feCheck(adapter.parseAndInitialize(argc, argv,
                                       solver->getWriteDir(),
                                       solver->getInitialMesh(),
                                       solver->getReadDir(),
                                       solver->getOrderForAdaptation()
                                       // , nIntervals,
                                       // nVertices,
                                       // nTimeStepsPerInterval,
                                       // nIntervalsReferenceSolution
                                       ));

    feInfo("Integrating from t0 = %f to tEnd = %f in %d intervals of %d time steps each = %d time steps total",
      adapter.t0, adapter.tEnd, adapter.nIntervals, adapter.nTimeStepsPerInterval, adapter.nIntervals * adapter.nTimeStepsPerInterval);

    // The errors for various fields for each fixed-point
    // iteration and for each sub-interval
    std::vector<std::vector<std::vector<double>>> allErrors(adapter.nFixedPoint, std::vector<std::vector<double>>(adapter.nIntervals));
    int effectiveNst = 1;

    //
    // Create the reference test case for error estimation
    // Ideally, should be done only once but the adapter is created in the loop for now
    //
    if(readReferenceSolution && !solver->useAnalyticalSolution && !generateNewReferenceSolution) {
      feInfo("Reading reference solution on each interval");
      setVerbose(0);
      // #if defined(HAVE_OMP)
      // #pragma omp parallel for
      // #endif
      for(int iI = 0; iI < nIntervalsReferenceSolution; ++iI)
      {
        // std::string referenceMeshName = adapter.writeDir + "referenceMesh1em2_10steps.msh";
        // std::string solutionFileName = adapter.writeDir + "referenceSolution1em2_10steps.txt";
        std::string referenceMeshName = adapter.writeDir + "referenceMesh_interval" + std::to_string(iI) + ".msh";
        std::string solutionFileName = adapter.writeDir + "referenceSolution_interval" + std::to_string(iI) + ".txt";
        solver->readReferenceTestCase(iI, adapter, referenceMeshName, solutionFileName);
      }
      setVerbose(verbosity);
    }

    // Fixed point loop
    //
    for(int ifp = 0; ifp < adapter.nFixedPoint; ++ifp, ++adapter.iFixedPoint)
    {
      // Update solver parameters as the mesh is adapted
      solver->updateParameters(ifp);
      if(solver->isCHNS() && ifp > 0)
      {
        solver->printParameters();
        
        // Lp norm
        adapter.p = 100.;
        for(int iI = 0; iI < nIntervals; ++iI)
        {
          adapter.allOptions[iI].LpNorm = adapter.p;
        }
      }

      // Reset
      adapter.currentTime = adapter.t0;
      adapter.currentStep = 0;

      for(int iI = 0; iI < adapter.nIntervals; ++iI)
        for(auto &pair : adapter.allHi[iI])
          pair.second.setToZero();

      // Read adapted mesh if not at the first fixed-point iteration
      if(ifp > 0)
        feCheck(adapter.readAdaptedMeshes());

      //
      // Compute solution on each subinterval
      // Compute metric every few time steps
      //
      for(int iI = 0; iI < adapter.nIntervals; ++iI)
      {
        feInfo("");
        feInfo("////////////////////////////////////////////////");
        feInfo("Fixed-point %d/%d - Interval %d", ifp, adapter.nFixedPoint, iI);
        feInfo("////////////////////////////////////////////////");
        feInfo("");

        if(iI > 0) {
          //
          // Project solution on the mesh of the next sub-interval
          //
          feCheck(solver->projectSolution(adapter.allMeshes[iI-1], adapter.allMeshes[iI],
                                         adapter.allContainers[iI-1],
                                         adapter.bufferContainer,
                                         adapter.allContainers[iI],
                                         adapter.allNumberings[iI-1],
                                         adapter.allSpaces[iI-1],
                                         adapter.allEssentialSpaces[iI-1]));

          feInfo("");
          feInfo("DOF in source container : %d", adapter.allContainers[iI-1]->getNbDOFs());
          feInfo("DOF in buffer container : %d", adapter.bufferContainer.getNbDOFs());
          if(adapter.allContainers[iI])
            feInfo("DOF in target container : %d", adapter.allContainers[iI]->getNbDOFs());
          else
            feInfo("DOF in target container is null");
          feInfo("");
        }

        //
        // Solve on sub-interval and compute metrics
        //
        feMesh *mesh = adapter.allMeshes[iI];
        std::vector<feSpace*> &spaces = adapter.allSpaces[iI];
        std::vector<feSpace*> &essentialSpaces = adapter.allEssentialSpaces[iI];
        feCheck(solver->createSpaces(mesh, spaces, essentialSpaces));
        feMetaNumber *numbering = new feMetaNumber(mesh, spaces, essentialSpaces);
        if(ifp > 0)
          delete adapter.allNumberings[iI];
        adapter.allNumberings[iI] = numbering;
        feSolution sol(numbering->getNbDOFs(), spaces, essentialSpaces);

        // Start from projected solution if not first sub-interval
        if(iI > 0) {
          // sol.setSolFromContainer(adapter.allContainers[iI-1], 0); // << ============= SUS
          sol.setSolFromContainer(adapter.allContainers[iI], 0);
          for(auto *s : spaces)
            s->setDOFInitialization(dofInitialization::PREVIOUS_SOL);
          for(auto *s : essentialSpaces)
            s->setDOFInitialization(dofInitialization::NODEWISE);
        }

        feExporter *exporter = nullptr;
        std::string vtkFileRoot = adapter.writeDir + "sol_iConv" + std::to_string(iConv);
        if(adapter.exportSolution) {
          feCheck(createVisualizationExporter(exporter, VTK, numbering, &sol, mesh, spaces));
        }
        feExportData exportData = {exporter, adapter.exportFrequency, vtkFileRoot};

        feCheck(solver->solve(adapter, iI, mesh, &sol, numbering, spaces, numbering->getNbUnknowns(), exportData));

        // Compute error on this sub-interval
        if(generateNewReferenceSolution)
          allErrors[ifp][iI].resize(4);
        else {
          if(readReferenceSolution && nIntervalsReferenceSolution % nIntervals != 0) {
            feErrorMsg(FE_STATUS_ERROR, "Number of intervals is not a divider of the reference number of intervals. Probable mismatch.");
            delete solver;
            finalize();
            return -1;
          }

          const int refToCurrentRatio = nIntervalsReferenceSolution / nIntervals;
          const int iIntervalReference = ((iI + 1) * refToCurrentRatio) - 1;

          feInfo("Computing error:");
          if(readReferenceSolution)
          {
            feInfo("If reference solution is used, comparing:");
            feInfo("Reference solution at end of interval: %d (end time = %1.4e)",
              iIntervalReference, adapter.referenceTestCase[iIntervalReference].sol->getCurrentTime());
            feInfo("Current   solution at end of interval: %d (end time = %1.4e)",
              iI, sol.getCurrentTime());
          }

          feCheck(solver->computeError(iIntervalReference, adapter, spaces, mesh, &sol, allErrors[ifp][iI]));
        }

        //
        // Post-process
        //
        feCheck(solver->computePostProcessing(adapter, spaces, mesh, &sol, postProcData));
        for(auto val : postProcData.res) {
          postProcFile << val << "\t";
        }
        postProcFile << "\n";
        postProcFile.flush();

        // Write the isoline 0 if it is a CHNS problem
        if(solver->isCHNS()) {
          std::ofstream isolineFile;
          std::string name = solver->getWriteDir() + "isoline_iConv" + std::to_string(iConv) + "_ifp" + std::to_string(ifp) + "_endOfInterval" + std::to_string(iI) + ".txt";
          isolineFile.open(name);

          // Print time
          isolineFile << sol.getCurrentTime() << "\n";

          // Print points on the isoline
          for(size_t iPt = 0; iPt < postProcData.isoline.size() / 3; ++iPt) {
            isolineFile << postProcData.isoline[3 * iPt + 0] << "\t" << postProcData.isoline[3 * iPt + 1] << "\n";
          }
          isolineFile.close();
        }

        delete exporter;

        for(int jI = 0; jI < adapter.nIntervals; ++jI) {
          if(adapter.allContainers[jI] != nullptr)
            feInfo("Computed interval %d - Container %d has %d DOF", iI, jI, adapter.allContainers[jI]->getNbDOFs());
        }
      }

      if(adapter.adapt) {
        if(adapter.steadyMetric)
        {
          // Adapt w.r.t. solution of the last time step
          std::vector<int> numVerticesInAdaptedMeshes(adapter.nIntervals, 0);
          for(int iI = 0; iI < adapter.nIntervals; ++iI)
          {
            feSolution sol(*adapter.allContainers[iI], 0, adapter.allSpaces[iI], adapter.allEssentialSpaces[iI]);
            feCheck(solver->computeLastMetric(adapter, iI, adapter.allMeshes[iI], &sol, adapter.allNumberings[iI], adapter.allSpaces[iI]));
            setVerbose(0);
            feCheck(adapter.allMeshes[iI]->adapt(adapter.allMetrics[iI]));
            setVerbose(verbosity);

            // Read the adapted mesh to check the number of vertices
            feMesh2DP1 mesh(adapter.allMetrics[iI]->getOptions().adaptedMeshName);
            numVerticesInAdaptedMeshes[iI] = mesh.getNumVertices();
          }
          
          feInfo("Number of vertices in adapted meshes :");
          for(int iI = 0; iI < adapter.nIntervals; ++iI)
            feInfo("%d : %d", iI, numVerticesInAdaptedMeshes[iI]);
          int totalVertices = std::accumulate(numVerticesInAdaptedMeshes.begin(), numVerticesInAdaptedMeshes.end(), 0);

          effectiveNst = totalVertices;

          feInfo("Total  number of space vertices : %d", effectiveNst);
          feInfo("Target number of space vertices : %d", nVertices);

        } else
        {
          // Adapt on each sub-interval :
          //
          // Compute metric scaling
          //
          // Space-time complexity from Nv :
          // const double targetNst = (double) adapter.nVertices * adapter.nIntervals * adapter.nTimeStepsPerInterval;
          const double targetNst = (double) adapter.nVertices * adapter.nIntervals;
          const double dim = adapter.dim;
          const double p = adapter.p;
          const double nT = (double) adapter.nTimeStepsPerInterval;
          const double expDeterminantGlobal =   p / (2. * p + dim);
          const double expDeterminantLocal  = -1. / (2. * p + dim);

          double sumKj = 0.;
          for(int iI = 0; iI < adapter.nIntervals; ++iI)
          {
            feMetric *metricField = adapter.allMetrics[iI];
            metricField->setMetrics(adapter.allHi[iI]);
            const double Kj = metricField->computeIntegralOfDeterminant(expDeterminantGlobal);
            // Use simplification for constant time step. Otherwise, multiply Kj
            // by integral of (time step)^-1 on sub-interval before summing.
            sumKj += Kj;
          }
          
          // Integral of (time step)^-1 is nT, exponents of nT reduce to 2/d
          const double globalScalingFactor = pow(targetNst / (nT * sumKj), 2./dim);

          //
          // Apply scaling and adapt all sub-interval meshes
          //
          std::vector<int> numVerticesInAdaptedMeshes(adapter.nIntervals, 0);
          for(int iI = 0; iI < adapter.nIntervals; ++iI)
          {
          #if defined(HAVE_GMSH)
            feMetric *metricField = adapter.allMetrics[iI];
            metricField->scaleMetricsByDeterminant(globalScalingFactor, expDeterminantLocal);
            initializeGmsh();
            feMetricOptions my_options = metricField->getOptions();
            gmsh::open(my_options.backgroundMeshfile);
            // openGmshModel(my_options.backgroundMeshfile);
            metricField->setDebugPrefix(std::to_string(iI));
            metricField->applyGradation();
            metricField->writeMetricField();
            setVerbose(0);
            adapter.allMeshes[iI]->adapt(metricField);
            setVerbose(verbosity);

            // Read the adapted mesh to check the number of vertices
            feMesh2DP1 mesh(my_options.adaptedMeshName);
            numVerticesInAdaptedMeshes[iI] = mesh.getNumVertices();
          #else
            UNUSED(globalScalingFactor, expDeterminantLocal);
          #endif
          }
          
          feInfo("Number of vertices in adapted meshes :");
          for(int iI = 0; iI < adapter.nIntervals; ++iI)
            feInfo("%d : %d", iI, numVerticesInAdaptedMeshes[iI]);
          int totalVertices = std::accumulate(numVerticesInAdaptedMeshes.begin(), numVerticesInAdaptedMeshes.end(), 0);

          effectiveNst = totalVertices * adapter.nTimeStepsPerInterval;

          feInfo("Total  number of space      vertices : %d", totalVertices);
          feInfo("Total  number of space-time vertices : %d", effectiveNst);
          feInfo("Target number of space-time vertices : %f", targetNst);
        }
      } // if adapt
    } // fixed-point loop

    //
    // Export the reference solution on each sub-interval
    //
    if(exportReferenceSolution)
    {
      for(int iI = 0; iI < adapter.nIntervals; ++iI)
      {
        // std::string solutionFileName = adapter.writeDir + "referenceSolution.txt";
        std::string solutionFileName = adapter.writeDir + "currentSolution_interval" + std::to_string(iI) + ".txt";
        feSolution sol(*adapter.allContainers[iI], 0, adapter.allSpaces[iI], adapter.allEssentialSpaces[iI]);
        feInfo("Writing solution with time %f and %d (%d) dofs", sol.getCurrentTime(), sol.getNumDOFs(), adapter.allContainers[iI]->getNbDOFs());
        sol.printSol(solutionFileName);
      }

      // Check: read new mesh and solution
      feInfo("Check : Reading reference solution all each interval");
      for(int iI = 0; iI < adapter.nIntervals; ++iI)
      {
        std::string referenceMeshName = adapter.allMeshes[iI]->getID();
        std::string solutionFileName = adapter.writeDir + "currentSolution_interval" + std::to_string(iI) + ".txt";
        solver->readReferenceTestCase(iI, adapter, referenceMeshName, solutionFileName);
      }
    }

    //
    // Space-time error
    //
    for(int ifp = 0; ifp < adapter.nFixedPoint; ++ifp)
    {
      // Space-time error (integral) for each field and this fixed-point iteration
      std::vector<double> e_st(solver->getNumFields(), 0.);

      for(int iI = 0; iI < adapter.nIntervals; ++iI) {
        for(int j = 0; j < solver->getNumFields(); ++j) {

          // L1 error from the paper?
          e_st[j] += adapter.dt * adapter.nTimeStepsPerInterval * allErrors[ifp][iI][j];

          // Linf error in time
          // e_st[j] = fmax(e_st[j], allErrors[ifp][iI][j]);

          feInfo("Fixed-point iteration %d - Interval %2d - Error on field %d = %+-1.6e", ifp, iI, j, allErrors[ifp][iI][j]);
        }
      }

      for(int j = 0; j < solver->getNumFields(); ++j) {
        feInfo("Fixed-point iteration %d - Space-time error for field %d = %+-1.6e", ifp, j, e_st[j]);
      }
      feInfo("");

      // Print error of the last fixed-point iteration to file
      if(ifp == adapter.nFixedPoint-1) {
        const int targetNst = adapter.nIntervals * adapter.nVertices * adapter.nTimeStepsPerInterval;

        errorFile << adapter.nVertices << "\t";
        errorFile << (effectiveNst / adapter.nTimeStepsPerInterval) << "\t"; // sum of all effective vertices
        errorFile << adapter.nIntervals << "\t";
        errorFile << adapter.nTimeStepsPerInterval << "\t";

        errorFile << targetNst << "\t";
        errorFile << effectiveNst << "\t";

        errorFile << adapter.dt;
        for(int j = 0; j < solver->getNumFields(); ++j) {
          errorFile << "\t" << std::setprecision(8) << e_st[j];
        }
        errorFile << "\n";
        errorFile.flush();
      }
    }

    postProcFile.close();
    
  } // Convergence loop

  errorFile.close();

  delete solver;

  finalize();
  return 0;
}