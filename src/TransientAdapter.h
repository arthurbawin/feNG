#ifndef _TRANSIENT_ADAPTER_
#define _TRANSIENT_ADAPTER_

#include "feMesh.h"
#include "feMetric.h"
#include "feNumber.h"
#include "feParameters.h"
#include "feSolution.h"
#include "feSpace.h"
#include "GenericSolver.h"

/*

  This class handles the unsteady mesh adaptation loop

*/

class TransientAdapter
{
public:
  Parameters::MeshAdaptation  _adapt_parameters;
  Parameters::TimeIntegration _time_parameters;
  Parameters::IO              _io_parameters;

  bool adapt           = true;
  bool exportSolution  = true;
  int  exportFrequency = 1;

  std::ofstream postProcFile;

  int iFixedPoint = 0;

  double currentTime;
  int currentTimeStep;

  // Time integration
  // double t0                    = 0.;
  // double tEnd                  = 0.1;
  // double currentTime           = t0;
  // int    currentStep           = 0;
  // int    nTimeStepsPerInterval = 4;
  // double dt = (tEnd - t0) / (nIntervals * nTimeStepsPerInterval);

  std::vector<feMesh2DP1 *>                allMeshes;
  std::vector<feMetric *>                  allMetrics;
  std::vector<feMetricOptions>             allOptions;
  std::vector<std::map<int, MetricTensor>> allHi;
  std::vector<feSolutionContainer *>       allContainers;
  std::vector<feMetaNumber *>              allNumberings;
  std::vector<std::vector<feSpace *>>      allSpaces;
  std::vector<std::vector<feSpace *>>      allEssentialSpaces;

  //
  // Error computed for each global fixed-point iteration, 
  // each sub-interval and for each field.
  //
  std::vector<std::vector<std::vector<double>>> allErrors;

  struct referenceTestCaseStruct
  {
    feMesh2DP1            *mesh;
    feMetaNumber          *numbering;
    feSolution            *sol;
    std::vector<feSpace *> spaces;
  };

  // A fine mesh solution for each interval, against which the error on coarser
  // meshes is evaluated
  std::vector<referenceTestCaseStruct> referenceTestCases;

  // A buffer container to use when transfering solution from one mesh to the
  // next
  feSolutionContainer bufferContainer;

public:
  TransientAdapter(const Parameters::MeshAdaptation  &adaptParameters,
                   const Parameters::TimeIntegration &timeParameters,
                   const Parameters::IO              &IOParameters);
  ~TransientAdapter()
  {
    for (auto *ptr : allMeshes)
      delete ptr;
    for (auto *ptr : allMetrics)
      delete ptr;
    for (auto *ptr : allContainers)
      delete ptr;
    for (auto *ptr : allNumberings)
      delete ptr;
    for (auto &vec : allSpaces)
      for (auto *ptr : vec)
        delete ptr;

    for (auto &refCase : referenceTestCases)
    {
      delete refCase.mesh;
      for (auto &s : refCase.spaces)
        delete s;
      delete refCase.numbering;
      delete refCase.sol;
    }
  }

  feStatus allocate(const SolverBase &solver);
  feStatus readAdaptedMeshes();

  //
  // Any updatein the transient adapter that should be done
  // before computing the next global fixed point iteration,
  // such as changes in the Lp norm.
  //
  feStatus updateBeforeFixedPointIteration();

  //
  // Compute the space-time metric fields for all meshes
  //
  feStatus computeSpaceTimeMetrics();

  //
  // Adapt the mesh of each sub-interval after metrics were computed
  //
  feStatus adaptAllMeshes();

  //
  // Perform a fixed point iteration:
  //
  // - (Optional updates in adapter and/or solver)
  // - Solution transfer from previous iteration
  // - Solver for the prescribed time interval
  // - Postprocessing
  // - Mesh adaptation
  //
  feStatus fixedPointIteration(SolverBase &solver);

  feStatus finalize();
};

#endif