#ifndef _TRANSIENT_ADAPTER_
#define _TRANSIENT_ADAPTER_

#include "GenericSolver.h"
#include "feMesh.h"
#include "feMetric.h"
#include "feNumber.h"
#include "feParameters.h"
#include "feSolution.h"
#include "feSpace.h"

/*

  This class handles the unsteady mesh adaptation loop

*/

class TransientAdapter
{
public:
  bool wasAllocated = false;

  Parameters::MeshAdaptation  _adapt_parameters;
  Parameters::TimeIntegration _time_parameters;
  Parameters::IO              _io_parameters;

  std::ofstream postProcFile;

  int iFixedPoint = 0;

  double currentTime;
  int    currentTimeStep;

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
  // Size is nFixedPoint x nIntervals x nUnknownFields
  //
  std::vector<std::vector<std::vector<double>>> allErrors;

  // Similar as allErrors, but this is the error computed at each
  // time step. Size is nFixedPoint x (nUnknownFields + 1) x nTotalTimeSteps.
  //
  // The first field is associated with the time.
  std::vector<std::vector<std::vector<double>>> allErrorsAllTimeSteps;

  // The effective space-only and space-time complexity computed after all
  // meshes are adapted, for each fixed-point iteration.
  std::vector<int> spaceComplexity;
  std::vector<int> spaceTimeComplexity;
  std::vector<int> prescribedSpaceTimeComplexity;

  struct referenceTestCaseStruct
  {
    feMesh                *mesh;
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

  std::vector<int> getSpaceComplexities() const
  {
    return spaceComplexity;
  }
  std::vector<int> getSpaceTimeComplexities() const
  {
    return spaceTimeComplexity;
  }
  std::vector<int> getPrescribedSpaceTimeComplexities() const
  {
    return prescribedSpaceTimeComplexity;
  }
  auto getErrorAllTimeSteps() const
  {
    return allErrorsAllTimeSteps;
  }

  //
  // Set the callback "derivatives" as the exact derivatives of order "order"
  // for the computation of all metric fields.
  //
  feStatus setExactDerivatives(const int               order,
                               const feVectorFunction *derivatives);

  //
  // Allocate the structures and perform some compatibility checks between
  // parameters, and between the solver and the adapter.
  //
  feStatus allocate(const SolverBase &solver);

  feStatus readAdaptedMeshes();

  //
  // Any update in the transient adapter that should be done
  // before computing the next global fixed point iteration,
  // such as changes in the Lp norm.
  //
  void updateBeforeFixedPointIteration();

  //
  // Compute the space-time metric fields for all meshes
  // Adapt the mesh of each sub-interval after metrics were computed
  //
  feStatus computeSpaceTimeMetrics();

  //
  // Perform a single fixed point iteration:
  //
  // - (Optional updates in adapter and/or solver)
  // - Solve the problem on all time sub-intervals:
  //
  //   For each interval [ti, ti+1]:
  //
  //   - Transfer solution from previous interval (if > 0)
  //   - Advance solution to ti+1
  //   - Compute metrics every few time steps and increment
  //     time integral of metric
  //   - Compute error at end of interval
  //   - Post-process solution at end of interval
  //
  // - Compute all space-time metric fields
  // - Adapt all meshes (one per sub-interval)
  //
  feStatus fixedPointIteration(SolverBase &solver);

  //
  // Compute the space-time error for all fixed-point iterations and all fields:
  //
  //         /T  /
  // e_st =  |   |      |u_ref(x,t) - u(x,t)| dx dt  (L1 norm in time),
  //         /0  /Omega
  //
  // or
  //                   /
  // e_st =  max_t_i   |      |u_ref(x,t) - u(x,t)| dx (Linf norm).
  //                   /Omega
  //
  // Returns a vector of size nFixedPoint x nFields.
  //
  feStatus
  computeSpaceTimeErrors(std::vector<std::vector<double>> &spaceTimeErrors);

  feStatus finalize();
};

#endif