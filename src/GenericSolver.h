#ifndef _GENERIC_SOLVER_
#define _GENERIC_SOLVER_

#include "feBoundaryConditions.h"
#include "feExporter.h"
#include "feFunction.h"
#include "feMessage.h"
#include "feNonLinearSolver.h"
#include "feNumber.h"
#include "feParameters.h"

//
// Generic class to store quantities computed during post-processing
//
class PostProcData
{
public:
  // A vector of post-processing results, which can be anything. Even a boat.
  std::vector<double> res;

  // The 0 isoline for CHNS test cases
  std::vector<double> isoline;

public:
  PostProcData(){};
};

//
// A callback used in executable scripts to specify the operations
// to be performed during post-processing. This is because the
// types of operations varies greatly from one application to another,
// so it would be awkward to try and have a generic built-in function
// to handle all kinds of post-processings.
//
typedef std::function<feStatus(const std::vector<feSpace *> &spaces,
                               feSolution                   *sol,
                               PostProcData                 &data)>
  PostProcCallback;

class TransientAdapter;

/*

  A generic abstract solver, which can be inserted into an unsteady
  mesh adaptation loop.

*/
class SolverBase
{
public:
  Parameters::NonLinearSolver _NLSolver_parameters;
  Parameters::TimeIntegration _TimeIntegration_parameters;

  std::vector<FEDescriptor>      _fieldDescriptors;
  std::vector<BoundaryCondition> _boundaryConditions;

  PostProcCallback *_postprocCallback;

public:
  SolverBase(const Parameters::NonLinearSolver    &NLSolver_parameters,
             const Parameters::TimeIntegration    &TimeIntegration_parameters,
             const std::vector<FEDescriptor>      &fieldDescriptors,
             const std::vector<BoundaryCondition> &boundaryConditions,
             PostProcCallback                     *postProcCallback = nullptr)
    : _NLSolver_parameters(NLSolver_parameters)
    , _TimeIntegration_parameters(TimeIntegration_parameters)
    , _fieldDescriptors(fieldDescriptors)
    , _boundaryConditions(boundaryConditions)
    , _postprocCallback(postProcCallback)
  {}
  virtual ~SolverBase(){};

  feStatus createSpaces(feMesh                 *mesh,
                        std::vector<feSpace *> &spaces,
                        std::vector<feSpace *> &essentialSpaces) const;

  //
  // Solve with time stepping information provided by the TimeIntegration
  // parameters. Used when using this solver as a "stand-alone" solver.
  //
  virtual feStatus solve(feMesh                 *mesh,
                         feSolution             *sol,
                         feMetaNumber           *numbering,
                         std::vector<feSpace *> &spaces,
                         feExportData           &exportData) const = 0;

  //
  // Solve with time stepping information provided by a TransientAdapter.
  // Used inside the unsteady adaptation loop when providing this solver
  // to a TransientAdapter.
  //
  virtual feStatus solve(feMesh                 *mesh,
                         feSolution             *sol,
                         feMetaNumber           *numbering,
                         std::vector<feSpace *> &spaces,
                         feExportData           &exportData,
                         TransientAdapter       &adapter,
                         const int               iInterval) const = 0;

  feStatus projectSolution(feMesh2DP1                  *currentMesh,
                           feMesh2DP1                  *nextMesh,
                           feSolutionContainer         *container,
                           feSolutionContainer         &bufferContainer,
                           feSolutionContainer        *&nextContainer,
                           feMetaNumber                *numbering,
                           const std::vector<feSpace *> spaces,
                           const std::vector<feSpace *> essentialSpaces) const;

  //
  // Create the structures for a reference test case
  // used as an ersatz of the exact solution.
  //
  feStatus readReferenceTestCase(const std::string       meshName,
                                 const std::string       solutionFileName,
                                 feMesh                *&mesh,
                                 std::vector<feSpace *> &spaces,
                                 feMetaNumber          *&numbering,
                                 feSolution            *&sol) const;

  virtual feStatus updateBeforeFixedPointIteration(const int iter) = 0;

  virtual feStatus computeError(const std::vector<feSpace *> &spaces,
                                feMesh                       *mesh,
                                feSolution                   *sol,
                                const int            iIntervalReferenceSolution,
                                TransientAdapter    &adapter,
                                std::vector<double> &errors) const = 0;

  feStatus computePostProcessing(const std::vector<feSpace *> &spaces,
                                 feSolution                   *sol,
                                 PostProcData                 &data) const
  {
    // Simply call the given callback
    if (_postprocCallback)
    {
      feCheckReturn((*_postprocCallback)(spaces, sol, data));
    }
    return FE_STATUS_OK;
  };
};

#endif