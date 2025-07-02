#ifndef _GENERIC_SOLVER_
#define _GENERIC_SOLVER_

#include "feBoundaryConditions.h"
#include "feExporter.h"
#include "feFunction.h"
#include "feMessage.h"
#include "feNonLinearSolver.h"
#include "feNorm.h"
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

class SolverBase;
class TransientAdapter;

//
// A callback used in executable scripts to specify the operations
// to be performed during post-processing. This is because the
// types of operations varies greatly from one application to another,
// so it would be awkward to try and have a generic built-in function
// to handle all kinds of post-processings.
//
typedef std::function<feStatus(const SolverBase             *solver,
                               const std::vector<feSpace *> &spaces,
                               feSolution                   *sol,
                               PostProcData                 &data)>
  PostProcCallback;

/*

  A generic abstract solver, which can be inserted into an unsteady
  mesh adaptation loop.

*/
class SolverBase
{
public:
  Parameters::NonLinearSolver _NLSolver_parameters;
  Parameters::TimeIntegration _TimeIntegration_parameters;

  std::vector<FEDescriptor> _fieldDescriptors;
  std::vector<const BoundaryConditions::BoundaryCondition *>
    _boundaryConditions;

  PostProcCallback *_postprocCallback;

  // Number of unknown (solved) scalar and/or vector fields
  int _numFields = 0;

  // Scalar and vecotr exact solutions, along with the index of the
  // associated space given in fieldDescriptors.
  std::vector<std::pair<int, const feFunction *>> _scalarExactSolutions;
  std::vector<std::pair<int, const feVectorFunction *>> _vectorExactSolutions;

  std::vector<std::pair<int, normType>> _normTypes;

  mutable std::vector<std::vector<double>> _postProcessingData;

public:
  SolverBase(const Parameters::NonLinearSolver &NLSolver_parameters,
             const Parameters::TimeIntegration &TimeIntegration_parameters,
             const std::vector<FEDescriptor>   &fieldDescriptors,
             const std::vector<const BoundaryConditions::BoundaryCondition *>
                              &boundaryConditions,
             PostProcCallback *postProcCallback = nullptr)
    : _NLSolver_parameters(NLSolver_parameters)
    , _TimeIntegration_parameters(TimeIntegration_parameters)
    , _fieldDescriptors(fieldDescriptors)
    , _boundaryConditions(boundaryConditions)
    , _postprocCallback(postProcCallback)
  {}
  virtual ~SolverBase(){};

  int  getNumFields() const { return _numFields; };
  void setScalarSolutions(
    const std::vector<std::pair<int, const feFunction *>> &solutions)
  {
    _scalarExactSolutions = solutions;
  }
  void setVectorSolutions(
    const std::vector<std::pair<int, const feVectorFunction *>> &solutions)
  {
    _vectorExactSolutions = solutions;
  }
  void setNormTypes(
    const std::vector<std::pair<int, normType>> &normTypes)
  {
    _normTypes = normTypes;
  }
  std::vector<std::vector<double>> getPostProcessingData() const
  {
    return _postProcessingData;
  }

private:
  feStatus addBoundarySpace(const BoundaryConditions::BoundaryCondition *bc,
                            const FEDescriptor     &spaceDescriptor,
                            feMesh                 *mesh,
                            std::vector<feSpace *> &spaces,
                            std::vector<feSpace *> &essentialSpaces) const;

public:
  //
  // Create the finite element spaces (domain and boundaries)
  //
  feStatus createSpaces(feMesh                 *mesh,
                        std::vector<feSpace *> &spaces,
                        std::vector<feSpace *> &essentialSpaces) const;

  //
  // Create the vector of linear and bilinear forms to assemble
  //
  virtual feStatus
  createBilinearForms(const std::vector<feSpace *>  &spaces,
                      std::vector<feBilinearForm *> &forms) const = 0;

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

  //
  // When used within a TransientAdapter, perform any update to the solver
  // parameters before computing the next global fixed-point iteration.
  //
  virtual void updateBeforeFixedPointIteration(const int /*iter*/) {}

  //
  // Solve with time stepping information provided by the TimeIntegration
  // parameters. Used when using this solver as a "stand-alone" solver.
  //
  // The generic solve method does the following:
  //
  //  - Create the linear and bilinear forms
  //  - Create a linear system
  //  - Create a time integrator.
  //  - Integrate N time steps
  //  - Free the allocated memory
  //
  // This method can be overriden by derived solvers to perform
  // more specific actions.
  //
  virtual feStatus solve(feMesh                 *mesh,
                         feSolution             *sol,
                         feMetaNumber           *numbering,
                         std::vector<feSpace *> &spaces,
                         feExportData           &exportData) const;

  //
  // Solve with time stepping information provided by a TransientAdapter.
  // Used inside the unsteady adaptation loop when providing this solver
  // to a TransientAdapter.
  //
  // The generic solve method does the following:
  //
  //  - Create the linear and bilinear forms
  //  - Create a linear system
  //  - Create a time integrator. The time integrator is restarted
  //    from the solution on the previous time sub-interval, after
  //    this solution has been transfered to the mesh of this interval.
  //  - Integrate N time steps
  //  - Compute a metric field based on the prescribed field every
  //    prescribed number of time steps (for now 1).
  //  - Increment the time integral of the metric field
  //  - Store the solution at the last N BDF time steps into the
  //    container of the adapter.
  //  - Free the allocated memory
  //
  // This method can be overriden by derived solvers to perform
  // more specific actions.
  //
  virtual feStatus solve(feMesh                 *mesh,
                         feSolution             *sol,
                         feMetaNumber           *numbering,
                         std::vector<feSpace *> &spaces,
                         feExportData           &exportData,
                         TransientAdapter       &adapter,
                         const int               iInterval) const;

  //
  // Transfer solution from current mesh to next mesh
  //
  feStatus projectSolution(feMesh2DP1                  *currentMesh,
                           feMesh2DP1                  *nextMesh,
                           feSolutionContainer         *container,
                           feSolutionContainer         &bufferContainer,
                           feSolutionContainer        *&nextContainer,
                           feMetaNumber                *numbering,
                           const std::vector<feSpace *> spaces,
                           const std::vector<feSpace *> essentialSpaces) const;

  //
  // Compute errors based on either an exact solution
  // or a finer reference solution.
  //
  virtual feStatus computeError(const std::vector<feSpace *> &spaces,
                                feMesh                       *mesh,
                                feSolution                   *sol,
                                const int            iIntervalReferenceSolution,
                                TransientAdapter    &adapter,
                                std::vector<double> &errors) const = 0;

  //
  // Call the provided post-processing callback.
  //
  feStatus computePostProcessing(const std::vector<feSpace *> &spaces,
                                 feSolution                   *sol,
                                 PostProcData                 &data) const
  {
    // Simply call the given callback
    if (_postprocCallback)
    {
      feCheckReturn((*_postprocCallback)(this, spaces, sol, data));
    }
    return FE_STATUS_OK;
  };
};

#endif