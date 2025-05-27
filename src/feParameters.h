#ifndef _FEPARAMETERS_
#define _FEPARAMETERS_

namespace Parameters
{
  struct TimeIntegration
  {
    double       dt;
    double       t_initial;
    double       t_final;
    unsigned int nTimeSteps;

    enum class TimeIntegrationMethod
    {
      stationary,
      bdf1,
      bdf2,
      dc2bdf1,
      dc3bdf1,
      dc3bdf2,
      eigensolver
    } method;

    enum class BDF2Starter
    {
      bdf1,
      exactSolution
    } bdf2starter = BDF2Starter::bdf1;
  };

  struct NonLinearSolver
  {
    // Nonlinear solver stops if the norm of the residual NL(u) passes below
    // this tolerance. This is the expected outcome if the nonlinear solver has
    // converged.
    double tolResidual = 1e-10;

    // Nonlinear solver stops if the norm of the solution correction du passes
    // below this tolerance (but issues a warning if the residual NL(u) has not
    // decreased enough).
    double tolCorrection = 1e-10;

    // Nonlinear solver stops if ||NL(u)|| > tolDivergence at any iteration and
    // returns an error.
    double tolDivergence = 1e4;

    // Nonlinear solver stops if numIter > maxIter and returns an error.
    double maxIter = 20;

    // For unsteady problems, Jacobian matrix is recomputed every N steps by
    // default. It is also recomputed if the nonlinear residual did not decrease
    // enough, see below.
    int recomputeJacobianEveryNsteps = 3;

    // For unsteady problems, Jacobian matrix is recomputed at next iteration if
    // the NL residual max norm did not decrease enough w.r.t. previous
    // iteration, that is, if:
    //
    //  ||NL(u)||_i > residualDecrease * ||NL(u)||_i-1
    //
    double residualDecrease = 1e-1;
  };

  struct FluidProperties
  {
    double density;
    double dynamic_viscosity;
  };

  struct Diffusion 
  {
    double placeholder; // Unused for now
  };

  /*
   * Parameters for the Cahn-Hilliard Navier-Stokes solver
   */
  struct CHNS
  {
    enum class Formulation
    {
      volumeAveraged,
      massAveraged
    } formulation;

    enum class MobilityType
    {
      constant,
      degenerate
    } mobilityType;

    bool limitPhaseMarker;

    FluidProperties fluids[2];

    double gravity = -9.81;
    double epsilon;
    double surfaceTension;
    double mobility;
    double mass_alpha;

    //
    // Control the CHNS parameters to modify
    // during mesh adaptation.
    //
    struct MeshAdaptation
    {
      // Decrease interface thickness every N global fixed point iterations
      int epsilonDecreaseFrequency = 2;

      // Divide epsilon by this value when condition is met
      double epsilonDecrease = 1.;
    } adaptation;
  };

  struct MeshAdaptation
  {
    bool adapt = true;
    int orderForAdaptation;
    int spaceIDForAdaptation = 0;
    int targetVertices = 5000;
    double Lp_norm = 2.;
    double dim = 2.;

    bool enableGradation = false;
    bool enableSpaceTimeGradation = false;
    double gradation = 2.;

    // Use the exact function derivatives to compute the metric field
    bool useExactDerivatives = false;

    // For unsteady adaptation
    int nFixedPoint;
    int nIntervals;
    int nTimeStepsPerInterval;
    int nIntervalsReferenceSolution;

    // The Lp norm used to compute the space-time error,
    // at the end of an unsteady fixed-point iteration
    enum class SpaceTimeErrorNorm
    {
      L1,
      Linf
    } spaceTimeNorm = SpaceTimeErrorNorm::L1;

    // Increase the number of target vertices every N fixed point iterations
    int verticesIncreaseFrequency = 1;

    // Multiply target number of vertices by this value when condition is met
    double targetVerticesIncrease = 1.;

    // Restart from prescribed meshes.
    // If not the provided initial mesh is used as starting mesh 
    // for all time sub-intervals.
    bool restartFromAdaptedMeshes = false;
    // If restarting from specified meshes, look for the mesh files:
    //
    // IO.writeDir + restartMeshesRoot + "_intervalXX.msh"
    std::string restartMeshesRoot = "";

    // When reading a reference mesh and solution, we look for the files:
    //
    // Mesh: IO.writeDir + referenceMeshRoot + "_intervalXX.msh"
    //
    // and
    // 
    // Solution file: IO.writeDir + referenceSolutionRoot + "_intervalXX.txt"
    bool readReferenceSolution = false;
    std::string referenceMeshRoot = "referenceMesh";
    std::string referenceSolutionRoot = "referenceSolution";
  };

  struct IO
  {
    std::string readDir = "";
    std::string writeDir = "";
    std::string initialMesh = "";

    bool exportVisualizationFile = true;
    int exportFrequency = 1;
  };
} // namespace Parameters

#endif