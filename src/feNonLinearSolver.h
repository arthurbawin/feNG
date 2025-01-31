#ifndef _FENONLINEARSOLVER_
#define _FENONLINEARSOLVER_

#include "feMessage.h"
#include "feLinearSystem.h"
#include "feSolution.h"
#include "feSolutionContainer.h"

//
// Options for the nonlinear solver
//
typedef struct
{
  // Nonlinear solver stops if the norm of the residual NL(u) passes below this tolerance.
  // This is the expected outcome if the nonlinear solver has converged.
  double tolResidual = 1e-10;

  // Nonlinear solver stops if the norm of the solution correction du passes below this
  // tolerance (but issues a warning if the residual NL(u) has not decreased enough).
  double tolCorrection = 1e-10;

  // Nonlinear solver stops if ||NL(u)|| > tolDivergence at any iteration and returns an error.
  double tolDivergence = 1e4;

  // // Nonlinear solver stops if numIter > maxIter and returns an error.
  double maxIter = 20;

  // For unsteady problems, Jacobian matrix is recomputed every N steps by default.
  // It is also recomputed if the nonlinear residual did not decrease enough, see below.
  int recomputeJacobianEveryNsteps = 3;

  // For unsteady problems, Jacobian matrix is recomputed at next iteration if
  // the NL residual max norm did not decrease enough w.r.t. previous iteration,
  // that is, if:
  //
  //  ||NL(u)||_i > residualDecrease * ||NL(u)||_i-1
  //
  double residualDecrease = 1e-1;

} feNLSolverOptions;

//
// Assembles and solves the nonlinear problem NL(u) = 0 at a given time step
// using the Newton-Raphson method.
//
feStatus solveNewtonRaphson(feLinearSystem *linearSystem,
                            feSolution *sol,
                            feSolutionContainer *solDot,
                            feNLSolverOptions &tol);



#endif