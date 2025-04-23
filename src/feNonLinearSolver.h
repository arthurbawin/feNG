#ifndef _FENONLINEARSOLVER_
#define _FENONLINEARSOLVER_

#include "feMessage.h"
#include "feLinearSystem.h"
#include "feSolution.h"
#include "feSolutionContainer.h"

//
// Options for the nonlinear solver
//
typedef struct feNLSolverOptionsStruct
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
// If solveForTimeDerivative = true, instead solve the nonlinear problem,
// but the unknown is the time derivative of the solution dudt.
// This is done by 
//   - replacing the coefficient c0 (e.g. first coefficient of BDF sum)
//     by 1 in the "transient" mass matrices, turning them into regular mass matrices,
//   - assembling only the residuals and the "transient" mass matrices
//
// For instance, the Navier-Stokes momentum equation:
//
// dudt + (u_0 dot grad)u_0 - div(sigma_0) - rho*f = 0
//
// becomes the equation for w := dudt:
//
// w + (u_0 dot grad)u_0 - div(sigma_0) - rho*f = 0
//
// and is solved for w. This effectively determines dudt|_t=0 weakly.
// This is a linear problem in w, which converges in a single iteration.
//
feStatus solveNewtonRaphson(feLinearSystem *linearSystem,
                            feSolution *sol,
                            feSolutionContainer *solDot,
                            const feNLSolverOptions &tol,
                            const bool solveForTimeDerivative = false);



#endif