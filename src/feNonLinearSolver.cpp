
#include "feNonLinearSolver.h"

extern int FE_VERBOSE;

#define NEWTON_OUTPUT_PRECISION "%10.6e"

//
// Decide if the Jacobian matrix should be recomputed
// using various heuristics
//
bool recomputeMatrixHeuristic(const feNLSolverOptions &tol,
                              const bool isStationary,
                              const int iter,
                              const double currentResidual,
                              const double previousResidual)
{
  // Recompute at each iteration when solving for steady-state solution
  if(isStationary)
    return true;

  // Transient problem
  // If at first or prescribed number of steps, recompute
  if(iter % tol.recomputeJacobianEveryNsteps == 0)
    return true;

  // Recompute if nonlinear residual did not decrease enough
  if(currentResidual > tol.residualDecrease * previousResidual)
    return true;

  return false;
}

feStatus solveNewtonRaphson(feLinearSystem *linearSystem,
                            feSolution *sol,
                            feSolutionContainer *container,
                            const feNLSolverOptions &tol,
                            const bool solveForTimeDerivative)
{
  feInfoCond(FE_VERBOSE > 0, "\t\tNONLINEAR SOLVER:");

  bool stop = false, isStationary = (container->getNbSol() == 1);
  int iter = 0, linearSystemIter;
  double normCorrection = 0., normResidual = 0., normAxb = 0., previousResidual = 1e22;

  // Always compute Jacobian matrix for first iteration
  linearSystem->setRecomputeStatus(true);

  if(solveForTimeDerivative) {
    sol->setC0(1.);
  }

  // Newton-Rapshon iteration
  while(!stop)
  {
    if(!solveForTimeDerivative) {
      // Update time derivative
      container->computeSolTimeDerivative(sol, linearSystem);
      // For clarity should be replaced by:
      // container->computeCurrentSolDot(linearSystem);
      // sol->setSolDotFromContainer(container, 0);
    }

    //
    // Reset, assemble and solve the linear system J(u) * du = -NL(u)
    //
    linearSystem->setToZero();

    // Check next residual norm and return if tolerance is reached
    linearSystem->assembleResiduals(sol);
    linearSystem->getRHSMaxNorm(&normResidual);
    if(iter > 0 && normResidual <= tol.tolResidual) {
      feInfoCond(FE_VERBOSE > 0,
               "\t\t\tStopping because next residual norm ||NL(u)|| = " NEWTON_OUTPUT_PRECISION " is"
               " below prescribed tolerance (%10.4e)", normResidual, tol.tolResidual);
      break;
    }

    // Assemble and constrain FE matrix
    if(linearSystem->getRecomputeStatus()) {
      linearSystem->assembleMatrices(sol, solveForTimeDerivative);
    }
    
    linearSystem->constrainEssentialComponents(sol);

    // Solve J(u) * du = -NL(u)
    bool successSolve = linearSystem->solve(&normCorrection, &normResidual, &normAxb, &linearSystemIter);

    if(!successSolve) {
      feWarning("It. %2d : "
                "||J*du - NL|| = " NEWTON_OUTPUT_PRECISION " (%4d iter.) \t "
                "||du|| = " NEWTON_OUTPUT_PRECISION " \t "
                "||NL(u)|| = " NEWTON_OUTPUT_PRECISION " (%s)",
                ++iter, normAxb, linearSystemIter, normCorrection, normResidual,
                linearSystem->getRecomputeStatus() ? "true" : "false");
      return feErrorMsg(FE_STATUS_ERROR, "Could not solve linear system at iter %2d )-:", iter);
    }

    // Abort if nonlinear term NL(u) is greater than divergence tolerance
    if(normResidual > tol.tolDivergence) {

      feErrorMsg(FE_STATUS_ERROR,
                      "Residual norm is greater than divergence tolerance: ||NL(u)|| = %10.4e > %10.4e",
                      normResidual, tol.tolDivergence);
      exit(-1);
      return feErrorMsg(FE_STATUS_ERROR,
                      "Residual norm is greater than divergence tolerance: ||NL(u)|| = %10.4e > %10.4e",
                      normResidual, tol.tolDivergence);
    }

    // Apply correction u += du and copy solution vector in container
    linearSystem->correctSolution(sol, solveForTimeDerivative);
    container->setCurrentSolution(sol);
    container->setCurrentSolutionDot(sol);

    // Status message for successful iteration
    feInfoCond(FE_VERBOSE > 0,
               "\t\t\tIt. %2d : "
               "||J*du - NL|| = " NEWTON_OUTPUT_PRECISION " (%4d iter.) \t "
               "||du|| = " NEWTON_OUTPUT_PRECISION " \t "
               "||NL(u)|| = " NEWTON_OUTPUT_PRECISION " (%s)",
               ++iter, normAxb, linearSystemIter, normCorrection, normResidual,
               linearSystem->getRecomputeStatus() ? "true" : "false");

    // Decide if the Jacobian matrix should be recomputed or not
    bool flag = recomputeMatrixHeuristic(tol, isStationary, iter, normResidual, previousResidual);
    linearSystem->setRecomputeStatus(flag);

    previousResidual = normResidual;

    stop = (normResidual <= tol.tolResidual) || (normCorrection <= tol.tolCorrection) || (iter > tol.maxIter);
  }

  if(normResidual <= tol.tolResidual)
  {
    // Solver converged as expected
    feInfoCond(
      FE_VERBOSE > 0,
      "\t\t\tConverged in %2d Newton iteration(s) - Residual converged: "
      "||du|| = " NEWTON_OUTPUT_PRECISION " \t "
      "||NL(u)|| = " NEWTON_OUTPUT_PRECISION "", iter,
      normCorrection, normResidual);
    return FE_STATUS_OK;
  } else if(normCorrection <= tol.tolCorrection)
  {
    // Increment du is low enough but residual NL(u) is not
    feInfoCond(
      FE_VERBOSE > 0,
      "\t\t\tConverged in %2d Newton iteration(s) (increment converged): ||du|| = %10.10e \t ||NL(u)|| = %10.10e", iter,
      normCorrection, normResidual);
    printf("\n");
    feWarning("Nonlinear solver converged because increment du = %10.4e is below prescribed tolerance (%1.4e),\n"
      "but the equation residual NL(u) = %10.4e is not (tol = %1.4e). Consider decreasing the tolerance for du.\n",
      normCorrection, tol.tolCorrection, normResidual, tol.tolResidual);
    return FE_STATUS_OK;
  } else
  {
    // Solver did not converge
    feErrorMsg(FE_STATUS_ERROR,
                      "Nonlinear solver did not converge at iter %2d : ||J*du - NL|| = %10.10e  "
                      "(%4d iter.) \t ||du|| = "
                      "%10.10e \t ||NL(u)|| = %10.10e",
                      iter, normAxb, linearSystemIter, normCorrection, normResidual);
    exit(-1);
    return feErrorMsg(FE_STATUS_ERROR,
                      "Nonlinear solver did not converge at iter %2d : ||J*du - NL|| = %10.10e  "
                      "(%4d iter.) \t ||du|| = "
                      "%10.10e \t ||NL(u)|| = %10.10e",
                      iter, normAxb, linearSystemIter, normCorrection, normResidual);
  }
}