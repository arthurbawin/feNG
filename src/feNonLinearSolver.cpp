
#include "feNonLinearSolver.h"

extern int FE_VERBOSE;

// #define DEBUG_CHNS

//
// Decide if the Jacobian matrix should be recomputed
// using various heuristics
//
bool recomputeMatrixHeuristic(feNLSolverOptions &tol,
                              bool isStationary,
                              int iter,
                              double currentResidual,
                              double previousResidual)
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

extern int currentBDFStep;

//////////////////////////////////////////////////////
void checkForNan(const std::vector<double> &vec)
{
  for(const auto &val : vec) {
    if(std::isnan(val)) {
      feInfo("Nan alert in NL solver");
      exit(-1);
    }
  }
}

#if defined(DEBUG_CHNS)
static void writeSolution(const feSolution *sol, const std::string &fileName)
{
  const std::vector<double> &solreg = sol->getSolution();
  FILE *myfile = fopen(fileName.data(), "w");
  for(size_t i = 0; i < solreg.size(); ++i) {
    fprintf(myfile, "%+-10.20e\n", solreg[i]);
  }
  fclose(myfile);
} 

static void writeSolutionDot(const feSolution *sol, const std::string &fileName)
{
  const std::vector<double> &soldot = sol->getSolutionDot();
  FILE *myfile = fopen(fileName.data(), "w");
  for(size_t i = 0; i < soldot.size(); ++i) {
    fprintf(myfile, "%+-10.20e\n", soldot[i]);
  }
  fclose(myfile);
}
#endif
//////////////////////////////////////////////////////

feStatus solveNewtonRaphson(feLinearSystem *linearSystem,
                            feSolution *sol,
                            feSolutionContainer *container,
                            feNLSolverOptions &tol)
{
  feInfoCond(FE_VERBOSE > 0, "\t\t\tNONLINEAR SOLVER:");

  bool stop = false, isStationary = (container->getNbSol() == 1);
  int iter = 0, linearSystemIter;
  double normCorrection = 0., normResidual = 0., normAxb = 0., previousResidual = 1e22;

  // Always compute Jacobian matrix for first iteration
  linearSystem->setRecomputeStatus(true);

  // Newton-Rapshon iteration
  while(!stop)
  {
    // Update time derivative
    container->computeSolTimeDerivative(sol, linearSystem);

    ////////////////////////////////////////////////////////
    checkForNan(sol->getSolution());
    checkForNan(sol->getSolutionDot());

    #if defined(DEBUG_CHNS)
    if(iter == 0) {
      std::string solName = "ref/sol_beforeIter0_BDF_" + std::to_string(currentBDFStep) + ".txt";
      writeSolution(sol, solName);
      solName = "ref/soldot_beforeIter0_BDF_" + std::to_string(currentBDFStep) + ".txt";
      writeSolutionDot(sol, solName);
    }
    #endif
    ////////////////////////////////////////////////////////

    // For clarity should be replaced by:
    // container->computeCurrentSolDot(linearSystem);
    // sol->setSolDotFromContainer(container, 0);

    //
    // Reset, assemble and solve the linear system J(u) * du = -NL(u)
    //
    linearSystem->setToZero();

    // Check next residual norm and return if tolerance is reached
    linearSystem->assembleResiduals(sol);
    linearSystem->getRHSMaxNorm(&normResidual);
    if(iter > 0 && normResidual <= tol.tolResidual) {
      feInfoCond(FE_VERBOSE > 0,
               "\t\t\t\tStopping because residual norm ||NL(u)|| = %10.10e is"
               " below prescribed tolerance (%10.4e)", normResidual, tol.tolResidual);
      break;
    }

    // Assemble and constrain FE matrix
    if(linearSystem->getRecomputeStatus()) {
      linearSystem->assembleMatrices(sol);
    }
    
    linearSystem->constrainEssentialComponents(sol);

    ////////////////////////////////////////////////////
    #if defined(DEBUG_CHNS)
    if(iter == 0) {
      // Write RHS
      std::string rhsName = "ref/rhs_atIter0_BDF_" + std::to_string(currentBDFStep) + ".txt";
      linearSystem->writeRHS(rhsName);
      // Write matrix
      std::string matName = "ref/mat_atIter0_BDF_" + std::to_string(currentBDFStep) + ".txt";
      linearSystem->writeMatrix(matName);
    }
    #endif
    ////////////////////////////////////////////////////////

    // Solve J(u) * du = -NL(u)
    bool successSolve = linearSystem->solve(&normCorrection, &normResidual, &normAxb, &linearSystemIter);

    ////////////////////////////////////////////////////////
    #if defined(DEBUG_CHNS)
    if(iter == 0) {
      // Write correction (solution of linear system)
      std::string corName = "ref/cor_atIter0_BDF_" + std::to_string(currentBDFStep) + ".txt";
      linearSystem->writeResidual(corName);
    }
    #endif
    ////////////////////////////////////////////////////////

    if(!successSolve) {
      feWarning("Iter %2d : ||J*du - NL(u)|| = %10.10e (%4d iter.) \t ||du|| = %10.10e \t "
                "||NL(u)|| = %10.10e (%s)",
                ++iter, normAxb, linearSystemIter, normCorrection, normResidual,
                linearSystem->getRecomputeStatus() ? "true" : "false");
      return feErrorMsg(FE_STATUS_ERROR, "Could not solve linear system at iter %2d )-:", iter);
    }

    // Abort if nonlinear term NL(u) is greater than divergence tolerance
    if(normResidual > tol.tolDivergence) {

      ///////////////////////////////////////////////////////
      // std::string solName = "sol_beforeIter0_BDF_" + std::to_string(currentBDFStep) + ".txt";
      // writeSolution(sol, solName);
      // solName = "soldot_beforeIter0_BDF_" + std::to_string(currentBDFStep) + ".txt";
      // writeSolutionDot(sol, solName);
      // linearSystem->writeRHS("RHS.txt", container->getTime()[0]);
      // linearSystem->writeMatrix("matrix.txt", container->getTime()[0]);
      // linearSystem->writeResidual("cor.txt", container->getTime()[0]);
      ///////////////////////////////////////////////////////

      feErrorMsg(FE_STATUS_ERROR,
                      "Residual norm is greater than divergence tolerance: ||NL(u)|| = %10.4e > %10.4e",
                      normResidual, tol.tolDivergence);
      exit(-1);
      return feErrorMsg(FE_STATUS_ERROR,
                      "Residual norm is greater than divergence tolerance: ||NL(u)|| = %10.4e > %10.4e",
                      normResidual, tol.tolDivergence);
    }

    // Apply correction u += du and copy solution vector in container
    linearSystem->correctSolution(sol);
    container->setCurrentSolution(sol);
    container->setCurrentSolutionDot(sol);

    ////////////////////////////////////////////////////////
    #if defined(DEBUG_CHNS)
    if(iter == 0) {
      std::string solName = "ref/sol_afterIter0_BDF_" + std::to_string(currentBDFStep) + ".txt";
      writeSolution(sol, solName);
      solName = "ref/soldot_afterIter0_BDF_" + std::to_string(currentBDFStep) + ".txt";
      writeSolutionDot(sol, solName);
    }
    #endif
    ////////////////////////////////////////////////////////

    // Status message for successful iteration
    feInfoCond(FE_VERBOSE > 0,
               "\t\t\t\tIter %2d : ||J*du - NL(u)|| = %10.10e (%4d iter.) \t ||du|| = %10.20e \t "
               "||NL(u)|| = %10.20e (%s)",
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
      "\t\t\t\tConverged in %2d Newton iteration(s) (Residual converged):        ||du|| = %10.10e \t ||NL(u)|| = %10.10e", iter,
      normCorrection, normResidual);
    return FE_STATUS_OK;
  } else if(normCorrection <= tol.tolCorrection)
  {
    // Increment du is low enough but residual NL(u) is not
    feInfoCond(
      FE_VERBOSE > 0,
      "\t\t\t\tConverged in %2d Newton iteration(s) (Increment converged):       ||du|| = %10.10e \t ||NL(u)|| = %10.10e", iter,
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
                      "Nonlinear solver did not converge at iter %2d : ||J*du - NL(u)|| = %10.10e  "
                      "(%4d iter.) \t ||du|| = "
                      "%10.10e \t ||NL(u)|| = %10.10e",
                      iter, normAxb, linearSystemIter, normCorrection, normResidual);
    exit(-1);
    return feErrorMsg(FE_STATUS_ERROR,
                      "Nonlinear solver did not converge at iter %2d : ||J*du - NL(u)|| = %10.10e  "
                      "(%4d iter.) \t ||du|| = "
                      "%10.10e \t ||NL(u)|| = %10.10e",
                      iter, normAxb, linearSystemIter, normCorrection, normResidual);
  }
}