#include "feLinearSystem.h"

extern int FE_VERBOSE;

#if defined(HAVE_MKL)
#include "mkl.h"

feLinearSystemMklPardiso::feLinearSystemMklPardiso(std::vector<feBilinearForm *> bilinearForms,
                                                   int numUnknowns)
  : feLinearSystem(bilinearForms)
{
  recomputeMatrix = true;

  matrixOrder = numUnknowns;

  _EZCRS = new feEZCompressedRowStorage(numUnknowns, _formMatrices, _numMatrixForms);

  nz = _EZCRS->getNumNNZ();
  Ap = new PardisoInt[matrixOrder + 1];
  Aj = new PardisoInt[nz];
  _EZCRS->get_ia_Pardiso(&Ap);
  _EZCRS->get_ja_Pardiso(&Aj);
  Ax = _EZCRS->allocateMatrixArray();
  _EZCRS->setMatrixArrayToZero(Ax);

  du = new double[matrixOrder];
  residu = new double[matrixOrder];
  symbolicFactorization = true;
  //=================================================================
  //  INITIALISATION - MATRICE NON SYMETRIQUE
  //=================================================================
  // int num_procs = (int) nbthreads;
  // mkl_set_num_threads(num_procs);
  IPIVOT = 13;

  MTYPE = 11;
  SOLVER = 0; // sparse direct solver

  // pardisoinit (PT,&MYTPE,&SOLVER,IPARM,DPARM,&ERROR);
  // VERSION MKL PARDISO
  for(feInt i = 0; i < 64; i++) IPARM[i] = 0;
  IPARM[0] = 1; /* No solver default */
  IPARM[1] = 3; /* Fill-in reordering: (0)Min Degree, (2)METIS, (3)OpenMP */
  IPARM[2] = 0; /* Reserved. Set to zero for MKL version. Number of threads otherwise */
  IPARM[3] = 0; /* No iterative-direct algorithm */
  IPARM[4] = 0; /* No user fill-in reducing permutation */
  IPARM[5] = 0; /* Write solution into (0)x, (1)b */
  IPARM[6] = 0; /* Not in use - Number of iterative refinement steps performed */
  IPARM[7] = 2; /* Max numbers of iterative refinement steps. (2)=2 iterations */
  IPARM[8] = 0; /* Reserved */
  IPARM[9] = IPIVOT; /* Perturb the pivot elements with 1E-13 */ //-AG
  IPARM[10] = 0; /* (0)Disable scaling, (1) Use nonsymmetric permutation and scaling MPS */ //--AG
  IPARM[11] = 0; /* Not in use - Solve AX=B */
  IPARM[12] =
    1; /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
  IPARM[13] = 0; /* Output: Number of perturbed pivots */
  IPARM[14] = 0; /* Not in use - Peak memory on symbolic factorization*/
  IPARM[15] = 0; /* Not in use - Permanent memory on symbolic factorization */
  IPARM[16] = 0; /* Not in use - Peak memory on numerical factorization */
  IPARM[17] = -1; /* Output: Number of nonzeros in the factor LU */
  IPARM[18] = 0; /* Output: Mflops for LU factorization */
  IPARM[19] = 0; /* Output: Numbers of CG Iterations */
  IPARM[23] = 1; /* Parallel factorization control - (0) sequential - (1) parallel*/
  IPARM[24] = 0; /* Parallel LU solve control - (0) parallel - (1) sequential */

  IPARM[34] = 1; /* 0-based indexing - (0) Fortran-style (1-based) - (1) C-style (0-based) */

  for(feInt i = 0; i < 64; i++) PT[i] = NULL;

  N = (feInt)matrixOrder;
  NRHS = 1;
  // IPARM[2]= num_procs; // nombre de processeurs
  MAXFCT = 1;
  MNUM = 1;
  MSGLVL = 0;
  ERROR = 0;
}

double vectorMaxNorm(feInt N, double *V)
{
  double t = 0.0;
  for(feInt i = 0; i < N; i++) {
    t = fmax(t, fabs(V[i]));
  }
  return t;
}

// Print the matrix with a layout similar to PETSc's
void feLinearSystemMklPardiso::viewMatrix()
{
  // for(feInt i = 0; i < nz; i++) printf("%ld %g \n", i, Ax[i]);
  for(feInt i = 0; i < matrixOrder; i++) {
    printf("Row %ld: ", i);
    int debut = Ap[i];
    int fin = Ap[i + 1];
    for(feInt j = 0; j < fin - debut; ++j) {
      printf("(%lld, %f)  ", Aj[debut + j], Ax[debut + j]);
    }
    printf("\n");
  }
}

void feLinearSystemMklPardiso::assembleMatrices(feSolution *sol)
{
  tic();
  for(feInt eq = 0; eq < _numMatrixForms; ++eq) {
    feBilinearForm *f = _formMatrices[eq];
    feCncGeo *cnc = f->getCncGeo();
    int numColors = cnc->getNbColor();
    const std::vector<int> &nbElmPerColor = cnc->getNbElmPerColor();
    const std::vector<std::vector<int> > &listElmPerColor = cnc->getListElmPerColor();
    int numElementsInColor;
    std::vector<int> listElmC;

    for(int iColor = 0; iColor < numColors; ++iColor) {
      numElementsInColor = nbElmPerColor[iColor];
      listElmC = listElmPerColor[iColor];

      int elm;
      double **Ae;
      feInt sizeI;
      feInt sizeJ;
      std::vector<feInt> niElm;
      std::vector<feInt> njElm;
      std::vector<feInt> adrI;
      std::vector<feInt> adrJ;

#if defined(HAVE_OMP)
#pragma omp parallel for private(elm, f, niElm, njElm, sizeI, sizeJ, adrI, adrJ, Ae)               \
  schedule(dynamic)
#endif
      for(int iElm = 0; iElm < numElementsInColor; ++iElm) {
#if defined(HAVE_OMP)
        int eqt = eq + omp_get_thread_num() * _numMatrixForms;
        f = _formMatrices[eqt];
#endif
        elm = listElmC[iElm];

        // Compute element-wise matrix
        f->computeMatrix(sol, elm);
        Ae = f->getAe();

        // Determine global assignment indices
        adrI = f->getAdrI();
        adrJ = f->getAdrJ();
        sizeI = adrI.size();
        niElm.reserve(sizeI);
        for(feInt i = 0; i < sizeI; ++i) {
          if(adrI[i] < matrixOrder) niElm.push_back(i);
        }
        sizeJ = adrJ.size();
        njElm.reserve(sizeJ);
        for(feInt i = 0; i < sizeJ; ++i) {
          if(adrJ[i] < matrixOrder) njElm.push_back(i);
        }

        adrI.erase(std::remove_if(adrI.begin(), adrI.end(),
                                  [this](const int &x) { return x >= this->matrixOrder; }),
                   adrI.end());
        adrJ.erase(std::remove_if(adrJ.begin(), adrJ.end(),
                                  [this](const int &x) { return x >= this->matrixOrder; }),
                   adrJ.end());

        // Flatten Ae at relevant indices
        sizeI = adrI.size();
        sizeJ = adrJ.size();

        for(feInt i = 0; i < sizeI; i++) {
          feInt I = adrI[i];
          feInt debut = Ap[I];
          feInt fin = Ap[I + 1];
          feInt numColumns = fin - debut;

          // For each entry of the local matrix,
          // find the matching column in the sparse matrix.
          // The entries of the local matrix are not sorted for P2+ elements.
          for(feInt j = 0; j < sizeJ; ++j) {
            for(feInt J = 0; J < numColumns; ++J) {
              if(Aj[debut + J] == adrJ[j]) {
                Ax[debut + J] += Ae[niElm[i]][njElm[j]];
              }
            }
          }
          niElm.clear();
          njElm.clear();
        }
      }
    }
  }

  if(_displayMatrixInConsole) viewMatrix();

  feInfoCond(FE_VERBOSE > 1, "\t\t\t\tAssembled jacobian matrix in %f s", toc());
}

void feLinearSystemMklPardiso::assembleResiduals(feSolution *sol)
{
  tic();
  for(feInt eq = 0; eq < _numResidualForms; eq++) {
    feBilinearForm *f = _formResiduals[eq];
    feCncGeo *cnc = f->getCncGeo();
    int numColors = cnc->getNbColor();
    const std::vector<int> &numElemPerColor = cnc->getNbElmPerColor();
    const std::vector<std::vector<int> > &listElmPerColor = cnc->getListElmPerColor();
    std::vector<int> listElmC;

    for(int iColor = 0; iColor < numColors; ++iColor) {
      int numElementsInColor = numElemPerColor[iColor];
      listElmC = listElmPerColor[iColor];

      int elm;
      double *Be;
      feInt sizeI;
      std::vector<feInt> niElm;
      std::vector<feInt> adrI;

#if defined(HAVE_OMP)
#pragma omp parallel for private(elm, niElm, Be, f, adrI, sizeI) schedule(dynamic)
#endif
      for(int iElm = 0; iElm < numElementsInColor; ++iElm) {
#if defined(HAVE_OMP)
        int eqt = eq + omp_get_thread_num() * _numResidualForms;
        f = _formResiduals[eqt];
#endif
        elm = listElmC[iElm];

        // Compute the element-wise residual
        f->computeResidual(sol, elm);
        Be = f->getBe();
        adrI = f->getAdrI();
        sizeI = adrI.size();
        for(feInt i = 0; i < sizeI; i++) {
          if(adrI[i] < matrixOrder) residu[adrI[i]] += Be[i];
        }
      }
    }
  }
  feInfoCond(FE_VERBOSE > 1, "\t\t\t\tAssembled global residual in %f s", toc());
  // for(int i = 0; i < matrixOrder; ++i){
  //   printf("%g \n",residu[i]);
  // }
}

void feLinearSystemMklPardiso::assignResidualToDCResidual(feSolutionContainer *solContainer)
{
  for(feInt i = 0; i < matrixOrder; ++i) solContainer->_fResidual[0][i] = residu[i];
}

void feLinearSystemMklPardiso::applyCorrectionToResidual(double coeff, std::vector<double> &d)
{
  for(int i = 0; i < matrixOrder; ++i) residu[i] += coeff * d[i];
}

void feLinearSystemMklPardiso::correctSolution(feSolution *sol)
{
  // Est-ce efficace?
  // Pourquoi ne pas avoir un pointeur sur le vecteur solution;
  // Conversion de feInt à int?
  for(feInt i = 0; i < matrixOrder; i++) sol->incrementSolAtDOF((int)i, du[i]);
}

void feLinearSystemMklPardiso::correctSolution(double *sol)
{
  for(feInt i = 0; i < matrixOrder; i++) sol[i] += du[i];
}

bool checkPardisoErrorCode(PardisoInt &errorCode, std::string &pardisoStep)
{
  switch(errorCode) {
    case 0:
      // No error
      return true;
    case -1:
      feErrorMsg(FE_STATUS_ERROR,
                 "PARDISO %s step failed with error message: input inconsistent (-1)",
                 pardisoStep.data());
      return false;
    case -2:
      feErrorMsg(FE_STATUS_ERROR,
                 "PARDISO %s step failed with error message: not enough memory (-2)",
                 pardisoStep.data());
      return false;
    case -3:
      feErrorMsg(FE_STATUS_ERROR,
                 "PARDISO %s step failed with error message: reordering problem (-3)",
                 pardisoStep.data());
      return false;
    case -4:
      feErrorMsg(
        FE_STATUS_ERROR,
        "PARDISO %s step failed with error message: zero pivot, numerical factorization or "
        "iterative refinement problem. If the error appears during the solution phase, "
        "try to change the pivoting perturbation (iparm[9]) and also increase the number "
        "of iterative refinement steps. If it does not help, consider changing the "
        "scaling, matching and pivoting options (iparm[10], iparm[12], iparm[20]) (-4)",
        pardisoStep.data());
      return false;
    case -5:
      feErrorMsg(FE_STATUS_ERROR,
                 "PARDISO %s step failed with error message: unclassified (internal) error (-5)",
                 pardisoStep.data());
      return false;
    case -6:
      feErrorMsg(FE_STATUS_ERROR,
                 "PARDISO %s step failed with error message: reordering failed "
                 "(matrix types 11 and 13 only) (-6)",
                 pardisoStep.data());
      return false;
    case -7:
      feErrorMsg(FE_STATUS_ERROR,
                 "PARDISO %s step failed with error message: diagonal matrix is singular (-7)",
                 pardisoStep.data());
      return false;
    case -8:
      feErrorMsg(FE_STATUS_ERROR,
                 "PARDISO %s step failed with error message: 32-bit integer overflow problem (-8)",
                 pardisoStep.data());
      return false;
    case -9:
      feErrorMsg(FE_STATUS_ERROR,
                 "PARDISO %s step failed with error message: not enough memory for OOC (-9)",
                 pardisoStep.data());
      return false;
    case -10:
      feErrorMsg(FE_STATUS_ERROR,
                 "PARDISO %s step failed with error message: error opening OOC files (-10)",
                 pardisoStep.data());
      return false;
    case -11:
      feErrorMsg(FE_STATUS_ERROR,
                 "PARDISO %s step failed with error message: read/write error with OOC files (-11)",
                 pardisoStep.data());
      return false;
    case -12:
      feErrorMsg(FE_STATUS_ERROR,
                 "PARDISO %s step failed with error message: (pardiso_64 only) "
                 "pardiso_64 called from 32-bit library (-12)",
                 pardisoStep.data());
      return false;
    case -13:
      feErrorMsg(FE_STATUS_ERROR,
                 "PARDISO %s step failed with error message: interrupted by the "
                 "(user-defined) mkl_progress function (-13)",
                 pardisoStep.data());
      return false;
    case -15:
      feErrorMsg(FE_STATUS_ERROR,
                 "PARDISO %s step failed with error message: internal error which can appear for "
                 "iparm[23]=10 "
                 "and iparm[12]=1. Try switch matching off (set iparm[12]=0 and rerun.) (-15)",
                 pardisoStep.data());
      return false;
    default:
      feErrorMsg(FE_STATUS_ERROR, "Unexpected error code from PARDISO! This should NOT happen.",
                 pardisoStep.data());
      return false;
  }
}

bool feLinearSystemMklPardiso::solve(double *normDx, double *normResidual, double *normAxb,
                                     int *nIter)
{
  std::string pardisoStep;

  tic();
  // Symbolic factorization and allocation
  // (only once or when the matrix sparsity changes)
  if(symbolicFactorization) mklSymbolicFactorization();
  pardisoStep = "symbolic factorization";
  if(!checkPardisoErrorCode(ERROR, pardisoStep)) return false;

  // Actual matrix factorization
  if(recomputeMatrix) mklFactorization();
  pardisoStep = "factorization";
  if(!checkPardisoErrorCode(ERROR, pardisoStep)) return false;

  mklSolve();
  pardisoStep = "solve";
  if(!checkPardisoErrorCode(ERROR, pardisoStep)) return false;
  feInfoCond(FE_VERBOSE > 1, "\t\t\t\tSolved linear system with MKL Pardiso in %f s", toc());

  symbolicFactorization = false;

  *normDx = vectorMaxNorm(matrixOrder, du);
  *normResidual = vectorMaxNorm(matrixOrder, residu);
  *normAxb = 0.0;
  *nIter = 1;

  return true;
}

void feLinearSystemMklPardiso::setToZero()
{
  if(recomputeMatrix) setMatrixToZero();
  setResidualToZero();
}

void feLinearSystemMklPardiso::setMatrixToZero()
{
  for(feInt i = 0; i < nz; i++) Ax[i] = 0;
}

void feLinearSystemMklPardiso::setResidualToZero()
{
  for(feInt i = 0; i < matrixOrder; i++) residu[i] = 0;
}

void feLinearSystemMklPardiso::assemble(feSolution *sol)
{
  // printf("ICI   feLinearSystemMklPardiso::assemble\n");
  if(recomputeMatrix) assembleMatrices(sol);
  assembleResiduals(sol);
}

// ====================================================================
// Les méthodes privées
// mklSymbolicFactorization : Factorisation de Pardiso (symbolique)
//							  préparation des structures et réduire
//                            le remplissage
// mklFactorization         : Factorisation de Pardiso (réelle)
// mklSolve 				: Descente et montée triangulaire
// ====================================================================
void feLinearSystemMklPardiso::mklSymbolicFactorization(void)
{
  // for(feInt i=0;i<matrixOrder;i++) Ax[i] = 1.0;  // dubitatif
  PHASE = 11;
  IPARM[12] = iparm12;
  pardiso_64(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, Ax, Ap, Aj, &IDUM, &NRHS, IPARM, &MSGLVL,
             &DDUM, &DDUM, &ERROR);
}

void feLinearSystemMklPardiso::mklFactorization(void)
{
  PHASE = 22;
  pardiso_64(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, Ax, Ap, Aj, &IDUM, &NRHS, IPARM, &MSGLVL,
             &DDUM, &DDUM, &ERROR);
}

void feLinearSystemMklPardiso::mklSolve(void)
{
  PHASE = 33;
  IPARM[7] = 1;

  // Reset solution vector
  for(feInt i = 0; i < matrixOrder; i++) du[i] = 0.0;

  pardiso_64(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, Ax, Ap, Aj, &IDUM, &NRHS, IPARM, &MSGLVL,
             residu, du, &ERROR);
}

// ====================================================================
// Destructeur/Constructeur de la classe de base dérivée  MKL PARDISO
// ====================================================================
feLinearSystemMklPardiso::~feLinearSystemMklPardiso(void)
{
  delete[] Ap;
  delete[] Aj;
  delete[] Ax;
  delete[] du;
  delete[] residu;
}

void feLinearSystemMklPardiso::setPivot(int pivot) { IPARM[9] = (feInt)pivot; }

#endif
