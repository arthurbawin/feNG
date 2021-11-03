#include "feLinearSystemMklPardiso.h"
#include <math.h>

#if defined(HAVE_MKL)
#include "mkl.h"

double vectorMaxNorm(feInt N, double *V) {
  double t = 0.0;

  for(feInt i = 0; i < N; i++) t = fmax(t, abs(V[i]));

  return t;
}

double vectorL2Norm(feInt N, double *V) {
  double t = 0.0;

  for(feInt i = 0; i < N; i++) t += V[i] * V[i];

  return pow(t, 0.5);
}

void feLinearSystemMklPardiso::print_matrix() {
  for(feInt i = 0; i < nz; i++) printf("%ld %g \n", i, Ax[i]);
}

void feLinearSystemMklPardiso::assembleMatrices(feSolution *sol) {
  feInt NumberOfBilinearForms = _formMatrices.size();
  // printf("ICI feLinearSystemMklPardiso::assembleMatrices  %ld\n", NumberOfBilinearForms);
  for(feInt eq = 0; eq < NumberOfBilinearForms; eq++) {
    feBilinearForm *equelm = _formMatrices[eq];
    feInt cncGeoTag = equelm->getCncGeoTag();
    feInt nbElems = _mesh->getNbElm(cncGeoTag);
    for(feInt e = 0; e < nbElems; e++) {
      equelm->initialize(_metaNumber, _mesh, sol, e);

      equelm->computeMatrix(_metaNumber, _mesh, sol, e);

      feInt nRow = equelm->getNiElm();
      std::vector<feInt> Row(equelm->getAdrI().begin(), equelm->getAdrI().end());
      feInt nColumn = equelm->getNjElm();
      std::vector<feInt> Column(equelm->getAdrJ().begin(), equelm->getAdrJ().end());
      double **Ae = equelm->getAe();

      crsMklPardiso->matrixAddValues(Ax, nRow, Row.data(), nColumn, Column.data(), Ae);
    }
  }
}

void feLinearSystemMklPardiso::assembleResiduals(feSolution *sol) {
  feInt NumberOfBilinearForms = _formResiduals.size();

  for(feInt eq = 0; eq < NumberOfBilinearForms; eq++) {
    feBilinearForm *equelm = _formResiduals[eq];
    feInt cncGeoTag = equelm->getCncGeoTag();
    feInt nbElems = _mesh->getNbElm(cncGeoTag);
    for(feInt e = 0; e < nbElems; e++) {
      equelm->initialize(_metaNumber, _mesh, sol, e);

      equelm->computeResidual(_metaNumber, _mesh, sol, e);

      feInt nRow = equelm->getNiElm();
      std::vector<feInt> Row(equelm->getAdrI().begin(), equelm->getAdrI().end());
      double *Be = equelm->getBe();

      for(feInt i = 0; i < nRow; i++)
        if(Row[i] < matrixOrder) residu[Row[i]] += Be[i];
    }
  }
}

void feLinearSystemMklPardiso::assignResidualToDCResidual(feSolutionContainer *solContainer) {
  for(feInt i = 0; i < matrixOrder; ++i) solContainer->_fResidual[0][i] = residu[i];
}

void feLinearSystemMklPardiso::applyCorrectionToResidual(double coeff, std::vector<double> &d) {
  for(int i = 0; i < matrixOrder; ++i) residu[i] += coeff * d[i];
}

void feLinearSystemMklPardiso::correctSolution(feSolution *sol) {
  // Est-ce efficace?
  // Pourquoi ne pas avoir un pointeur sur le vecteur solution;
  // Conversion de feInt à int?
  for(feInt i = 0; i < matrixOrder; i++) sol->incrementSolAtDOF((int)i, du[i]);
}

void feLinearSystemMklPardiso::correctSolution(double *sol) {
  for(feInt i = 0; i < matrixOrder; i++) sol[i] += du[i];
}

void feLinearSystemMklPardiso::solve(double *normDx, double *normResidual, double *normAxb,
                                     int *nIter) {
  if(symbolicFactorization) mklSymbolicFactorization();
  if(recomputeMatrix) {
    tic();
    mklFactorization();
    toc();
  }
  mklSolve();

  symbolicFactorization = false;

  *normDx = vectorMaxNorm(matrixOrder, du);
  *normResidual = vectorMaxNorm(matrixOrder, residu);
  *normAxb = 0.0;
  *nIter = 0;
}

void feLinearSystemMklPardiso::setToZero() {
  if(recomputeMatrix) setMatrixToZero();
  setResidualToZero();
}

void feLinearSystemMklPardiso::setMatrixToZero() {
  for(feInt i = 0; i < nz; i++) Ax[i] = 0;
}

void feLinearSystemMklPardiso::setResidualToZero() {
  for(feInt i = 0; i < matrixOrder; i++) residu[i] = 0;
}

void feLinearSystemMklPardiso::assemble(feSolution *sol) {
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
void feLinearSystemMklPardiso::mklSymbolicFactorization(void) {
  // for(feInt i=0;i<matrixOrder;i++) Ax[i] = 1.0;  // dubitatif
  PHASE = 11;
  IPARM[12] = iparm12;
  pardiso_64(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, Ax, Ap, Aj, &IDUM, &NRHS, IPARM, &MSGLVL,
             &DDUM, &DDUM, &ERROR);
  if(ERROR != 0) {
    printf("feLinearSystemMklPardiso::mklSymbolicFactorization - erreur %ld\n", ERROR);
  }
}

void feLinearSystemMklPardiso::mklFactorization(void) {
  PHASE = 22;
  pardiso_64(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, Ax, Ap, Aj, &IDUM, &NRHS, IPARM, &MSGLVL,
             &DDUM, &DDUM, &ERROR);
  if(ERROR != 0) { printf("eLinearSystemMklPardiso::mklFactorization- erreur %ld\n", ERROR); }
}

void feLinearSystemMklPardiso::mklSolve(void) {
  PHASE = 33;
  IPARM[7] = 1;
  for(feInt i = 0; i < matrixOrder; i++) du[i] = 0.0;

  pardiso_64(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, Ax, Ap, Aj, &IDUM, &NRHS, IPARM, &MSGLVL,
             residu, du, &ERROR);
  if(ERROR != 0) {
    printf("feLinearSystemMklPardiso::mklSolve - erreur %ld\n", ERROR);
    exit(-1);
  }
}

// ====================================================================
// Destructeur/Constructeur de la classe de base dérivée  MKL PARDISO
// ====================================================================
feLinearSystemMklPardiso::~feLinearSystemMklPardiso(void) {
  delete crsMklPardiso;
  delete[] Ax;
  delete[] du;
  delete[] residu;
}

feLinearSystemMklPardiso::feLinearSystemMklPardiso(std::vector<feBilinearForm *> &formMatrices,
                                                   std::vector<feBilinearForm *> &formResiduals,
                                                   feMetaNumber *metaNumber, feMesh *mesh)
  : feLinearSystem(formMatrices, formResiduals, metaNumber, mesh) {
  // long int cnt = 0;
  // #pragma omp parallel for private(cnt)
  // for(int i = 0; i < 100; ++i){
  //   // printf("Printing %3d from thread %d\n", i, omp_get_thread_num());
  //   printf("Thread %d has printed %2d times\n", omp_get_thread_num(), cnt++);
  // }

  recomputeMatrix = true;
  //=================================================================
  // Structure Creuse CSR de MKL
  //=================================================================
  // tic();
  crsMklPardiso = new feCompressedRowStorageMklPardiso(metaNumber, mesh, _formMatrices);
  // toc();
  nz = crsMklPardiso->getNz();
  Ap = (feMKLPardisoInt *)crsMklPardiso->getAp();
  Aj = (feMKLPardisoInt *)crsMklPardiso->getAj();
  Ax = crsMklPardiso->allocateMatrix();
  crsMklPardiso->zeroMatrix(Ax);
  matrixOrder = (feMKLPardisoInt)crsMklPardiso->getMatrixOrder();
  du = new double[matrixOrder];
  residu = new double[matrixOrder];
  symbolicFactorization = true;
  //=================================================================
  // 	INITIALISATION - MATRICE NON SYMETRIQUE
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

  for(feInt i = 0; i < 64; i++) PT[i] = NULL;

  N = (feMKLPardisoInt)matrixOrder;
  NRHS = 1;
  // IPARM[2]= num_procs; // nombre de processeurs
  MAXFCT = 1;
  MNUM = 1;
  MSGLVL = 0; // ou 1
  ERROR = 0;
  //=================================================================
  // 	Factorisation symbolique
  //=================================================================
  // mklSymbolicFactorization();
}

void feLinearSystemMklPardiso::setPivot(int pivot) { IPARM[9] = (feMKLPardisoInt)pivot; }

#endif