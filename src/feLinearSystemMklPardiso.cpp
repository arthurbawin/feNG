#include "feLinearSystemMklPardiso.h"
#include <math.h>

#if defined(HAVE_MKL)
#include "mkl.h"

feLinearSystemMklPardiso::feLinearSystemMklPardiso(std::vector<feBilinearForm *> bilinearForms,
                                                   feMetaNumber *metaNumber, feMesh *mesh)
  : feLinearSystem(bilinearForms, metaNumber, mesh)
{
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
  crsMklPardiso = new feCompressedRowStorageMklPardiso(metaNumber, mesh, _formMatrices, _numMatrixForms);
  // toc();
  nz = crsMklPardiso->getNz();
  Ap = (feMKLPardisoInt *) crsMklPardiso->getAp();
  Aj = (feMKLPardisoInt *) crsMklPardiso->getAj();
  Ax = crsMklPardiso->allocateMatrix();
  crsMklPardiso->zeroMatrix(Ax);
  matrixOrder = (feMKLPardisoInt)crsMklPardiso->getMatrixOrder();
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

  for(feInt i = 0; i < 64; i++) PT[i] = NULL;

  N = (feMKLPardisoInt)matrixOrder;
  NRHS = 1;
  // IPARM[2]= num_procs; // nombre de processeurs
  MAXFCT = 1;
  MNUM = 1;
  MSGLVL = 0; // ou 0
  ERROR = 0;
  //=================================================================
  //  Factorisation symbolique
  //=================================================================
  // mklSymbolicFactorization();
}

double vectorMaxNorm(feInt N, double *V)
{
  double t = 0.0;
  for(feInt i = 0; i < N; i++) {
    t = fmax(t, fabs(V[i]));
  }

  return t;
}

double vectorL2Norm(feInt N, double *V)
{
  double t = 0.0;

#if defined(HAVE_OMP)
#pragma omp parallel for reduction(+ : t) schedule(dynamic)
#endif
  for(feInt i = 0; i < N; i++){
    t += V[i] * V[i];
  }
  return pow(t, 0.5);
}

void feLinearSystemMklPardiso::print_matrix()
{
  for(feInt i = 0; i < nz; i++) printf("%ld %g \n", i, Ax[i]);
}

void feLinearSystemMklPardiso::assembleMatrices(feSolution *sol)
{
  feInfo("Assembling the Matrix...");
  tic();

  // feInt NumberOfBilinearForms = _formMatrices.size();

  feInfo("Repartition des NNZ: ");
  for(int i = 0; i < matrixOrder; ++i){
    feInfo("NNZ[%d] = %d", i, crsMklPardiso->getNnz()[i]);
  }
  feInfo("Repartition de Ap: ");
  for(int i = 0; i < matrixOrder+1; ++i){
    feInfo("Ap[%d] = %d", i, Ap[i]);
  }
  // feInfo("Repartition de Aj: ");
  // for(int i = 0; i < matrixOrder; ++i){
  //   feInfo("NNZ[%d] = %d", i, crsMklPardiso->getNnz()[i]);
  // }

  feInfo("PRINTING INFO");
  crsMklPardiso->print_info();


  for(feInt eq = 0; eq < _numMatrixForms; eq++) {
    feBilinearForm *f = _formMatrices[eq];
    feCncGeo *cnc = f->getCncGeo();
    int nbColor = cnc->getNbColor();
    std::vector<int> &nbElmPerColor = cnc->getNbElmPerColor();
    std::vector<std::vector<int> > &listElmPerColor = cnc->getListElmPerColor();

    int nbElmC; // nb elm of the same color
    std::vector<int> listElmC; // list elm of the same color;

    for(int iColor = 0; iColor < nbColor; iColor++) {
      nbElmC = nbElmPerColor[iColor]; // nbElmC : nombre d'elm de meme couleur
      listElmC = listElmPerColor[iColor];
      // nbElmC = cnc->getNbElmPerColorI(iColor);
      // listElmC = cnc -> getListElmPerColorI(iColor);

      feInt numThread = 0;
      int elm;
      int eqt;

      double **Ae;
      feInt nRow;
      feInt nColumn;
      std::vector<feInt> Row;
      std::vector<feInt> Column;
      feInt I;
      feInt J;

      feInt debut;
      feInt fin;
      feInt ncf;
      std::vector<feInt> irangee;

#if defined(HAVE_OMP)
#pragma omp parallel for private(numThread, elm, eqt, f, nRow, nColumn, Row, Column, Ae, I, J, debut, fin, ncf,       \
                                 irangee) schedule(dynamic)
#endif
      for(int iElm = 0; iElm < nbElmC; ++iElm) {
#if defined(HAVE_OMP)
        numThread = omp_get_thread_num();
        eqt = eq + numThread * _numMatrixForms;
        f = _formMatrices[eqt];
#endif
        elm = listElmC[iElm];
        f->computeMatrix(_metaNumber, _mesh, sol, elm);

        nRow = f->getNiElm();
        Row = f->getAdrI();

        feInfo("elm %d SIZE OF ROW is %d - matrixORder = %d", iElm, Row.size(), matrixOrder);
        for(auto val : Row)
          feInfo("val = %d", val);

        nColumn = f->getNjElm();
        Column = f->getAdrJ();
        Ae = f->getAe();

        feInfo("ELEM ADRI & ADRJ");
        feInfo("%d - %d - %d", Row[0], Row[1], Row[2]);
        feInfo("%d - %d - %d", Column[0], Column[1], Column[2]);

        for(feInt i = 0; i < nRow; i++) {
          irangee.resize(matrixOrder);
          I = Row[i];
          debut = 0;
          fin = 0;
          ncf = 0;
          if(I < matrixOrder) {
            debut = Ap[I] - 1;
            fin = Ap[I + 1] - 1;
            ncf = fin - debut;

            feInfo("Ligne %d", I);
            feInfo("DEBUT = %d", debut);
            feInfo("Fin = %d", fin);
            feInfo("NCF = %d", ncf);
            feInfo("NZ = %d au total", nz);

            for(feInt j = 0; j < ncf; j++) {
              // feInfo("Accessing entry %d in vector of size %d", debut + j, nz);
              // feInfo("Result in Aj is %d ", Aj[debut + j]);
              // feInfo("irangee.size = %d ", irangee.size());
              // feInfo("irangee[0] = %d ", irangee[0]);
              // feInfo("accessing %d ", Aj[debut + j] - 1);
              int indexx = Aj[debut + j] - 1;
              // feInfo("Result in irangee is %d ", irangee[indexx]);
              // irangee[Aj[debut + j] - 1] = debut + j;
              irangee[indexx] = debut + j;
              // feInfo("Result in irangee is %d ", irangee[Aj[debut + j] - 1]);
            }

            for(feInt j = 0; j < nColumn; j++) {
              J = Column[j];
              if(J < matrixOrder) {
                Ax[irangee[J]] += Ae[i][j];
              }
            }
          }
        }
      }
    }
  }
  // print_matrix();
  // double res = 0.0;
  // for(feInt i = 0; i < nz; i++) res += fabs(Ax[i]);
  // feInfo("sumMatrix = %f", res);
  feInfo("Done...");
  toc();
}

void feLinearSystemMklPardiso::assembleResiduals(feSolution *sol)
{
  feInfo("Assembling the residual...");
  tic();

  for(feInt eq = 0; eq < _numResidualForms; eq++) {
    feBilinearForm *f = _formResiduals[eq];
    feCncGeo *cnc = f->getCncGeo();
    int nbColor = cnc->getNbColor();
    std::vector<int> &nbElmPerColor = cnc->getNbElmPerColor();
    std::vector<std::vector<int> > &listElmPerColor = cnc->getListElmPerColor();
    // feInfo("Looping over %d colors", nbColor);

    int nbElmC;
    std::vector<int> listElmC;

    for(int iColor = 0; iColor < nbColor; ++iColor) {
      nbElmC = nbElmPerColor[iColor]; // nbElmC : nombre d'elm de meme couleur
      listElmC = listElmPerColor[iColor];

      feInt numThread = 0;
      int elm;
      int eqt;
      // feBilinearForm *f_t;

      double *Be;
      feInt nRow;
      std::vector<feInt> Row;

#if defined(HAVE_OMP)
#pragma omp parallel for private(numThread, elm, eqt, nRow, Row, Be, f) schedule(dynamic)
#endif
      for(int iElm = 0; iElm < nbElmC; ++iElm) {
#if defined(HAVE_OMP)
        numThread = omp_get_thread_num();
        eqt = eq + numThread * _numResidualForms;
        f = _formResiduals[eqt];
#endif
        elm = listElmC[iElm];
        f->computeResidual(_metaNumber, _mesh, sol, elm);

        nRow = f->getNiElm();
        Row = f->getAdrI();

        Be = f->getBe();

        for(feInt i = 0; i < nRow; i++) {
          if(Row[i] < matrixOrder) residu[Row[i]] += Be[i];
        }
      }
    }
  }
  toc();
  feInfo("Done...");
  // toc();
  // for(int i=0;i<matrixOrder;i++){
  //   printf("%g \n",residu[i]);
  // }
  // double res = 0.0;
  // for(feInt i = 0; i < matrixOrder; i++){
  //   res += fabs(residu[i]);
  // }
  // feInfo("sumResidu = %f", res);
  // feInfo("nz = %d", nz);
  // feInfo("matrixOrder = %d", matrixOrder);
  // feInfo("max threads = %d", omp_get_max_threads());
  // feInfo("num threads = %d", omp_get_num_threads());
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

void feLinearSystemMklPardiso::solve(double *normDx, double *normResidual, double *normAxb,
                                     int *nIter)
{
  feInfo("Solving ...");
  if(symbolicFactorization) mklSymbolicFactorization();
  if(recomputeMatrix) {
    // tic();
    mklFactorization();
    // toc();
  }
  mklSolve();
  feInfo("Done.");
  symbolicFactorization = false;

  // feInfo("matrixOrder : %d",matrixOrder);

  // for (int i=0;i<matrixOrder;i++)
  //   feInfo("%g",du[i]);

  *normDx = vectorMaxNorm(matrixOrder, du);
  *normResidual = vectorMaxNorm(matrixOrder, residu);
  *normAxb = 0.0;
  *nIter = 0;
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
  if(ERROR != 0) {
    printf("feLinearSystemMklPardiso::mklSymbolicFactorization - erreur %ld\n", ERROR);
  }
}

void feLinearSystemMklPardiso::mklFactorization(void)
{
  PHASE = 22;
  pardiso_64(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, Ax, Ap, Aj, &IDUM, &NRHS, IPARM, &MSGLVL,
             &DDUM, &DDUM, &ERROR);
  if(ERROR != 0) {
    printf("eLinearSystemMklPardiso::mklFactorization- erreur %ld\n", ERROR);
  }
}

void feLinearSystemMklPardiso::mklSolve(void)
{
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
feLinearSystemMklPardiso::~feLinearSystemMklPardiso(void)
{
  delete crsMklPardiso;
  delete[] Ax;
  delete[] du;
  delete[] residu;
}

void feLinearSystemMklPardiso::setPivot(int pivot) { IPARM[9] = (feMKLPardisoInt)pivot; }

#endif
