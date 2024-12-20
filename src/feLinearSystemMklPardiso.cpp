#include "feLinearSystem.h"

extern int FE_VERBOSE;

#if defined(HAVE_MKL)
#include "mkl.h"

#if defined(HAVE_MPI)
  #include "mpi.h"
  #include "mkl_cluster_sparse_solver.h"
#endif

//
// Determine the matrix, unknown vector and RHS distribution
// across MPI processes once the sparsity pattern is known
//
void feLinearSystemMklPardiso::setOwnershipAndAllocate(void)
{
#if defined(HAVE_MPI)

  // If _ownershipSplit = -1, then the linear system is stored
  // on a single process, not distributed.
  if(_ownershipSplit >= 0) {

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Number of rows
    PardisoInt N = _nInc;
    PardisoInt low = 0, high;

    if(_ownershipSplit != 0 && _ownershipSplit != 1 && _ownershipSplit != 2) {
      _ownershipSplit = 0;
    }

    switch(_ownershipSplit) {
      case 0:
        // Even number of rows per proc according to PetscSplitOwnership()'s formula
        {
          if(N < size) {
            // There are more MPI processes than unknowns.
            // Cluster Pardiso supports idle processes (with 0 matrix rows),
            // but the situation is unlikely to happen so we'll just ask 
            // to run with less processes
            feErrorMsg(FE_STATUS_ERROR, "There are more MPI processes (%d) than FE unknowns! (%d)\n"
              "Idle processes storing 0 rows are not supported at the moment, please run again"
              " with less MPI processes.", size, N);
            exit(-1);
          }
          _numOwnedRows = N/size + ((N % size) > rank);
          PardisoInt sendbuf[1] = { _numOwnedRows };

          // Synchronize the number of rows on each process
          PardisoInt *rowsOnEachProc = (PardisoInt*) malloc(size*sizeof(PardisoInt));
          MPI_Allgather(sendbuf, 1, MPI_INT, rowsOnEachProc, 1, MPI_INT, MPI_COMM_WORLD);

          // Determine lower and upper bounds on this process
          high = -1;
          for(int i = 0; i <= rank; ++i) {
            low = high+1;
            high += rowsOnEachProc[i];
          }
          free(rowsOnEachProc);

          // feInfo("proc %d : nrows_total = %d - nrows_onproc = %d - low = %d - high = %d", rank, N, _numOwnedRows, low, high);

          // Synchronize the bounds on each process
          PardisoInt localLow[1] = {low}, localHigh[1] = {high};
          _ownedUpperBounds = (PardisoInt*) malloc(size*sizeof(PardisoInt));
          _ownedLowerBounds = (PardisoInt*) malloc(size*sizeof(PardisoInt)); 
          MPI_Allgather(localLow , 1, MPI_INT, _ownedLowerBounds, 1, MPI_INT, MPI_COMM_WORLD);
          MPI_Allgather(localHigh, 1, MPI_INT, _ownedUpperBounds, 1, MPI_INT, MPI_COMM_WORLD);

          _numOwnedRowsOnAllProc = (PardisoInt*) malloc(size*sizeof(PardisoInt));
          for(int i = 0; i < size; ++i) {
            _numOwnedRowsOnAllProc[i] = (_ownedUpperBounds[i] - _ownedLowerBounds[i]) + 1;
          }

          // Determine the number of nonzero on this proc knowing _numOwnedRows
          _numOwnedNNZ = _mat_ia[high+1] - _mat_ia[low];

          // Create the local ia, ja and values arrays
          PardisoInt *ia_loc = (PardisoInt*) malloc((_numOwnedRows + 1)*sizeof(PardisoInt));
          PardisoInt *ja_loc = (PardisoInt*) malloc(_numOwnedNNZ*sizeof(PardisoInt));
          double  *a_loc = (double*) malloc(_numOwnedNNZ*sizeof(double));

          // ia
          // printf("proc %d ia = ", rank);
          for(int i = 0; i < _numOwnedRows; ++i) {
            ia_loc[i] = _mat_ia[low + i] - _mat_ia[low];
            // printf("%d ", ia_loc[i]);
          }
          ia_loc[_numOwnedRows] = _numOwnedNNZ;// nnz
          // printf("%d", ia_loc[_numOwnedRows]);
          // printf("\n");

          // ja and matrix entries are simply the local portion of the global arrays
          PardisoInt offset = _mat_ia[low];
          for(int i = 0; i < _numOwnedNNZ; ++i) {
            ja_loc[i] = _mat_ja[offset + i];
             a_loc[i] = _mat_values[offset + i];
          }

          // Owned part of the solution array for update
          _ownedSolution = (double*) malloc(_numOwnedRows * sizeof(double));

          // printf("proc %d ja = ", rank);
          // for(int i = 0; i < _numOwnedNNZ; ++i) {
          //   printf("%d ", ja_loc[i]);
          // }
          // printf("\n");
          // printf("proc %d  a = ", rank);
          // for(int i = 0; i < _numOwnedNNZ; ++i) {
          //   printf("%d ", a_loc[i]);
          // }
          // printf("\n");

          // Set the matrix arrays
          // delete[] _mat_ia;
          // delete[] _mat_ja;
          // delete[] _mat_values;
          free(_mat_ia);
          free(_mat_ja);
          free(_mat_values);
          _mat_ia = ia_loc;
          _mat_ja = ja_loc;
          _mat_values = a_loc;

          // Create the vector solution and RHS arrays
          du = new double[_numOwnedRows];
          residu = new double[_numOwnedRows];

          break;
        }
      case 1:
        // TODO: Even number of nnz per proc, possible row overlap
        // The number of nnz is split according to PetscSplitOwnership()'s formula
        // PardisoInt nNNZPerProc = _nnz/size + ((_nnz % size) > rank);
        feErrorMsg(FE_STATUS_ERROR, "Even split of nonzero not yet implemented. Choose ownershipSplit = 0.");
        exit(-1);
        break;
      case 2:
        // TODO: Even number of nnz per proc (mostly), without row overlap
        // Il faut:
        // Estimer le nombre de nnz per proc comme dans le cas 1
        // Synchroniser (sûrement) sur chaque proc la valeur de sur les autres
        // Progresser dans ia/ja pour voir s'il y a des sauts de lignes
        // Si la valeur précédente dans ja est plus grande, alors
        // il y a eu un saut de ligne et c'est suffisant pour accepter le découpage.
        // Sinon, on ne peut pas conclure et il faut sûrement regarder la valeur dans ia.
        feErrorMsg(FE_STATUS_ERROR, "Even split of nonzero not yet implemented. Choose ownershipSplit = 0.");
        exit(-1);
        break;
    }
  } else { // if _ownershipSplit > 0 
#endif
  _numOwnedNNZ = _nnz;
  _numOwnedRows = _nInc;
  _numOwnedRowsOnAllProc = (PardisoInt*) malloc(1*sizeof(PardisoInt));
  _numOwnedRowsOnAllProc[0] = _nInc;
  _ownedLowerBounds = (PardisoInt*) malloc(1*sizeof(PardisoInt));
  _ownedLowerBounds[0] = 0; 
  _ownedUpperBounds = (PardisoInt*) malloc(1*sizeof(PardisoInt));
  _ownedUpperBounds[0] = _nInc-1;
  _ownedSolution = (double*) malloc(_nInc * sizeof(double));
  // Leave the matrix arrays as-is (ia, ja, values),
  // and allocate solution vector (du) and RHS (residu)
  du = new double[_nInc];
  residu = new double[_nInc];

#if defined(HAVE_MPI)
  }
#endif
}

feLinearSystemMklPardiso::feLinearSystemMklPardiso(std::vector<feBilinearForm *> bilinearForms,
                                                   int numUnknowns,
                                                   int ownershipSplit)
  : feLinearSystem(bilinearForms), _nInc(numUnknowns), _ownershipSplit(ownershipSplit)
{
  recomputeMatrix = true;
  symbolicFactorization = true;

  //
  // Determine sparsity pattern, allocate MKL Pardiso arrays
  // for sequential or distributed matrix.
  //
  feEZCompressedRowStorage _EZCRS(numUnknowns, _formMatrices, _numMatrixForms);
  _nnz = _EZCRS.getNumNNZ();
  // _mat_ia = new PardisoInt[_nInc + 1]; // ia in pardiso doc (beginning of row in ja)
  // _mat_ja = new PardisoInt[_nnz]; // ja in pardiso doc (columns)
  // _mat_values = new double[_nnz]; // a in pardiso doc (matrix coefficients)
  _mat_ia = (PardisoInt*) malloc((_nInc + 1)*sizeof(PardisoInt)); // ia in pardiso doc (beginning of row in ja)
  _mat_ja = (PardisoInt*) malloc( _nnz      *sizeof(PardisoInt)); // ja in pardiso doc (columns)
  _mat_values = (double*) malloc( _nnz      *sizeof(double)); // a in pardiso doc (matrix coefficients)  
  _EZCRS.get_ia_Pardiso(&_mat_ia);
  _EZCRS.get_ja_Pardiso(&_mat_ja);

  for(feInt i = 0; i < _nnz; i++) _mat_values[i] = 0.;

  this->setOwnershipAndAllocate();

  //
  // Initialize MKL Pardiso parameters
  //
  // Matrix type: 11 = real nonsymmetric
  MTYPE = 11;
  // Zero pivots are perturbed to 10^(-IPIVOT) 
  IPIVOT = 13;

  for(feInt i = 0; i < 64; i++) {
    IPARM[i] = 0;
    PT[i] = 0;
  }

#if defined(HAVE_MPI)
  //
  // Option for distributed Cluster MKL Pardiso
  //
  IPARM[ 0] =  1; /* Solver default parameters overriden with provided by iparm */

  /* Reordering algorithm to reduce fill-in:
    2 : Nested dissection provided by METIS (default)
    3 : Multithreaded version of 2
    10: MPI distributed version of 2. Only for distributed matrix (iparm[39] > 0),
        Not working if either of the following is true:
        iparm[10] = 1 (Scaling)
        iparm[12] = 1 (Maximum Weighted Matching algorithm)
        + [iparm[30] > 0 (Partial solve)            ]
          [iparm[35] > 0 (Schur complement control) ]
          [iparm[36] > 1 (BSR format)               ]
  */
  IPARM[ 1] =  3;

  if(IPARM[1] == 10 && _ownershipSplit < 0) {
    feErrorMsg(FE_STATUS_ERROR, "iparm[1] = 10 is for distributed matrices only - Incompatible with ownershipSplit = -1");
    exit(-1);
  }

  IPARM[ 5] =  0; /* Write solution into x */
  IPARM[ 7] =  1; /* Max number of iterative refinement steps */
  IPARM[ 9] = IPIVOT; /* Perturb the pivot elements with 1E-13 */
  IPARM[10] =  0; /* Use or not nonsymmetric permutation and scaling MPS */

  // !!! Currently, turning IPARM[12] to 0 causes crashes with more than 1 OMP thread :/ !!!
  IPARM[12] =  1; /* Switch on Maximum Weighted Matching algorithm (default for non-symmetric) */
    
  IPARM[17] = -1; /* Output: Number of nonzeros in the factor LU */
  IPARM[18] = -1; /* Output: Mflops for LU factorization */
  IPARM[26] =  1; /* Check input data for correctness */
  IPARM[34] =  1; /* 0-based indexing - (0) Fortran-style (1-based) - (1) C-style (0-based) */

  if(_ownershipSplit < 0) {
    IPARM[39] =  0; /* Input: matrix/rhs/solution are stored on master process  */
  } else {
    IPARM[39] =  2; /* Input: matrix/rhs/solution are distributed between MPI processes  */

    /* Beginning and ending of input domain. 
    The number of the matrix A row, RHS element, and, for iparm[39]=2, solution vector
    that begins the input domain belonging to this MPI process.
    Only applicable to the distributed assembled matrix input format (iparm[39]> 0).

    If IPARM[40] > IPARM[41], then this rank does not have any part of the domain.
    This is currently unused, i.e., all MPI processes store at least 1 row.
    (see cl_solver_unsym_distr_c.c example file) */
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    IPARM[40] = _ownedLowerBounds[rank];
    IPARM[41] = _ownedUpperBounds[rank];
  }

  NRHS = 1;
  MAXFCT = 1;
  MSGLVL = FE_VERBOSE > VERBOSE_MODERATE;
  MNUM = 1;
  ERROR = 0;
#else
  //
  // Option for non-distributed (but multithreaded) MKL Pardiso
  //
  IPARM[0] = 1; /* No solver default */
  IPARM[1] = 3; /* Fill-in reordering: (0)Min Degree, (2)METIS, (3)OpenMP */
  IPARM[2] = 0; /* Reserved. Set to zero for MKL version. Number of threads otherwise */
  IPARM[3] = 0; /* No iterative-direct algorithm */
  IPARM[4] = 0; /* No user fill-in reducing permutation */
  IPARM[5] = 0; /* Write solution into (0)x, (1)b */
  IPARM[6] = 0; /* Not in use - Number of iterative refinement steps performed */
  IPARM[7] = 1; /* Max numbers of iterative refinement steps. (2)=2 iterations */
  IPARM[8] = 0; /* Reserved */
  IPARM[9] = IPIVOT; /* Perturb the pivot elements with 1E-13 */ //-AG
  IPARM[10] = 0; /* (0)Disable scaling, (1) Use nonsymmetric permutation and scaling MPS */ //--AG
  IPARM[11] = 0; /* Not in use - Solve AX=B */

  // !!! Turning IPARM[12] to 0 causes crashes with more than 1 OMP thread :/ !!!
  IPARM[12] = 1; /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */

  IPARM[13] = 0; /* Output: Number of perturbed pivots */
  IPARM[14] = 0; /* Not in use - Peak memory on symbolic factorization*/
  IPARM[15] = 0; /* Not in use - Permanent memory on symbolic factorization */
  IPARM[16] = 0; /* Not in use - Peak memory on numerical factorization */
  IPARM[17] = -1; /* Output: Number of nonzeros in the factor LU */
  IPARM[18] = 0; /* Output: Mflops for LU factorization */
  IPARM[19] = 0; /* Output: Numbers of CG Iterations */
  IPARM[23] = 1; /* Parallel factorization control - (0) sequential - (1) parallel*/
  IPARM[24] = 0; /* Parallel LU solve control - (0) parallel - (1) sequential */
  IPARM[26] = 1; /* Check input data for correctness */
  IPARM[34] = 1; /* 0-based indexing - (0) Fortran-style (1-based) - (1) C-style (0-based) */
  NRHS = 1;
  MAXFCT = 1;
  MSGLVL = FE_VERBOSE > VERBOSE_MODERATE;
  MNUM = 1;
  ERROR = 0;
#endif
  
  feInfoCond(FE_VERBOSE > 0, "\t\tCreated a MKL Pardiso linear system of size %d x %d", _nInc, _nInc);
  feInfoCond(FE_VERBOSE > 0, "\t\tNumber of nonzero elements: %d", _nnz);
#if defined(HAVE_MPI)
  feInfoCond(FE_VERBOSE > 0, "\t\tSolving with cluster version of MKL Pardiso because MPI is enabled");
  if(_ownershipSplit < 0) {
    feInfoCond(FE_VERBOSE > 0, "\t\tLinear system storage: sequential (stored on master process)");
  } else {
    feInfoCond(FE_VERBOSE > 0, "\t\tLinear system storage: distributed among MPI processes (strategy %d)", _ownershipSplit);
  }
#else
  feInfoCond(FE_VERBOSE > 0, "\t\tSolving with non-cluster version of MKL Pardiso because MPI is disabled");
#endif
  feInfoCond(FE_VERBOSE > 0, "\t\t\tTODO: Add Pardiso options");
}

//
// Get the maximum absolute value of a vector on all MPI processes
//
double vectorMaxNorm(int numOwnedRows, double *vec)
{
  double my_max = 0.;
  for(feInt i = 0; i < numOwnedRows; i++) {
    my_max = fmax(my_max, fabs(vec[i]));
  }
#if defined(HAVE_MPI)
  double global_max = 0;
  MPI_Allreduce(&my_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return global_max;
#else
  return my_max;
#endif
}

//
// Print the matrix on each process with a layout similar to PETSc's
//
void feLinearSystemMklPardiso::viewMatrix()
{
  if(_displayMatrixInConsole) {
#if defined(HAVE_MPI)
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int ranktoprint = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    while(ranktoprint < size) {
      if(rank == ranktoprint) {
        // feInfo("Matrix on process %d:", rank);
#else
        feInfo("Matrix:");
#endif
    // for(feInt i = 0; i < _nnz; i++) printf("%ld %g \n", i, _mat_values[i]);
    for(feInt i = 0; i < _numOwnedRows; i++) {
      printf("Row %ld: ", i);
      int debut = _mat_ia[i];
      int fin = _mat_ia[i + 1];
      for(feInt j = 0; j < fin - debut; ++j) {
        printf("(%d, %g)  ", _mat_ja[debut + j], _mat_values[debut + j]);
      }
      printf("\n");
    }
#if defined(HAVE_MPI)
      }
      ranktoprint++;
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
  }
}

//
// Print the RHS on each process
//
void feLinearSystemMklPardiso::viewRHS()
{
  if(_displayRHSInConsole) {
#if defined(HAVE_MPI)
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int ranktoprint = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    while(ranktoprint < size) {
      if(rank == ranktoprint) {
        feInfo("RHS on process %d:", rank);
#else
        feInfo("RHS:");
#endif
        for(int i = 0; i < _numOwnedRows; ++i) {
          feInfo("RHS[%d] = %+-1.16e", i, residu[i]);
        }
#if defined(HAVE_MPI)
      }
      ranktoprint++;
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
  }
}

//
// Print the solution vector (the Newton residual du) on each process
//
void feLinearSystemMklPardiso::viewResidual()
{
#if defined(HAVE_MPI)
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int ranktoprint = 0;
  MPI_Barrier(MPI_COMM_WORLD);
  while(ranktoprint < size) {
    if(rank == ranktoprint) {
      feInfo("Residual (solution vector) on process %d:", rank);
#else
      feInfo("Residual (solution vector) :");
#endif
      for(int i = 0; i < _numOwnedRows; ++i) {
        feInfo("du[%d] = %+-1.16e", i, du[i]);
      }
#if defined(HAVE_MPI)
    }
    ranktoprint++;
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif
}

void feLinearSystemMklPardiso::assembleMatrices(feSolution *sol)
{
  tic();

#if defined(HAVE_MPI)
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Assemble on all procs if system is distributed
  // OR
  // Assemble on master proc only if system is sequential
  if(_ownershipSplit >= 0 || (_ownershipSplit < 0 && rank == 0)) {
    feInfo("Assembling mat with %d and proc %d", _ownershipSplit, rank);
#else
  int rank = 0;
#endif

  PardisoInt low = _ownedLowerBounds[rank];
#if defined(HAVE_MPI)
  PardisoInt high = _ownedUpperBounds[rank];
#endif

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
          if(adrI[i] < _nInc) niElm.push_back(i);
        }
        sizeJ = adrJ.size();
        njElm.reserve(sizeJ);
        for(feInt i = 0; i < sizeJ; ++i) {
          if(adrJ[i] < _nInc) njElm.push_back(i);
        }

        adrI.erase(std::remove_if(adrI.begin(), adrI.end(),
                                  [this](const int &x) { return x >= this->_nInc; }),
                   adrI.end());
        adrJ.erase(std::remove_if(adrJ.begin(), adrJ.end(),
                                  [this](const int &x) { return x >= this->_nInc; }),
                   adrJ.end());

#if defined(HAVE_MPI)
        // Remove row indices in adrI that are not owned by this proc
        auto it_nielm = niElm.begin();
        for(auto it_i = adrI.begin();  it_i != adrI.end(); )
        {
          if(low <= (*it_i) && (*it_i) <= high) {
            // In owned range of rows : do nothing
            ++it_i;
            ++it_nielm;
          } else {
            // Erase in adrI only
            it_i = adrI.erase(it_i);
            it_nielm = niElm.erase(it_nielm);
          }
        }
#endif

        // Flatten Ae at relevant indices
        sizeI = adrI.size();
        sizeJ = adrJ.size();

        for(feInt i = 0; i < sizeI; i++) {
          feInt I = adrI[i] - low;
          feInt debut = _mat_ia[I];
          feInt fin = _mat_ia[I + 1];
          feInt numColumns = fin - debut;

          // For each entry of the local matrix,
          // find the matching column in the sparse matrix.
          // The entries of the local matrix are not sorted for P2+ elements.
          for(feInt j = 0; j < sizeJ; ++j) {
            for(feInt J = 0; J < numColumns; ++J) {
              if(_mat_ja[debut + J] == adrJ[j]) {
                _mat_values[debut + J] += Ae[niElm[i]][njElm[j]];
                break;
              }
            }
          }
        }
        niElm.clear();
        njElm.clear();
      }
    }
  }

  feInfoCond(FE_VERBOSE > 1, "\t\t\t\tAssembled jacobian matrix in %f s", toc());
  
  viewMatrix();

#if defined(HAVE_MPI)
  }
#endif
}

void feLinearSystemMklPardiso::assembleResiduals(feSolution *sol)
{
  tic();

#if defined(HAVE_MPI)
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Assemble on all procs if system is distributed
  // OR
  // Assemble on master proc only if system is sequential
  if(_ownershipSplit >= 0 || (_ownershipSplit < 0 && rank == 0)) {
    feInfo("Assembling RHS with %d and proc %d", _ownershipSplit, rank);
#else
  int rank = 0;
#endif

  PardisoInt low = _ownedLowerBounds[rank];
  PardisoInt high = _ownedUpperBounds[rank];

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

        // Determine global assignment indices
        adrI = f->getAdrI();
        sizeI = adrI.size();
        for(feInt i = 0; i < sizeI; i++) {
          if(adrI[i] < _nInc) {
            // Remove indices in adrI that are not owned by this proc
            if(low <= adrI[i] && adrI[i] <= high) {
              residu[adrI[i] - low] += Be[i];
            }
          }
        }
      }
    }
  }
  feInfoCond(FE_VERBOSE > 1, "\t\t\t\tAssembled global residual in %f s", toc());

  viewRHS();

#if defined(HAVE_MPI)
  }
#endif
}

void feLinearSystemMklPardiso::assignResidualToDCResidual(feSolutionContainer *solContainer)
{
  for(feInt i = 0; i < _nInc; ++i) solContainer->_fResidual[0][i] = residu[i];
}

void feLinearSystemMklPardiso::applyCorrectionToResidual(double coeff, std::vector<double> &d)
{
  for(int i = 0; i < _nInc; ++i) residu[i] += coeff * d[i];
}

void feLinearSystemMklPardiso::correctSolution(feSolution *sol)
{
  std::vector<double> &solArray = sol->getSolutionReference();

#if defined(HAVE_MPI)
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Correct solution on all procs if system is distributed
  if(_ownershipSplit >= 0) {

    PardisoInt low  = _ownedLowerBounds[rank];
    PardisoInt high = _ownedUpperBounds[rank];

    for(feInt i = 0; i < high - low + 1; i++) {
      _ownedSolution[i] = solArray[low + i] + du[i];
    }

    // Synchronize the solution vector on all processes
    MPI_Allgatherv(_ownedSolution, high-low+1, MPI_DOUBLE, solArray.data(),
      _numOwnedRowsOnAllProc, _ownedLowerBounds, MPI_DOUBLE, MPI_COMM_WORLD);

  }
  // OR
  // Correct solution on master proc only if system is sequential
  else if(rank == 0) {
#endif
    // Increment the solution array
    for(feInt i = 0; i < _nInc; i++) { solArray[i] += du[i]; }
#if defined(HAVE_MPI)
  }
#endif
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

void feLinearSystemMklPardiso::permute()
{
  // Do nothing, at least for now (permutation is done by Pardiso)
}

bool feLinearSystemMklPardiso::solve(double *normDx, double *normResidual, double *normAxb,
                                     int *nIter)
{
  std::string pardisoStep;

  *normAxb = 0.0;
  *nIter = 1;

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

#if defined(HAVE_MPI)
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0)
#endif
  {
    feInfoCond(FE_VERBOSE > 1, "\t\t\t\tSolved linear system with MKL Pardiso in %f s", toc());
  }

  // Keep the data structures from first symbolic factorization
  // for subsequent solves
  symbolicFactorization = false;

  *normDx = vectorMaxNorm(_numOwnedRows, du);
  *normResidual = vectorMaxNorm(_numOwnedRows, residu);

  return true;
}

void feLinearSystemMklPardiso::setToZero()
{
  if(recomputeMatrix) setMatrixToZero();
  setResidualToZero();
}

void feLinearSystemMklPardiso::setMatrixToZero()
{
  // for(feInt i = 0; i < _nnz; i++) _mat_values[i] = 0;
  for(feInt i = 0; i < _numOwnedNNZ; i++) _mat_values[i] = 0;
}

void feLinearSystemMklPardiso::setResidualToZero()
{
  // for(feInt i = 0; i < _nInc; i++) residu[i] = 0;
  for(feInt i = 0; i < _numOwnedRows; i++) residu[i] = 0;
}

void feLinearSystemMklPardiso::assemble(feSolution *sol)
{
  if(recomputeMatrix) {
    this->assembleMatrices(sol);
  }
  this->assembleResiduals(sol);
}

void feLinearSystemMklPardiso::constrainEssentialComponents(feSolution *sol)
{
  tic();
  // Initialize data at first pass
  if(_initializeConstrainData) {

    _initializeConstrainData = false;

    _rowsToConstrain.reserve(_nInc);

    std::vector<feInt> adr;
    for(auto space : sol->_spaces) {
      int nComponents = space->getNumComponents();
      if(nComponents > 1) {
        adr.resize(space->getNumFunctions(), 0);
        for(int i = 0; i < nComponents; ++i) {
          if(space->isEssentialComponent(i)) {
            for(int iElm = 0; iElm < space->getNumElements(); ++iElm) {
              // Constrain matrix and RHS
              space->initializeAddressingVector(iElm, adr);
              for(int j = 0; j < space->getNumFunctions(); ++j) {
                if(j % nComponents == i) {
                  feInt DOF = adr[j];
                  if(adr[j] < _nInc) {

                    #if defined(HAVE_MPI)
                      feErrorMsg(FE_STATUS_ERROR, "constrainEssentialComponents must be modified for a distributed matrix!");
                      finalize();
                      exit(0);
                    #endif
                    // A DOF shared between this space and another
                    // essential space has a tag higher than _nInc,
                    // hence  we cannot constrain it, but it's
                    // already an essential DOF for which there is no
                    // need to solve (its value will be imposed by the
                    // other space though, so check for spaces overlap).
                    _rowsToConstrain.push_back(DOF);
                  }
                }
              }
            }
          }
        }
      }
    }

    _numOccurencesInJa.resize(_rowsToConstrain.size(), 0);
    _posOccurencesInJa.resize(_rowsToConstrain.size());

    #if defined(HAVE_OMP)
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < _rowsToConstrain.size(); ++i)
    {
      feInt row = _rowsToConstrain[i];
      // Find every occurence of 'row' in ja
      for(feInt j = 0; j < _nnz; ++j) {
        if(_mat_ja[j] == row) {
          _numOccurencesInJa[i]++;
          _posOccurencesInJa[i].insert(j);
        }
      }
    }
  }

  // Constrain matrix and RHS
  for(size_t iRow = 0; iRow < _rowsToConstrain.size(); ++iRow)
  {
    feInt row = _rowsToConstrain[iRow];
    // Constrain columns. For every occurence of 'row' in ja,
    // set the associated value to zero
    for(auto jRow : _posOccurencesInJa[iRow]) {
      _mat_values[jRow] = 0.;
    }

    // Then constrain rows. Assign 0 to all columns of the row
    // and set 1 at (row,row)
    int debut = _mat_ia[row];
    int fin = _mat_ia[row + 1];
    for(feInt j = 0; j < fin - debut; ++j) {
      _mat_values[debut + j] = 0.;
      if(_mat_ja[debut + j] == row)
        _mat_values[debut + j] = 1.;
    }

    // Constrain RHS
    residu[row] = 0.;
  }

  feInfoCond(FE_VERBOSE > 1, "\t\t\t\tConstrained essential DOFs in %f s", toc());
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
  PHASE = 11;

#if defined(HAVE_MPI)
  int comm =  MPI_Comm_c2f( MPI_COMM_WORLD );
  // Cluster Pardiso automatically falls back to Pardiso if there is a single MPI process
    cluster_sparse_solver(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &_nInc, _mat_values, _mat_ia, _mat_ja, &IDUM, &NRHS, IPARM, &MSGLVL,
               &DDUM, &DDUM, &comm, &ERROR);
#else
    pardiso(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &_nInc, _mat_values, _mat_ia, _mat_ja, &IDUM, &NRHS, IPARM, &MSGLVL,
               &DDUM, &DDUM, &ERROR);
#endif
}

void feLinearSystemMklPardiso::mklFactorization(void)
{
  PHASE = 22;

#if defined(HAVE_MPI)
  int comm =  MPI_Comm_c2f( MPI_COMM_WORLD );
    cluster_sparse_solver(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &_nInc, _mat_values, _mat_ia, _mat_ja, &IDUM, &NRHS, IPARM, &MSGLVL,
               &DDUM, &DDUM, &comm, &ERROR);
#else
    pardiso(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &_nInc, _mat_values, _mat_ia, _mat_ja, &IDUM, &NRHS, IPARM, &MSGLVL,
               &DDUM, &DDUM, &ERROR);
#endif
}

void feLinearSystemMklPardiso::mklSolve(void)
{
  PHASE = 33;

  // Reset solution vector
  for(feInt i = 0; i < _numOwnedRows; i++) du[i] = 0.0;

#if defined(HAVE_MPI)
  int comm =  MPI_Comm_c2f( MPI_COMM_WORLD );
    cluster_sparse_solver(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &_nInc, _mat_values, _mat_ia, _mat_ja, &IDUM, &NRHS, IPARM, &MSGLVL,
               residu, du, &comm, &ERROR);
#else
    pardiso(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &_nInc, _mat_values, _mat_ia, _mat_ja, &IDUM, &NRHS, IPARM, &MSGLVL,
               residu, du, &ERROR);
#endif
}

feLinearSystemMklPardiso::~feLinearSystemMklPardiso(void)
{
  /* Release internal memory. */
  PHASE = -1;

#if defined(HAVE_MPI)
  int comm =  MPI_Comm_c2f( MPI_COMM_WORLD );
    cluster_sparse_solver(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE,
             &_nInc, &DDUM, _mat_ia, _mat_ja, &IDUM, &NRHS,
             IPARM, &MSGLVL, &DDUM, &DDUM, &comm, &ERROR);
#else
    pardiso(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE,
             &_nInc, &DDUM, _mat_ia, _mat_ja, &IDUM, &NRHS,
             IPARM, &MSGLVL, &DDUM, &DDUM, &ERROR);
#endif

  free(_ownedLowerBounds);
  free(_ownedUpperBounds);
#if defined(HAVE_MPI)
  free(_numOwnedRowsOnAllProc);
  free(_ownedSolution);
#endif

  // delete[] _mat_ia;
  // delete[] _mat_ja;
  // delete[] _mat_values;
  free(_mat_ia);
  free(_mat_ja);
  free(_mat_values);
  delete[] du;
  delete[] residu;
}

void feLinearSystemMklPardiso::setPivot(int pivot) { IPARM[9] = (feInt) pivot; }

#endif
