#include "feLinearSystem.h"

#if defined(HAVE_PETSC)
#include "petscksp.h"
#include "petscdraw.h"
#endif

extern int FE_VERBOSE;

//
// Create and allocate the FE matrix and vectors (Newton residual, RHS, linear system residual (Ax-b))
// Create the KSP (PETSc solver) object.
// Should be removed once the initializeMPI() version is thoroughly tested.
//
void feLinearSystemPETSc::initializeSequential()
{
#if defined(HAVE_PETSC)
  PetscMPIInt size;
  PetscErrorCode ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  if(size != 1) {
    feErrorMsg(FE_STATUS_ERROR, "Cannot run on more than 1 processor :/");
    exit(-1);
  }

  // Allocate vectors and matrix
  ierr = VecCreate(PETSC_COMM_WORLD, &_du);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSetSizes(_du, PETSC_DECIDE, _nInc);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSetFromOptions(_du);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = VecDuplicate(_du, &_rhs);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDuplicate(_du, &_linSysRes);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSet(_du, 1.0);

  ierr = VecDuplicate(_du, &_constrainedDOFValue);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSet(_constrainedDOFValue, 0.0);

  // Determine the nonzero structure
  feEZCompressedRowStorage _EZCRS(_nInc, _formMatrices, _numMatrixForms);

  feInt num_nnz = _EZCRS.getNumNNZ();

  std::vector<PetscInt> nnz(_nInc, 0);
  for(int i = 0; i < _nInc; ++i) {
    nnz[i] = _EZCRS.getNnzAt(i);
  }

  bool withPrealloc = true;
  bool printInfos = true;

  if(withPrealloc) {
    ierr = MatCreate(PETSC_COMM_WORLD, &_A);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, _nInc, _nInc);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetFromOptions(_A);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSeqAIJSetPreallocation(_A, 0, nnz.data());
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetUp(_A);

    if(printInfos) {
      MatInfo info;
      double mal, nz_a, nz_u, mem, nz_un;

      MatGetInfo(_A, MAT_LOCAL, &info);
      mal = info.mallocs;
      nz_a = info.nz_allocated;
      nz_u = info.nz_used;
      nz_un = info.nz_unneeded;
      mem = info.memory;

      feInfo("\t\tAdditional info from PETSc for matrix creation:");
      feInfo("\t\tNumber of mallocs during MatSetValues() = %f", mal);
      feInfo("\t\tMemory allocated                        = %f", mem );
      feInfo("\t\tNumber of nonzero allocated             = %f", nz_a);
      feInfo("\t\tNumber of nonzero used                  = %f", nz_u);
      feInfo("\t\tNumber of nonzero unneeded              = %f", nz_un);
      feInfo("");
    }
  } else {
    // Without allocation (bad) :
    ierr = MatCreate(PETSC_COMM_WORLD, &_A);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, _nInc, _nInc);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetFromOptions(_A);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetUp(_A);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  }

  // Insert elements on the diagonal and set to zero to
  // avoid using a weak form for a block of zeros
  PetscInt indiceDiag;
  PetscScalar val = 0.;
  for(int i = 0; i < _nInc; ++i) {
    indiceDiag = i;
    ierr = MatSetValues(_A, 1, &indiceDiag, 1, &indiceDiag, &val, INSERT_VALUES);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  }
  ierr = MatAssemblyBegin(_A, MAT_FLUSH_ASSEMBLY);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatAssemblyEnd(_A, MAT_FLUSH_ASSEMBLY);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // Trivial row ownership arrays for a single process, used in correctSolution()
  _ownedLowerBounds = (int*) malloc(1*sizeof(int));
  _ownedUpperBounds = (int*) malloc(1*sizeof(int));
  _numOwnedRows = (int*) malloc(1*sizeof(int));
  _ownedLowerBounds[0] = 0;
  _ownedUpperBounds[0] = _nInc;
  _numOwnedRows[0] = _nInc;

  // Initialize trivial solution buffer owned on this process (= full solution)
  _ownedSolution = (double*) malloc(_nInc * sizeof(double));

  // Create the Krylov solver (default is GMRES)
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = KSPSetType(ksp, KSPGMRES);
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = KSPSetOperators(ksp, _A, _A);
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // Set preconditioner
  ierr = KSPGetPC(ksp, &preconditioner);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = PCSetType(preconditioner, PCILU);
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = KSPSetTolerances(ksp, (PetscReal)_rel_tol, (PetscReal)_abs_tol, (PetscReal)_div_tol,
                          (PetscInt)_max_iter);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetFromOptions(ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  PetscInt M, N, Nres, Ndx, NlinSysRes;
  MatGetSize(_A, &M, &N);
  VecGetSize(_rhs, &Nres);
  VecGetSize(_du, &Ndx);
  VecGetSize(_linSysRes, &NlinSysRes);

  // ierr = MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // // Initialize diagonal to zero
  // // Vec diag;
  // // VecDuplicate(_du, &diag);
  // // MatGetDiagonal(_A, diag);
  // // VecSet(diag, 0.);
  // std::vector<PetscInt> d(_nInc, 0);
  // std::vector<PetscScalar> dval(_nInc, 0.);
  // for(int i = 0; i < _nInc; ++i)
  //   d[i] = i;
  // ierr = MatSetValues(_A, d.size(), d.data(), d.size(), d.data(), dval.data(), INSERT_VALUES);
  // // ierr = MatDiagonalSet(_A, diag, INSERT_VALUES);
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);

  feInfoCond(FE_VERBOSE > 0, "\t\tCreated a PETSc linear system of size %d x %d", M, N);
  feInfoCond(FE_VERBOSE > 0, "\t\tNumber of nonzero elements: %d", num_nnz);
  MatType   matType;
  ierr = MatGetType(_A, &matType);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  feInfoCond(FE_VERBOSE > 0, "\t\tPETSc matrix type: %s", matType);
  feInfoCond(FE_VERBOSE > 0, "\t\tIterative linear solver info:");

  KSPType krylovType;
  ierr = KSPGetType(ksp, &krylovType);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tKrylov method: %s", krylovType);

  PCType preconditionerType;
  ierr = PCGetType(preconditioner, &preconditionerType);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tPreconditioner: %s", preconditionerType);

  feInfoCond(FE_VERBOSE > 0, "\t\t\tRelative tolerance: %1.4e", _rel_tol);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tAbsolute tolerance: %1.4e", _abs_tol);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tDivergence criterion: %1.4e", _div_tol);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tMax number of iterations: %d", _max_iter);

#endif
}

int feLinearSystemPETSc::initializeMPI()
{
#if defined(HAVE_PETSC)

  feInt num_nnz = -1;

  PetscMPIInt size, rank;
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

  //
  // Create matrix
  //
  PetscCall(MatCreate(PETSC_COMM_WORLD, &_A));
  PetscCall(MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, _nInc, _nInc));
  PetscCall(MatSetFromOptions(_A));

  //
  // Determine sparsity pattern and preallocate matrix
  //
  // Get and share the upper bound of the owned (global) indices on each process
  PetscInt low, high;
  PetscCall(MatGetOwnershipRange(_A, &low, &high));

  int localLow[1] = {low}, localHigh[1] = {high};
  _ownedUpperBounds = (int*) malloc(size*sizeof(int));
  _ownedLowerBounds = (int*) malloc(size*sizeof(int)); 
  PetscCallMPI(MPI_Allgather(localLow , 1, MPI_INT, _ownedLowerBounds, 1, MPI_INT, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Allgather(localHigh, 1, MPI_INT, _ownedUpperBounds, 1, MPI_INT, PETSC_COMM_WORLD));

  _numOwnedRows = (int*) malloc(size*sizeof(int));
  for(int i = 0; i < size; ++i) {
    _numOwnedRows[i] = _ownedUpperBounds[i] - _ownedLowerBounds[i];
  }

  // Initialize solution buffer owned on this process
  _ownedSolution = (double*) malloc(_numOwnedRows[rank] * sizeof(double));

  PetscInt numRowsOnProc = high - low;
  std::vector<feInt> diagNNZOnAllProcs, offdiagNNZOnAllProcs;
  int *diagNNZ = (int *) malloc(numRowsOnProc*sizeof(int));
  int *offdiagNNZ = (int *) malloc(numRowsOnProc*sizeof(int));

  // For sequential preallocation if running on single process?
  std::vector<PetscInt> nnz(_nInc, 0);

  if(rank == 0) {

    // Determine the nonzero structure
    feEZCompressedRowStorage _EZCRS(_nInc, _formMatrices, _numMatrixForms, _ownedUpperBounds);
    num_nnz = _EZCRS.getNumNNZ();

    for(int i = 0; i < _nInc; ++i) {
      nnz[i] = _EZCRS.getNnzAt(i);
    }

    diagNNZOnAllProcs = _EZCRS.getDiagNNZForAllMPIProcs();
    offdiagNNZOnAllProcs = _EZCRS.getOffDiagNNZForAllMPIProcs();
  }

  // Get the nonzero structure on each proc
  PetscCallMPI(MPI_Scatterv(diagNNZOnAllProcs.data(), _numOwnedRows, _ownedLowerBounds, MPI_INT, diagNNZ, numRowsOnProc, MPI_INT, 0, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Scatterv(offdiagNNZOnAllProcs.data(), _numOwnedRows, _ownedLowerBounds, MPI_INT, offdiagNNZ, numRowsOnProc, MPI_INT, 0, PETSC_COMM_WORLD));

  // Allocate parallel matrix, or sequential matrix if running on a single process
  // Only one of these preallocations will do something
  PetscCall(MatMPIAIJSetPreallocation(_A, 0, diagNNZ, 0, offdiagNNZ));
  PetscCall(MatSeqAIJSetPreallocation(_A, 0, nnz.data()));
  
  bool printInfos = true;

  if(printInfos) {
    int ranktoprint = 0;
    while(ranktoprint < size) {
      if(rank == ranktoprint) {
        MatInfo info;
        double mal, nz_a, nz_u, mem, nz_un;

        PetscCall(MatGetInfo(_A, MAT_LOCAL, &info));
        mal = info.mallocs;
        nz_a = info.nz_allocated;
        nz_u = info.nz_used;
        nz_un = info.nz_unneeded;
        mem = info.memory;

        feInfo("\t\tAdditional info from PETSc for matrix creation:");
        feInfo("\t\t\tNumber of mallocs during MatSetValues() on proc %d : %f", rank, mal);
        feInfo("\t\t\tMemory allocated                        on proc %d : %f", rank, mem );
        feInfo("\t\t\tNumber of nonzero allocated             on proc %d : %f", rank, nz_a);
        feInfo("\t\t\tNumber of nonzero used                  on proc %d : %f", rank, nz_u);
        feInfo("\t\t\tNumber of nonzero unneeded              on proc %d : %f", rank, nz_un);
      }

      ranktoprint++;
      MPI_Barrier(PETSC_COMM_WORLD);
    }
  }

  // Insert elements on the diagonal and set to zero to
  // avoid using a weak form for a block of zeros
  PetscInt indiceDiag;
  PetscScalar val = 0.;
  for(int i = 0; i < _nInc; ++i) {
    indiceDiag = i;
    PetscCall(MatSetValues(_A, 1, &indiceDiag, 1, &indiceDiag, &val, INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(_A, MAT_FLUSH_ASSEMBLY));
  PetscCall(MatAssemblyEnd(_A, MAT_FLUSH_ASSEMBLY));

  //
  // Create vectors
  //
  // Create solution vector du
  PetscCall(VecCreate(PETSC_COMM_WORLD, &_du));
  PetscCall(VecSetSizes(_du, PETSC_DECIDE, _nInc));
  PetscCall(VecSetFromOptions(_du));

  // Duplicate into RHS and residual
  PetscCall(VecDuplicate(_du, &_rhs));
  PetscCall(VecDuplicate(_du, &_linSysRes));
  PetscCall(VecSet(_du, 1.0));

  PetscCall(VecDuplicate(_du, &_constrainedDOFValue));
  PetscCall(VecSet(_constrainedDOFValue, 0.0));

  //
  // Create the Krylov solver (default is GMRES)
  //
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  // PetscCall(KSPSetType(ksp, KSPGMRES));
  // PetscCall(KSPSetOperators(ksp, _A, _A));

  // Set preconditioner
  PetscCall(KSPGetPC(ksp, &preconditioner));
  // PetscCall(PCSetType(preconditioner, PCILU));

  PetscCall(KSPSetTolerances(ksp, (PetscReal)_rel_tol, (PetscReal)_abs_tol, (PetscReal)_div_tol,
                          (PetscInt)_max_iter));
  // PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));

  // Set runtime options, overriding the options above
  PetscCall(KSPSetFromOptions(ksp));

  PetscInt M, N, Nres, Ndx, NlinSysRes;
  PetscCall(MatGetSize(_A, &M, &N));
  PetscCall(VecGetSize(_rhs, &Nres));
  PetscCall(VecGetSize(_du, &Ndx));
  PetscCall(VecGetSize(_linSysRes, &NlinSysRes));

  if(rank == 0) {
    feInfoCond(FE_VERBOSE > 0, "");
    feInfoCond(FE_VERBOSE > 0, "\t\tCreated a PETSc linear system of size %d x %d", M, N);
    feInfoCond(FE_VERBOSE > 0, "\t\tNumber of nonzero elements: %d", num_nnz);
    MatType   matType;
    PetscCall(MatGetType(_A, &matType));
    feInfoCond(FE_VERBOSE > 0, "\t\tPETSc matrix type: %s", matType);
    feInfoCond(FE_VERBOSE > 0, "\t\tIterative linear solver info:");

    KSPType ksptype;
    PetscCall(KSPGetType(ksp, &ksptype));
    feInfoCond(FE_VERBOSE > 0, "\t\t\tKrylov method: %s", ksptype);

    PCType pctype;
    PetscCall(PCGetType(preconditioner, &pctype));
    feInfoCond(FE_VERBOSE > 0, "\t\t\tPreconditioner: %s", pctype);

    feInfoCond(FE_VERBOSE > 0, "\t\t\tRelative tolerance: %1.4e", _rel_tol);
    feInfoCond(FE_VERBOSE > 0, "\t\t\tAbsolute tolerance: %1.4e", _abs_tol);
    feInfoCond(FE_VERBOSE > 0, "\t\t\tDivergence criterion: %1.4e", _div_tol);
    feInfoCond(FE_VERBOSE > 0, "\t\t\tMax number of iterations: %d", _max_iter);
  }
#endif

  return 0;
}

void feLinearSystemPETSc::initialize()
{
  //   this->initializeSequential();

  // This allocates a sequential matrix if #MPI proc = 1
  // or if given the runtime option -mpi_linear_solver_server
  this->initializeMPI();
}

feLinearSystemPETSc::feLinearSystemPETSc(int argc, char **argv,
                                         std::vector<feBilinearForm *> bilinearForms,
                                         int numUnknowns)
  : feLinearSystem(bilinearForms), _argc(argc), _argv(argv)
#if defined(HAVE_PETSC)
    , _nInc(numUnknowns)
#endif
{
  this->initialize();
}

void feLinearSystemPETSc::viewMatrix()
{
#if defined(HAVE_PETSC)
  if(_displayMatrixInConsole) {
    MatView(_A, PETSC_VIEWER_STDOUT_WORLD);
  }

  if(_displayMatrixInWindow) {
    PetscErrorCode ierr = 0;
    // ierr = MatGetOrdering(_A, MATORDERINGRCM, &_rowMap, &_colMap); CHKERRV(ierr);
    // ierr = MatPermute(_A, _rowMap, _colMap, &_Ap); CHKERRV(ierr);

    PetscViewer viewer;
    PetscDraw draw;
    ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1200, 1200, &viewer);
    CHKERRV(ierr);
    ierr = MatView(_A, viewer);
    CHKERRV(ierr);
    ierr = PetscViewerDrawGetDraw(viewer, 0, &draw);
    CHKERRV(ierr);
    ierr = PetscDrawSetPause(draw, -1);
    CHKERRV(ierr); // Wait for user
    PetscDrawPause(draw);
    ierr = PetscViewerDestroy(&viewer);
    CHKERRV(ierr);
  }

#endif
}

void feLinearSystemPETSc::viewRHS()
{
#if defined(HAVE_PETSC)
  if(_displayRHSInConsole) {
    VecView(_rhs, PETSC_VIEWER_STDOUT_WORLD);
  }
#endif
}

void feLinearSystemPETSc::viewResidual()
{
#if defined(HAVE_PETSC)
  VecView(_du, PETSC_VIEWER_STDOUT_WORLD);
#endif
}

void feLinearSystemPETSc::getRHSMaxNorm(double *norm)
{
#if defined(HAVE_PETSC)
  PetscErrorCode ierr = VecNorm(_rhs, NORM_MAX, norm);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif
}

void feLinearSystemPETSc::getResidualMaxNorm(double *norm)
{
#if defined(HAVE_PETSC)
  PetscErrorCode ierr = VecNorm(_du, NORM_MAX, norm);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif
}

void feLinearSystemPETSc::setToZero()
{
  if(recomputeMatrix) {
    this->setMatrixToZero();
  }
  this->setResidualToZero();
}

void feLinearSystemPETSc::setMatrixToZero()
{
#if defined(HAVE_PETSC)
  PetscErrorCode ierr = MatZeroEntries(_A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif
}

void feLinearSystemPETSc::setResidualToZero()
{
#if defined(HAVE_PETSC)
  PetscErrorCode ierr = VecZeroEntries(_rhs);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif
}

void feLinearSystemPETSc::assembleMatrices(feSolution *sol)
{
#if defined(HAVE_PETSC)

  PetscErrorCode ierr = 0;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  PetscInt low, high;
  ierr = MatGetOwnershipRange(_A, &low, &high);

  tic();
  for(feInt eq = 0; eq < _numMatrixForms; ++eq) {
    feBilinearForm *f = _formMatrices[eq];
    feCncGeo *cnc = f->getCncGeo();
    int numColors = cnc->getNbColor();
    const std::vector<int> &numElemPerColor = cnc->getNbElmPerColor();
    const std::vector<std::vector<int> > &listElmPerColor = cnc->getListElmPerColor();
    int numElementsInColor;
    std::vector<int> listElmC;

    for(int iColor = 0; iColor < numColors; ++iColor) {
      numElementsInColor = numElemPerColor[iColor];
      listElmC = listElmPerColor[iColor];

      int elm;
      double **Ae;
      feInt sizeI;
      feInt sizeJ;
      std::vector<feInt> niElm;
      std::vector<feInt> njElm;
      std::vector<feInt> adrI;
      std::vector<feInt> adrJ;
      std::vector<PetscScalar> values;

#if defined(HAVE_OMP)
#pragma omp parallel for private(elm, f, niElm, njElm, sizeI, sizeJ, adrI, adrJ, Ae, ierr, values)
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

        // //////////////////////////////////////////////////////////////////
        // MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        // if(iColor == 0 && iElm == 4) {
        //   feInfo("sizes = %d - %d", adrI.size(), adrJ.size());
        //   for(size_t i = 0; i < adrI.size(); ++i) {
        //     feInfo("From P%d - A adri = %d - adrj = %d", rank, adrI[i], adrJ[i]);
        //   }
        // }
        // MPI_Barrier(PETSC_COMM_WORLD);

        // In addition remove row indices in adrI that are not owned by this proc
        
        auto it_nielm = niElm.begin();
        for(auto it_i = adrI.begin();  it_i != adrI.end(); )
        {
          if(low <= (*it_i) && (*it_i) < high) {
            // In owned range of rows : do nothing
            ++it_i;
            ++it_nielm;
          } else {
            // Erase in adrI only
            it_i = adrI.erase(it_i);
            it_nielm = niElm.erase(it_nielm);
          }
        }

        // if(iColor == 0 && iElm == 4) {
        //   int ranktoprint = 0;
        //   while(ranktoprint < size) {
        //     if(rank == ranktoprint) {
        //       feInfo("sizes after erase on P%d = %d - %d", rank, adrI.size(), adrJ.size());
        //       for(size_t i = 0; i < adrI.size(); ++i) {
        //         for(size_t j = 0; j < adrJ.size(); ++j) {
        //           feInfo("proc %d - range %d - %d - assigning at (%d, %d)", rank, Istart, Iend, adrI[i], adrJ[j]);
        //         }
        //       }
        //     }

        //     ranktoprint++;
        //     MPI_Barrier(PETSC_COMM_WORLD);
        //   }
        // }
        // //////////////////////////////////////////////////////////////////

        // Flatten Ae at relevant indices
        sizeI = adrI.size();
        sizeJ = adrJ.size();

        if(!(adrI.size() == niElm.size())){ feErrorMsg(FE_STATUS_ERROR, "assert failed : adrI size = %d - niElm size = %d", adrI.size(), niElm.size()); };
        if(!(adrJ.size() == njElm.size())){ feErrorMsg(FE_STATUS_ERROR, "assert failed : adrJ size = %d - njElm size = %d", adrJ.size(), njElm.size()); };

        values.resize(sizeI * sizeJ);
        for(feInt i = 0; i < sizeI; ++i) {
          for(feInt j = 0; j < sizeJ; ++j) {
            values[sizeJ * i + j] = Ae[niElm[i]][njElm[j]];
          }
        }

// Increment global matrix
// Seems to be issues without the critical, even though it should not race 
#if defined(HAVE_OMP)
#pragma omp critical
#endif
        ierr = MatSetValues(_A, adrI.size(), adrI.data(), adrJ.size(), adrJ.data(), values.data(), ADD_VALUES);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
        niElm.clear();
        njElm.clear();
      }
    }
  }

  ierr = MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);   CHKERRABORT(PETSC_COMM_WORLD, ierr);

  feInfoCond(FE_VERBOSE > 1, "\t\t\t\tAssembled jacobian matrix in %f s", toc());

  // double normMat;
  // ierr = MatNorm(_A, NORM_FROBENIUS, &normMat);
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);

  viewMatrix();

  if(_exportMatrixMatlab) {
    feInfoCond(FE_VERBOSE > 1, "\t\t\tExporting global matrix to Matlab");
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrix.m", &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    MatView(_A, viewer);
    PetscViewerPopFormat(viewer);
    PetscViewerDestroy(&viewer);
  }

#endif
}

void feLinearSystemPETSc::assembleResiduals(feSolution *sol)
{
#if defined(HAVE_PETSC)

  PetscErrorCode ierr = 0;

  PetscInt low, high;
  ierr = VecGetOwnershipRange(_rhs, &low, &high);

  tic();
  for(feInt eq = 0; eq < _numResidualForms; ++eq) {
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
      std::vector<PetscScalar> values;

#if defined(HAVE_OMP)
#pragma omp parallel for private(elm, niElm, Be, f, adrI, ierr, values, sizeI)
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
        niElm.reserve(sizeI);
        for(feInt i = 0; i < sizeI; ++i) {
          if(adrI[i] < _nInc) niElm.push_back(i);
        }

        adrI.erase(std::remove_if(adrI.begin(), adrI.end(),
                                  [this](const int &x) { return x >= this->_nInc; }),
                   adrI.end());

        // //////////////////////////////////////////////////////////////////
        // Remove indices in adrI that are not owned by this proc
        auto it_nielm = niElm.begin();
        for(auto it_i = adrI.begin();  it_i != adrI.end(); )
        {
          if(low <= (*it_i) && (*it_i) < high) {
            // In owned range of rows : do nothing
            ++it_i;
            ++it_nielm;
          } else {
            // Erase in adrI only
            it_i = adrI.erase(it_i);
            it_nielm = niElm.erase(it_nielm);
          }
        }
        // //////////////////////////////////////////////////////////////////

        sizeI = adrI.size();
        values.resize(sizeI);
        for(feInt i = 0; i < sizeI; ++i) {
          values[i] = Be[niElm[i]];
        }

// Increment global residual
#if defined(HAVE_OMP)
#pragma omp critical
#endif
        ierr = VecSetValues(_rhs, adrI.size(), adrI.data(), values.data(), ADD_VALUES);
        niElm.clear();
      }
    }
  }

  ierr = VecAssemblyBegin(_rhs);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecAssemblyEnd(_rhs);    CHKERRABORT(PETSC_COMM_WORLD, ierr);

  feInfoCond(FE_VERBOSE > 1, "\t\t\t\tAssembled global residual in %f s", toc());

  viewRHS();

#endif
}

void feLinearSystemPETSc::assemble(feSolution *sol)
{
  if(recomputeMatrix) {
    this->assembleMatrices(sol);
  }
  this->assembleResiduals(sol);
}

void feLinearSystemPETSc::constrainEssentialComponents(feSolution *sol)
{
#if defined(HAVE_PETSC)
  tic();
  std::vector<PetscInt> rowsToConstrain;
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
                PetscInt DOF = adr[j];
                if(adr[j] < _nInc) {
                  // A DOF shared between this space and another
                  // essential space has a tag higher than _nInc,
                  // hence  we cannot constrain it, but it's
                  // already an essential DOF for which there is no
                  // need to solve (its value will be imposed by the
                  // other space though, so check for spaces overlap).
                  rowsToConstrain.push_back(DOF);
                }
              }
            }
          }
        }
      }
    }
  }

  PetscErrorCode ierr = MatZeroRowsColumns(_A, rowsToConstrain.size(), rowsToConstrain.data(), 1.,
                                           _constrainedDOFValue, _rhs);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  feInfoCond(FE_VERBOSE > 1, "\t\t\t\tConstrained essential DOFs in %f s", toc());
#endif
}

//
// Apply RCMK permutation
//
void feLinearSystemPETSc::permute()
{
#if defined(HAVE_PETSC)

  PetscErrorCode ierr = 0;

  ierr = MatGetOrdering(_A, MATORDERINGRCM, &_rowMap, &_colMap); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatPermute(_A, _rowMap, _colMap, &_Ap); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecPermute(_rhs, _rowMap, PETSC_FALSE); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatDestroy(&_A); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // _A = _Ap; /* Replace original operator with permuted version */

  // Write matrix and RHS to binary file
  // const char* fileName = "matrixBinary.txt";
  // PetscViewer viewer;
  // ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName, FILE_MODE_WRITE, &viewer); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatView(_A, viewer); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = VecView(_rhs, viewer); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = PetscViewerDestroy(&viewer); CHKERRABORT(PETSC_COMM_WORLD, ierr);

  if(_displayMatrixInWindow) {
    PetscViewer viewer;
    PetscDraw draw;
    ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1200, 1200, &viewer);
    CHKERRV(ierr);
    ierr = MatView(_Ap, viewer);
    CHKERRV(ierr);
    ierr = PetscViewerDrawGetDraw(viewer, 0, &draw);
    CHKERRV(ierr);
    ierr = PetscDrawSetPause(draw, -1);
    CHKERRV(ierr); // Wait for user
    PetscDrawPause(draw);
    ierr = PetscViewerDestroy(&viewer);
    CHKERRV(ierr);
  }
#endif
}

// Solve the linear system and compute max norms of solution and residuals
bool feLinearSystemPETSc::solve(double *normSolution, double *normRHS, double *normResidualAxMinusb,
                                int *nIter)
{
#if defined(HAVE_PETSC)

  PetscErrorCode ierr = 0;

  // PetscCall(KSPSetOptionsPrefix(ksp, "prefix_test_"));

  if(_permute) {
    // Use permuted matrix as operator
    ierr = KSPSetOperators(ksp, _Ap, _Ap);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    // MatSetOptionsPrefix(_Ap, "prefix_test_");
  } else {
    // Use FE matrix as operator
    ierr = KSPSetOperators(ksp, _A, _A);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    // MatSetOptionsPrefix(_A, "prefix_test_");
  }

  // Set tolerances, in case the have changed since creation
  ierr = KSPSetTolerances(ksp, (PetscReal) _rel_tol, (PetscReal) _abs_tol, (PetscReal) _div_tol, (PetscInt) _max_iter);

  // Override with command line options
  ierr = KSPSetFromOptions(ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = VecNorm(_rhs, NORM_MAX, normRHS);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // Solve
  tic();
  ierr = KSPSolve(ksp, _rhs, _du);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  KSPGetIterationNumber(ksp, nIter);
  // PetscPrintf(PETSC_COMM_WORLD, "Iterations %" PetscInt_FMT "\n", *nIter);
  feInfoCond(FE_VERBOSE > 1, "\t\t\t\tSolved linear system with PETSc in %f s", toc());

  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(ksp, &reason);

  if(reason < 0) {
    // Solve failed
    const char *reasonString;
    KSPGetConvergedReasonString(ksp, &reasonString);
    feWarning("PETSc solve failed with KSP error code: %s", reasonString);
    if(reason == KSP_DIVERGED_PC_FAILED) {
      PCFailedReason reasonpc;
      ierr = PCGetFailedReason(preconditioner, &reasonpc);
      switch(reasonpc) {
        case PC_SETUP_ERROR:
          feWarning("\nPETSc solve failed with preconditioner error code: %s", "PC_SETUP_ERROR");
          break;
        case PC_NOERROR:
          feWarning("\nPETSc solve failed with preconditioner error code: %s", "PC_NOERROR");
          break;
        case PC_FACTOR_STRUCT_ZEROPIVOT:
          feWarning("\nPETSc solve failed with preconditioner error code: %s",
                    "PC_FACTOR_STRUCT_ZEROPIVOT");
          break;
        case PC_FACTOR_NUMERIC_ZEROPIVOT:
          feWarning("\nPETSc solve failed with preconditioner error code: %s",
                    "PC_FACTOR_NUMERIC_ZEROPIVOT");
          break;
        case PC_FACTOR_OUTMEMORY:
          feWarning("\nPETSc solve failed with preconditioner error code: %s",
                    "PC_FACTOR_OUTMEMORY");
          break;
        case PC_FACTOR_OTHER:
          feWarning("\nPETSc solve failed with preconditioner error code: %s", "PC_FACTOR_OTHER");
          break;
#if PETSC_VERSION_GE(3, 16, 0)
        case PC_INCONSISTENT_RHS:
          feWarning("\nPETSc solve failed with preconditioner error code: %s",
                    "PC_INCONSISTENT_RHS");
          break;
#endif
        case PC_SUBPC_ERROR:
          feWarning("\nPETSc solve failed with preconditioner error code: %s", "PC_SUBPC_ERROR");
          break;
        default:
          feWarning("\nPETSc solve failed with unknown preconditioner error code");
          break;
      }
    }
  }

  VecSet(_linSysRes, 0.0);
  if(_permute) {
    MatMult(_Ap, _du, _linSysRes);
    ierr = VecPermute(_du, _colMap, PETSC_TRUE);
    // ierr = VecPermute(_du, _rowMap, PETSC_TRUE);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  } else {
    MatMult(_A, _du, _linSysRes);
  }
  VecAXPY(_linSysRes, -1.0, _rhs);
  ierr = VecNorm(_rhs, NORM_MAX, normRHS);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecNorm(_du, NORM_MAX, normSolution);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecNorm(_linSysRes, NORM_MAX, normResidualAxMinusb);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  if(reason == KSP_DIVERGED_ITS) {
    feWarning("Max number of iteration reached");
    return true;
  }
  if(reason == KSP_DIVERGED_BREAKDOWN) {
    feWarning("Diverged breakdown");
    return true;
  }

  return reason > 0;
#endif
  return false;
}

void feLinearSystemPETSc::correctSolution(feSolution *sol)
{
#if defined(HAVE_PETSC)

  std::vector<double> &solArray = sol->getSolutionReference();

  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscInt low  = _ownedLowerBounds[rank];
  PetscInt high = _ownedUpperBounds[rank];

  const PetscScalar *array;
  VecGetArrayRead(_du, &array);
  // for(int i = 0; i < _nInc; ++i) solArray[i] += array[i];
  for(int i = 0; i < high - low; ++i) {
    _ownedSolution[i] = solArray[low + i] + array[i];
  }
  VecRestoreArrayRead(_du, &array);

  // Synchronize the solution vector on all processes
  MPI_Allgatherv(_ownedSolution, high-low, MPI_DOUBLE, solArray.data(),
    _numOwnedRows, _ownedLowerBounds, MPI_DOUBLE, PETSC_COMM_WORLD);

#endif
}

void feLinearSystemPETSc::assignResidualToDCResidual(feSolutionContainer *solContainer)
{
#if defined(HAVE_PETSC)
  PetscScalar *array;
  VecGetArray(_rhs, &array);
  for(int i = 0; i < _nInc; ++i) solContainer->_fResidual[0][i] = array[i];
  VecRestoreArray(_rhs, &array);
#endif
}

void feLinearSystemPETSc::applyCorrectionToResidual(double coeff, std::vector<double> &d)
{
#if defined(HAVE_PETSC)
  PetscScalar *array;
  VecGetArray(_rhs, &array);
  for(int i = 0; i < _nInc; ++i) array[i] += coeff * d[i];
  VecRestoreArray(_rhs, &array);
#endif
}

feLinearSystemPETSc::~feLinearSystemPETSc()
{
#if defined(HAVE_PETSC)
  free(_ownedLowerBounds);
  free(_ownedUpperBounds);
  free(_numOwnedRows);
  free(_ownedSolution);
  PetscErrorCode ierr = KSPDestroy(&ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatDestroy(&_A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&_linSysRes);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&_rhs);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&_du);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&_constrainedDOFValue);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif
}
